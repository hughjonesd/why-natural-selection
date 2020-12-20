
suppressPackageStartupMessages({
  library(drake)
  library(rlang)
  library(purrr)
  library(tidyr)
  library(readr)
  loadNamespace("future")
})

source("make-data.R") 
# infrastructure available at github.com/hughjonesd-private/import-ukbb-data
source("~/import-ukbb-data/import-ukbb-data.R")
source("make-census-weights.R")
source("regressions.R")
source("weight-data.R")


nomis_file         <- file.path(data_dir, "nomis-2011-statistics.csv")
ghs_file           <- file.path(data_dir, "UKDA-5804-stata8", "stata8", 
                        "Ghs06client.dta")
DC1108EW_file      <- file.path(data_dir, "DC1108EW.csv")
DC1108EW_msoa_file <- file.path(data_dir, "DC1108EW-msoas.csv")
DC5202SC_file      <- file.path(data_dir, "DC5202SC.csv")
msoa_shapefile     <- file.path(data_dir, "infuse_msoa_lyr_2011_clipped", 
                        "infuse_msoa_lyr_2011_clipped.shp")
dep_data_dir       <- file.path(data_dir, "deprivation-data")

weighting_scheme_syms <- rlang::syms(c("flb_weights", "age_qual_weights", 
                        "msoa_weights"))
period_weight_syms <- rlang::syms((c("age_qual_weights", "parent_weights")))

plan <- drake_plan(
  
  score_names  = import_score_names(pgs_dir),
  
  relatedness = target(
                  make_relatedness(file_in(!! relatedness_file)),
                  format = "fst_tbl"
                ),
  
  sib_groups = target(make_sib_groups(relatedness), format = "fst_tbl"),
  
  famhist_raw  = target(
                   import_famhist(famhist_files, pgs_dir),
                   format = "fst_tbl"
                 ),
  
  pcs          = import_pcs(file_in(!! pcs_file)),
  
  ashe_income  = import_ashe_income(file_in(!! ashe_income_file)),
  
  famhist_no_resid      =  target({
                      fnr <- clean_famhist(famhist_raw, score_names, sib_groups)

                      fnr$kids_ss <- fnr$age_at_recruitment >= 45
                      fnr
                    },
                    format = "fst_tbl"
                  ),
  
  famhist = target({
                     famhist <- join_famhist_resid_scores(famhist_no_resid, 
                                                            resid_scores)
                     # use residualized scores as defaults:
                     names(famhist)[names(famhist) %in% score_names] %<>%
                       paste0("_raw")
                     names(famhist) <- gsub("_resid$", "", names(famhist))
                     famhist
                   },
                    format = "fst_tbl"
                  ),

  famhist_msoa = target(
                   find_containing_msoas(famhist, msoa_shapefile), 
                   format = "fst_tbl"
                 ),
  
  famhist_townsend_71 = target(
                          add_deprivation_data(famhist, dep_data_dir),
                          format = "fst_tbl"
                        ),
  
  # mlogit now uses dfidx to create a "long" dataset, but mnlogit 
  # (used below) doesn't update. Skipping for now.
  # fhl_mlogit   =  make_famhist_long_mlogit(famhist, score_names),
  
  resid_scores = target({
                     resid_scores <- compute_resid_scores(famhist_no_resid, 
                                                            pcs, score_names)
                     resid_scores <- subset_resid_scores(resid_scores, 
                                                           famhist_no_resid, 
                                                           score_names)
                     resid_scores
                   },
                   format = "fst_tbl"
                 ),

  ghs_subset   = target(
                   make_ghs_subset(file_in(!! ghs_file)), 
                   format = "fst_tbl"
                 ),
  
  pgs_over_time = calc_pgs_over_time(famhist, score_names),
  
  census_age_qual =  make_census_age_qual(
                       file_in(!! DC1108EW_file),
                       file_in(!! DC5202SC_file)
                     ),
  
  census_msoa = target(
                  make_census_msoa(file_in(!! DC1108EW_msoa_file)), 
                  format = "fst_tbl"
                ),
  
  flb_weights = {
    weight_by_ghs(
            ~ age_at_recruitment + factor(flb_cat) + factor(edu_qual), 
            ghs_subset, 
            famhist
          )
  },
  
  age_qual_weights = weight_by_census_age_qual(famhist, census_age_qual),
  
  parent_weights = weight_parents(famhist, input_weights = age_qual_weights),
  
  msoa_weights = weight_by_census_msoa(famhist, census_msoa, famhist_msoa),
  
  mf_pairs =  target({
                mf_pairs <- make_mf_pairs(file_in(!! mf_pairs_file), famhist, 
                                            resid_scores = NULL, ashe_income)
                mf_pairs <- filter_mf_pairs(mf_pairs)
                mf_pairs
              },
                format = "fst_tbl"
              ), 
  
  rgs = make_rgs(file_in(!! rgs_file)),
  
  res_all = {
    res_sibs <- run_regs_basic(
            dep_var     = "n_sibs", 
            score_names = score_names,
            famhist     = famhist
          )
    res_chn <- run_regs_basic(
            dep_var     = "n_children", 
            score_names = score_names, 
            famhist     = famhist,
            subset      = quote(kids_ss)
          )
    dplyr::bind_rows(
            "N siblings" = res_sibs, 
            "N children" = res_chn, 
            .id = "dep.var"
          )
  },
  
  res_quadratic = {
    famhist %<>% inner_join(age_qual_weights, by = "f.eid")
    map_dfr(score_names, 
              ~run_regs_fml(
                             "n_children ~ {score_name} + I({score_name}^2)", 
                             score_name = .x,
                             famhist    = famhist,
                             weights    = quote(weights)
                           ),
              .id = "score_name"
            )
  },
  
  res_pcs = {
    pc_names <- grep("PC", names(pcs), value = TRUE)
    res_sibs_pcs <- run_regs_pcs( 
      dep_var     = "n_sibs", 
      famhist     = famhist,
      pcs         = pcs
    )
    res_chn_pcs <- run_regs_pcs(
      dep_var     = "n_children", 
      famhist     = famhist,
      pcs         = pcs,
      subset      = quote(kids_ss)
    )
    dplyr::bind_rows(
      "N siblings" = res_sibs_pcs, 
      "N children" = res_chn_pcs, 
      .id = "dep.var"
    ) 
  },
  
  res_wt =  target(
              map_dfr(score_names, 
                run_regs_weighted,
                famhist     = famhist,
                weight_data = weighting_scheme,
                dep.var     = "n_children",
                subset      = quote(kids_ss)
              ),
              transform = map(
                weighting_scheme = !! weighting_scheme_syms
              )
            ),
  
  res_sibs_parent_weights = map_dfr(score_names, 
                              run_regs_weighted, 
                              famhist     = famhist, 
                              weight_data = parent_weights,
                              dep.var     = "n_sibs"
                            ),
  
  res_children_comparison = map_dfr(score_names,
                              run_regs_weighted,
                              famhist     = famhist %>% filter(n_children > 0),
                              weight_data = age_qual_weights,
                              dep.var     = "n_children",
                              subset      = quote(kids_ss)
                            ),
  
  res_period =  target(
                  map_dfr(score_names, 
                    run_regs_period,
                    famhist     = famhist,
                    weight_data = weight_data,
                    children    = children
                  ),
                  transform = map(
                    children    = !! (c(TRUE, FALSE)),
                    weight_data = !! period_weight_syms,
                    .names      = c("res_period_children", "res_period_parents")
                  )
                ),
    
  res_sex = {
    sexes <- rlang::exprs(
            sex == 0, 
            sex == 1
          )
    res_sex <- expand_grid(
            score_name = score_names, 
            subset     = sexes
          ) %>% 
            pmap_dfr(
              run_regs_subset, 
              famhist = famhist %>% filter(kids_ss)
            )
    res_sex$sex <- ifelse(res_sex$subset == "sex == 0", "Female", "Male") 
    
    int_pval <- map_dfr(score_names, ~ run_regs_fml(
                    fml        = "n_children ~ {score_name}*I(sex==0)",
                    score_name = .x,
                    famhist    = famhist,
                    subset     = quote(kids_ss)
                  ),
                  .id = "score_name"
                ) %>% 
                filter(grepl(":", term)) %>% 
                select(score_name, diff.p.value = p.value)
    
    res_sex %<>% left_join(int_pval, by = c("term" = "score_name"))
    res_sex
  },
  
  res_age_flb = {
    res <- map_dfr(score_names, 
            ~run_regs_fml(
              fml        = "n_children ~ {score_name} + age_flb",
              score_name = .x, 
              famhist    = famhist,
              subset     = quote(kids_ss) 
            ), 
            .id = "score_name"
          )
    res %<>% filter(term != "(Intercept)")
    
    res
  },
  
  res_age_flb_cross = {
    res <- map_dfr(score_names, 
            ~run_regs_fml(
              fml        = "n_children ~ age_flb_cat + {score_name}:age_flb_cat",
              score_name = .x,
              famhist    = famhist,
              subset     = quote(kids_ss)
            ),
            .id = "score_name"
          )
    res %<>% filter(grepl(":", term))
          
  },
  
  res_age_flb_dv = {
    res <- map_dfr(score_names,
            ~run_regs_fml(
              fml        = "age_flb ~ {score_name}",
              score_name = .x, 
              famhist    = famhist,
              subset     = quote(kids_ss)
            )
          )
    res %<>% filter(term != "(Intercept)")
    
    res
  },
  
  res_age_birth_parents = {
    vars <- expand_grid(
            score_name = unname(score_names), # just easier to avoid hassle... 
            control    = c("fath_age_birth", "moth_age_birth")
          )
    res <- pmap_dfr(vars,
            run_regs_fml,
            fml = "n_sibs ~ {score_name} + {control}",
            subset = quote(birth_order == 1),
            famhist = famhist,
            .id = "id"
          )
    
    res$score_name <- vars$score_name[as.numeric(res$id)]
    res$control    <- vars$control[as.numeric(res$id)]
    res$id <- NULL
    res %<>% filter(term != "(Intercept)")
    
    res
  },
  
  res_townsend_parents = {
    res <-  map_dfr(score_names,
              ~run_regs_fml(
                "n_sibs ~ Quin71 + {score_name}:Quin71",
                score_name = .x,
                famhist    = famhist_townsend_71
              ),
              .id = "score_name"
            )
    res %<>% filter(grepl(":", term))
    res
  },
  
  res_age_birth_parents_dv = {
    vars <- expand_grid(
            score_name = unname(score_names), 
            dep.var    = c("fath_age_birth", "moth_age_birth")
          )
    res <- pmap_dfr(vars,
            run_regs_fml,
            fml = "{dep.var} ~ {score_name}",
            subset = quote(birth_order == 1),
            famhist = famhist,
            .id = "id"
          )
    
    res$score_name <- vars$score_name[as.numeric(res$id)]
    res$dep.var    <- vars$dep.var[as.numeric(res$id)]
    res$id <- NULL
    res %<>% filter(term != "(Intercept)")
    
    res
  },
  
  res_partners = {
    sexes <- rlang::exprs(sex == 0, sex == 1)
    pars <- expand_grid(subset = sexes, score_name = score_names)
    res <- pmap_dfr(pars,
            ~ run_regs_fml(
              fml        =
                "n_children ~ lo_partners + {score_name}:lo_partners",
              score_name = .y,
              famhist    = famhist %>% filter(kids_ss),
              subset     = .x
            ),
              .id = "row_number"
          )
    pars$row_number <- as.character(seq_len(nrow(pars)))
    
    left_join(res, pars, by = "row_number") %>%
          mutate(sex = ifelse(subset == "sex == 1", "Male", "Female")) %>%
          select(-row_number, -subset)
  },
  
  res_with_partner = {
    res <- map_dfr(score_names,
            ~run_regs_fml(
              "n_children ~ with_partner + {score_name}:with_partner",
              score_name = .x,
              famhist    = famhist,
              subset     = quote(kids_ss)
            ),
            .id = "score_name"
          )
    res %<>% filter(term != "(Intercept)")
    res
  },
  
  res_with_partner_sex = {
    sexes <- rlang::exprs(sex == 0, sex == 1)
    pars <- expand_grid(subset = sexes, score_name = score_names)
    res <- pmap_dfr(pars,
            ~ run_regs_fml(
              fml        =
                "n_children ~ with_partner + {score_name}:with_partner",
              score_name = .y,
              famhist    = famhist %>% filter(kids_ss),
              subset     = .x
            ),
              .id = "row_number"
          )
    pars$row_number <- as.character(seq_len(nrow(pars)))
    
    left_join(res, pars, by = "row_number") %>%
          mutate(sex = ifelse(subset == "sex == 1", "Male", "Female")) %>%
          select(-row_number, -subset)
  },
  
  res_with_partner_narrow = {
    mf_pairs_tmp <- mf_pairs %>% select(ID.m, ID.f, n_children.m, n_children.f)
    famhist_tmp <- famhist %>% 
                      filter(sex == 1) %>% 
                      select(f.eid, with_partner, n_children) %>% 
                      left_join(mf_pairs_tmp, by = c("f.eid" = "ID.m"))
    famhist_tmp2 <- famhist %>% 
                      filter(sex == 0) %>% 
                      select(f.eid, with_partner, n_children) %>% 
                      left_join(mf_pairs_tmp, by = c("f.eid" = "ID.f"))
    famhist_tmp  %<>% bind_rows(famhist_tmp2)
    
    famhist_tmp %<>% 
          filter(n_children.f == n_children.m, with_partner) %>% 
          group_by(f.eid) %>% 
          filter(n() == 1)
    
    famhist_tmp$with_partner_narrow <- TRUE
    famhist_tmp %<>% select(f.eid, with_partner_narrow)
    
    famhist %<>% left_join(famhist_tmp)
    famhist %<>% mutate(with_partner_narrow = ! is.na(with_partner_narrow))
    
    res <- map_dfr(score_names,
            ~run_regs_fml(
              "n_children ~ with_partner_narrow + {score_name}:with_partner_narrow",
              score_name = .x,
              famhist    = famhist,
              subset     = quote(kids_ss)
            ),
            .id = "score_name"
          )
    res %<>% filter(term != "(Intercept)")
    res
  },
  
  res_edu = {
    subsets <- rlang::exprs(
            age_fte_cat == "< 16",
            age_fte_cat == "16-18",
            age_fte_cat == "> 18"
          )
    # putting subset first matters below, because it is unnamed:
    pars <- expand_grid(subset = subsets, score_name = score_names) 
    res_edu <- pmap_dfr(pars,
            ~ run_regs_fml(
              "n_children ~ {score_name}", 
              score_name = .y, 
              subset     = .x,
              famhist    = famhist %>% filter(kids_ss)
            ),
            .id = "row_number"
          ) 
    pars$row_number <- as.character(seq_len(nrow(pars)))
    res_edu %>% left_join(pars, by = "row_number") %>% 
          filter(term != "(Intercept)") %>% 
          mutate(
            age_fte_cat = sub("age_fte_cat == \"(.*)\"", "\\1", subset)
          ) %>% 
          select(-row_number, -subset)
  },
  
  res_income = {
    subsets = rlang::exprs(
            income_cat == 1,
            income_cat == 2,
            income_cat == 3,
            income_cat == 4,
            income_cat == 5
          ) 
    # keep variables in this order (it matters to .id):
    pars <- expand_grid(subset = subsets, score_name = score_names) 
    res_income <- pmap_dfr(pars,
            ~ run_regs_fml(
              "n_children ~ {score_name}", 
              score_name = .y, 
              subset     = .x,
              famhist    = famhist %>% filter(kids_ss)
            ),
            .id = "row_number"
          ) 
    pars$row_number <- as.character(seq_len(nrow(pars)))
    res_income %>% left_join(pars, by = "row_number") %>% 
          filter(term != "(Intercept)") %>% 
          mutate(
            income_cat = sub("income_cat == (.*)", "\\1", subset)
          ) %>% 
          select(-row_number, -subset)
  },
  
  res_income_controlled = {
    res <- map_dfr(score_names, run_age_anova, 
                     control = "factor(income_cat)", 
                     famhist = famhist, 
                    .id = "score_name"
                   )
    res %>% filter(grepl(":factor\\(income_cat\\)", term))
  },
  
  res_edu_controlled = {
    res <- map_dfr(score_names, run_age_anova, 
                     control = "factor(age_fte_cat)", 
                     famhist = famhist, 
                    .id = "score_name"
                   )
    res %>% filter(grepl(":factor\\(age_fte_cat\\)", term))
  },
  
  res_ee_control = {
    famhist <- add_ashe_income(famhist, ashe_income)
    map_dfr(score_names, 
              ~ run_regs_fml(
                fml = "n_children ~ {score_name} + age_fte_cat + first_job_pay",
                score_name = .x,
                famhist    = famhist,
                subset     = quote(kids_ss)
              ), 
            .id = "score_name")
  },
  
    res_risk_control = {
      # f.2040.0.0 is risk attitude (via questionnaire)
      map_dfr(score_names, 
                ~ run_regs_fml(
                  fml = "n_children ~ {score_name} + f.2040.0.0",
                  score_name = .x,
                  famhist    = famhist,
                  subset     = quote(kids_ss)
                ), 
              .id = "score_name")
  },
  
  res_ineq = {
    famhist %>% 
          left_join(age_qual_weights, by = "f.eid") %>% 
          filter(! is.na(n_children), ! is.na(income_cat), ! is.na(weights)) %>% 
          group_by(income_cat) %>% 
          mutate(
            parent_weights = weights/n_sibs,
            parent_weights = ifelse(is.na(parent_weights), 0, parent_weights)
          ) %>% 
          summarize(across(all_of(score_names), 
                    list(
                      unw     = ~ weighted.mean(.x, weights, na.rm = TRUE),
                      wtd     = ~ weighted.mean(.x, weights * n_children, 
                                              na.rm = TRUE),
                      parents = ~ weighted.mean(.x, parent_weights, 
                                              na.rm = TRUE)
                    )
          )) %>% 
          tidyr::pivot_longer(-income_cat, 
                                names_to = c("score", ".value"), 
                                names_pattern = "(.*)_(.+?)$"
                              ) 
  },
  
  # res_margins = {
  #    res_extensive <- map_dfr(score_names, 
  #                       ~ run_regs_fml(
  #                         fml = "I(n_children > 0) ~ factor(sex) + {score_name}:factor(sex)",
  #                         score_name = .x,
  #                         famhist    = famhist,
  #                         subset     = quote(kids_ss)
  #                       ), 
  #                     .id = "score_name")
  #    res_intensive <- map_dfr(score_names, 
  #                       ~ run_regs_fml(
  #                         fml = "n_children ~ factor(sex) + {score_name}:factor(sex)",
  #                         score_name = .x,
  #                         famhist    = famhist,
  #                         subset     = quote(kids_ss & n_children > 0)
  #                       ), 
  #                     .id = "score_name")
  #    res <- bind_rows(
  #                      extensive = res_extensive, 
  #                      intensive = res_intensive, 
  #                      .id = "type"
  #                    )
  #    res %<>% 
  #      filter(term != "(Intercept)") %>% 
  #      arrange(score_name, type)
  #    
  #    res
  # },
  
  # res_together = {
  #   most_score_names <- setdiff(score_names, 
  #         c("EA2_noUKB", "whr_combined", "wc_combined", "hip_combined"))
  #   all_pgs <- paste0(most_score_names, collapse = " + ")
  #   dep_vars <- c("n_children", "n_sibs")
  #   names(dep_vars) <- dep_vars
  #   
  #   fml <- glue("{{dep_var}} ~ {all_pgs}") # double to quote and use below
  #   res_all_pgs <- map_dfr(dep_vars,
  #           ~run_regs_fml(fml, 
  #                           dep_var = .x, 
  #                           famhist = famhist, 
  #                           subset  = if (.x == "n_children") {
  #                                       quote(kids_ss)
  #                                     } else {
  #                                       NULL
  #                                     }
  #                         ),
  #           .id = "dep.var"
  #         )
  #   famhist <- join_famhist_pcs(famhist, pcs)
  #   all_pcs <- grep("PC", names(pcs), value = TRUE)
  #   all_pcs <- paste(all_pcs, collapse = " + ")
  #   fml_pcs <- glue("{{dep_var}} ~ {all_pgs} + {all_pcs}")
  #   res_all_pgs_pcs <- map_dfr(dep_vars,
  #           ~run_regs_fml(
  #             fml_pcs, 
  #             dep_var = .x, 
  #             famhist = famhist,
  #             subset  = if (.x == "n_children") quote(kids_ss) else NULL
  #             ),
  #           .id = "dep.var"
  #         )
  #   
  #   bind_rows("No" = res_all_pgs, "Yes" = res_all_pgs_pcs, .id = "PCs") %>% 
  #         filter(term != "(Intercept)", ! grepl("^PC", term))
  # },
  
  # crashes out of memory if we run the loop within the plan?
  # res_mnl = run_regs_mnlogit(score_names, fhl_mlogit),
  
  report =  {
              rmarkdown::render(
                input       = knitr_in("why-natural-selection.Rmd"), 
                output_file = file_out("why-natural-selection.pdf"), 
                quiet       = TRUE
              )
            }
  
)

drake_config(plan, history = FALSE, log_build_times = FALSE)
