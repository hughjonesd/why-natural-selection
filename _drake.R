
suppressPackageStartupMessages({
  library(drake)
  library(rlang)
  library(purrr)
  library(tidyr)
  library(readr)
  library(dplyr)
  loadNamespace("future")
  loadNamespace("santoku")
})

source("make-data.R") 
# infrastructure available at github.com/hughjonesd-private/import-ukbb-data
source("~/import-ukbb-data/import-ukbb-data.R")
source("make-census-weights.R")
source("regressions.R")
source("weight-data.R")


# file locations for data ====

nomis_file         <- file.path(data_dir, "nomis-2011-statistics.csv")
ghs_file           <- file.path(data_dir, "UKDA-5804-stata8", "stata8", 
                        "Ghs06client.dta")
DC1108EW_file      <- file.path(data_dir, "DC1108EW.csv")
DC1108EW_msoa_file <- file.path(data_dir, "DC1108EW-msoas.csv")
DC5202SC_file      <- file.path(data_dir, "DC5202SC.csv")
msoa_shapefile     <- file.path(data_dir, "infuse_msoa_lyr_2011_clipped", 
                        "infuse_msoa_lyr_2011_clipped.shp")
dep_data_dir       <- file.path(data_dir, "deprivation-data")
fertility_data_dir <- file.path(data_dir, "fertility_PRS")


# plan ====

weighting_scheme_syms <- rlang::syms(c("flb_weights", "age_qual_weights", 
                        "msoa_weights"))
period_weight_syms <- rlang::syms((c("age_qual_weights", "parent_aq_weights")))

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
                     
                     famhist <- add_fertility_prs(famhist, fertility_data_dir)
                     
                     # we put these here because they're tightly tied to this
                     # particular project
                     mab_terciles <- quantile(famhist$moth_age_birth[famhist$birth_order == 1], 
                                              1:2/3, na.rm = TRUE)
                     famhist$moth_age_birth_cat <- santoku::chop(famhist$moth_age_birth, 
                                                                 mab_terciles, 
                                                                 labels = lbl_discrete("-"))
                     fab_terciles <- quantile(famhist$fath_age_birth[famhist$birth_order == 1], 
                                              1:2/3, na.rm = TRUE)
                     famhist$fath_age_birth_cat <- santoku::chop(famhist$fath_age_birth, 
                                                                 fab_terciles, 
                                                                 labels = lbl_discrete("-"))
                     
                     famhist
                   },
                    format = "fst_tbl"
                  ),
  famhist_kids = target(
                   famhist %>% filter(kids_ss),
                   format = "fst_tbl"
                 ),
  
  famhist_pw = target(
                 inner_join(famhist, parent_weights, by = "f.eid"),
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
  
  parent_weights = weight_parents(famhist),
  
  parent_aq_weights = weight_parents(famhist, input_weights = age_qual_weights),
  
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
            famhist     = famhist_pw,
            weights     = famhist_pw$weights
          )
    res_chn <- run_regs_basic(
            dep_var     = "n_children", 
            score_names = score_names, 
            famhist     = famhist_kids
          )
    dplyr::bind_rows(
            "N siblings" = res_sibs, 
            "N children" = res_chn, 
            .id = "dep.var"
          )
  },
  
  res_quadratic = {
    map_dfr(score_names, 
              ~run_regs_fml(
                             "n_children ~ {score_name} + I({score_name}^2)", 
                             score_name = .x,
                             famhist    = famhist_kids
                           ),
              .id = "score_name"
            )
  },
  
  res_sibs_quadratic = {
    map_dfr(score_names, 
            ~run_regs_fml(
              "n_sibs ~ {score_name} + I({score_name}^2)", 
              score_name = .x,
              famhist    = famhist_pw,
              weights    = quote(weights)
            ),
            .id = "score_name"
    )
  },
  
  res_pcs = {
    pc_names <- grep("PC", names(pcs), value = TRUE)
    res_sibs_pcs <- run_regs_pcs( 
      dep_var     = "n_sibs", 
      famhist     = famhist_pw,
      pcs         = pcs,
      weights     = quote(famhist$weights) # famhist will be famhist_pw in the function
    )
    
    res_chn_pcs <- run_regs_pcs(
      dep_var     = "n_children", 
      famhist     = famhist_kids,
      pcs         = pcs,
      weights     = NULL 
    )
    dplyr::bind_rows(
      "N siblings" = res_sibs_pcs, 
      "N children" = res_chn_pcs, 
      .id = "dep.var"
    ) 
  },
  
  res_wt =  target({
                map_dfr(score_names, 
                  run_regs_weighted,
                  famhist     = famhist_kids,
                  weight_data = weighting_scheme,
                  dep.var     = "n_children"
                )
              },
              transform = map(
                weighting_scheme = !! weighting_scheme_syms
              )
            ),
  
  res_sibs_parent_weights = map_dfr(score_names, 
                              run_regs_weighted, 
                              famhist     = famhist, 
                              weight_data = parent_aq_weights,
                              dep.var     = "n_sibs"
                            ),
  
  res_children_comparison = map_dfr(score_names,
                              run_regs_weighted,
                              famhist     = famhist_kids %>% filter(n_children > 0),
                              weight_data = age_qual_weights,
                              dep.var     = "n_children"
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
              famhist = famhist_kids
            )
    res_sex$sex <- ifelse(res_sex$subset == "sex == 0", "Female", "Male") 
    
    int_pval <- map_dfr(score_names, ~ run_regs_fml(
                    fml        = "n_children ~ {score_name}*I(sex==0)",
                    score_name = .x,
                    famhist    = famhist_kids
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
              famhist    = famhist_kids
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
              famhist    = famhist_kids
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
              famhist    = famhist_kids
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
            famhist = famhist_pw,
            weights = weights,
            .id = "id"
          )
    
    res$score_name <- vars$score_name[as.numeric(res$id)]
    res$control    <- vars$control[as.numeric(res$id)]
    res$id <- NULL
    res %<>% filter(term != "(Intercept)")
    
    res
  },
  
  res_age_flb_mothers_cross = {
    run_regs_age_flb_parents_cross(famhist_pw, score_names, "moth_age_birth_cat")
  },
  
  res_age_flb_fathers_cross = {
    run_regs_age_flb_parents_cross(famhist_pw, score_names, "fath_age_birth_cat")
  },
  
  res_townsend_parents = {
    famhist_townsend_71 %<>% inner_join(parent_weights, by = "f.eid")
    res <-  map_dfr(score_names,
              ~run_regs_fml(
                "n_sibs ~ Quin71 + {score_name}:Quin71",
                score_name = .x,
                famhist    = famhist_townsend_71,
                weights    = weights
              ),
              .id = "score_name"
            )
    res %<>% filter(grepl(":", term))
    quintiles <- famhist_townsend_71 %>%
                   filter(! is.na(n_sibs), ! is.na(whr_combined)) %>%
                   group_by(Quin71) %>%
                   count()
    res %<>% 
      mutate(
        quintile = sub("Quin71(\\d):.*", "\\1", term)
      ) %>%
      left_join(quintiles, by = c("quintile" = "Quin71"))
    
    res
  },
  
  res_age_birth_parents_dv = {
    vars <- expand_grid(
            score_name = unname(score_names), 
            dep.var    = c("fath_age_birth", "moth_age_birth")
          )
    res <- pmap_dfr(vars,
            run_regs_fml,
            fml     = "{dep.var} ~ {score_name}",
            subset  = quote(birth_order == 1),
            famhist = famhist_pw,
            weights = weights,
            .id     = "id"
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
              famhist    = famhist_kids,
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
              famhist    = famhist_kids
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
              famhist    = famhist_kids,
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
    
    famhist_kids %<>% left_join(famhist_tmp)
    famhist_kids %<>% mutate(with_partner_narrow = ! is.na(with_partner_narrow))
    
    res <- map_dfr(score_names,
            ~run_regs_fml(
              "n_children ~ with_partner_narrow + {score_name}:with_partner_narrow",
              score_name = .x,
              famhist    = famhist_kids
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
              famhist    = famhist_kids
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
              famhist    = famhist_kids
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
                     famhist = famhist_kids, 
                    .id = "score_name"
                   )
    res %>% filter(grepl(":factor\\(income_cat\\)", term))
  },
  
  res_edu_controlled = {
    res <- map_dfr(score_names, run_age_anova, 
                     control = "factor(age_fte_cat)", 
                     famhist = famhist_kids, 
                    .id = "score_name"
                   )
    res %>% filter(grepl(":factor\\(age_fte_cat\\)", term))
  },
  
  res_mediation = run_mediation(famhist_kids, score_names),
  
  res_risk_control = {
    # f.2040.0.0 is risk attitude (via questionnaire)
    map_dfr(score_names, 
              ~ run_regs_fml(
                fml = "n_children ~ {score_name} + f.2040.0.0",
                score_name = .x,
                famhist    = famhist_kids
              ), 
            .id = "score_name")
  },
  
  res_ineq = run_cor_income(famhist_kids, score_names, age_qual_weights),
  
  res_cor_income = {
    famhist <- add_ashe_income(famhist, ashe_income)
    cor(famhist[score_names], famhist$first_job_pay, use = "pairwise")
  },
  
  res_cor_educ = {
    cor(famhist[score_names], famhist$age_fulltime_edu, use = "pairwise")
  },
  
  res_ineq_ea3 = run_ineq_ea3_calcs(famhist_kids, age_qual_weights, h2 = 0.4),
  
  res_fe_fertility = run_reg_fe_fertility(famhist_kids, score_names),
  
  report =  {
              rmarkdown::render(
                input       = knitr_in("why-natural-selection.Rmd"), 
                output_file = file_out("why-natural-selection.pdf"), 
                quiet       = TRUE
              )
            }
  
)

drake_config(plan, history = FALSE, log_build_times = FALSE, keep_going = TRUE)
