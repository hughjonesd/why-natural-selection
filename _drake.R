
suppressPackageStartupMessages({
  library(drake)
  library(rlang)
  library(purrr)
  library(tidyr)
  library(readr)
})


source("make-data.R")
source("make-census-weights.R")
source("regressions.R")
source("weight-data.R")

data_dir           <- "../negative-selection-data"
pgs_dir            <- file.path(data_dir, "polygenic_scores/")
ph_file            <- file.path(data_dir, "UKB.EA_pheno.coordinates.QC.csv")
pcs_file           <- file.path(data_dir, "UKB.HM3.100PCs.40310.txt")
famhist_file       <- file.path(data_dir, "david.family_history.traits.out.csv")
famhist2_file      <- file.path(data_dir, "david.family_history.traits.20042020.out.csv")
famhist3_file      <- file.path(data_dir, "david.family_history.traits.05052020.out.csv")
famhist4_file      <- file.path(data_dir, "david.family_history.traits.16052020.out.csv")
famhist5_file      <- file.path(data_dir, "david.family_history.traits.18052020.out.csv")
famhist6_file      <- file.path(data_dir, "david.family_history.traits.17062020.out.csv")
famhist7_file      <- file.path(data_dir, "david.birthinfo.traits.14072020.out.csv")
famhist8_file      <- file.path(data_dir, "david.traits.03112020.out.csv")
rgs_file           <- file.path(data_dir, "EA3_rgs.10052019.rgs.csv")
mf_pairs_file      <- file.path(data_dir, "spouse_pair_info", 
                        "UKB_out.mf_pairs_rebadged.csv")
nomis_file         <- file.path(data_dir, "nomis-2011-statistics.csv")
ghs_file           <- file.path(data_dir, "UKDA-5804-stata8", "stata8", 
                        "Ghs06client.dta")
DC1108EW_file      <- file.path(data_dir, "DC1108EW.csv")
DC1108EW_msoa_file <- file.path(data_dir, "DC1108EW-msoas.csv")
DC5202SC_file      <- file.path(data_dir, "DC5202SC.csv")
msoa_shapefile     <- file.path(data_dir, "infuse_msoa_lyr_2011_clipped", 
                        "infuse_msoa_lyr_2011_clipped.shp")
ashe_income_file   <- file.path(data_dir, 
                        "SOC-income", 
                        "Occupation (4) Table 14.7a   Annual pay - Gross 2007.xls") 


weighting_schemes <- rlang::syms(c("flb_weights", "age_qual_weights", 
                        "msoa_weights"))
period_weight_syms <- rlang::syms((c("age_qual_weights", "parent_weights")))

plan <- drake_plan(
  score_names  = {
    score_names <- sub(
      ".*UKB\\.AMC\\.(.*?)\\..*", 
      "\\1", 
      list.files(file_in(!! pgs_dir), pattern = "csv$"), 
      perl = TRUE
    )
    setNames(score_names, score_names)
  },
  
  famhist_raw  = target(
                   make_famhist( 
                      file_in(!! ph_file), 
                      file_in(!! famhist_file),
                      file_in(!! famhist2_file),
                      file_in(!! famhist3_file),
                      file_in(!! famhist4_file),
                      file_in(!! famhist5_file),
                      file_in(!! famhist6_file),
                      file_in(!! famhist7_file),
                      file_in(!! famhist8_file),
                      file_in(!! pgs_dir)
                    ), 
                    format = "fst"
                  ),
  
  pcs          = read_table2(file_in(!! pcs_file)),
  
  ashe_income = make_ashe_income(file_in(!! ashe_income_file)),
  
  famhist      =  target(
                   edit_famhist(famhist_raw, score_names, ashe_income), 
                   format = "fst"
                  ),
  
  famhist_msoa = target(
                   find_containing_msoas(famhist, msoa_shapefile), 
                   format = "fst"
                 ),
  
  # mlogit now uses dfidx to create a "long" dataset, but mnlogit 
  # (used below) doesn't update. Skipping for now.
  # fhl_mlogit   =  make_famhist_long_mlogit(famhist, score_names),
  
  resid_scores = target(
                   make_resid_scores(famhist, pcs, score_names), 
                   format = "fst"
                 ),

  ghs_subset   = target(make_ghs_subset(file_in(!! ghs_file)), format = "fst"),
  
  pgs_over_time = calc_pgs_over_time(famhist, score_names),
  
  census_age_qual = target(
                      make_census_age_qual(
                        file_in(!! DC1108EW_file),
                        file_in(!! DC5202SC_file)
                      ), 
                      format = "fst"
                    ),
  
  census_msoa = target(
                  make_census_msoa(file_in(!! DC1108EW_msoa_file)), 
                  format = "fst"
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
  
  mf_pairs =  target(
                make_mf_pairs(file_in(!! mf_pairs_file), famhist), 
                format = "fst"
              ), 
  
  rgs = make_rgs(file_in(!! rgs_file)),
  
  res_all = {
    famhist <- join_famhist_resid_scores(famhist, resid_scores)
    res_sibs <- run_regs_basic(
            dep_var     = "n_sibs", 
            score_names = score_names,
            famhist     = famhist
          )
    res_chn <- run_regs_basic(
            dep_var     = "n_children", 
            score_names = score_names, 
            famhist     = famhist
          )
    dplyr::bind_rows(
            "N siblings" = res_sibs, 
            "N children" = res_chn, 
            .id = "dep.var"
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
      pcs         = pcs
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
                famhist = famhist,
                weight_data = weighting_scheme,
                dep.var = "n_children"
              ),
              transform = map(
                weighting_scheme = !! weighting_schemes
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
            pmap_dfr(run_regs_subset, famhist = famhist)
    res_sex$sex <- ifelse(res_sex$subset == "sex == 0", "Female", "Male") 
    
    int_pval <- map_dfr(score_names, ~ run_regs_fml(
                    fml        = "n_children ~ {score_name}*I(sex==0)",
                    score_name = .x,
                    famhist    = famhist
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
              famhist    = famhist
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
              famhist    = famhist
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
              famhist    = famhist
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
            subset = quote(n_older_sibs == 0),
            famhist = famhist,
            .id = "id"
          )
    
    res$score_name <- vars$score_name[as.numeric(res$id)]
    res$control    <- vars$control[as.numeric(res$id)]
    res$id <- NULL
    res %<>% filter(term != "(Intercept)")
    
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
            subset = quote(n_older_sibs == 0),
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
              famhist    = famhist,
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
              famhist    = famhist
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
              famhist    = famhist
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
              famhist    = famhist
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
    res <- map_dfr(score_names, 
            ~run_regs_fml(
              "n_children ~ {score_name}*(income_cat + YOB +I(YOB^2))", 
              score_name = .x, famhist = famhist
            ),
            .id = "score_name"
          )
    res %>% filter(grepl(":income_cat", term))
  },
  
  res_edu_controlled = {
    res <- map_dfr(score_names, 
      ~run_regs_fml(
        "n_children ~ {score_name}*(as.numeric(age_fte_cat) + YOB + I(YOB^2))", 
        score_name = .x, famhist = famhist
      ),
      .id = "score_name"
    )
    res %>% filter(grepl(":.*age_fte_cat", term))
  },
  
  res_partners_controlled = {
    sexes <- rlang::exprs(sex == 0, sex == 1)
    pars <- expand_grid(subset = sexes, score_name = score_names)
    res <- pmap_dfr(pars,
            ~ run_regs_fml(
              fml        =
                "n_children ~ lo_partners + lo_partners:({score_name} + YOB + I(YOB^2))",
              score_name = .y,
              famhist    = famhist,
              subset     = .x
            ),
              .id = "row_number"
          )
    pars$row_number <- as.character(seq_len(nrow(pars)))
    
    left_join(res, pars, by = "row_number") %>%
          mutate(sex = ifelse(subset == "sex == 1", "Male", "Female")) %>%
          select(-row_number, -subset) %>% 
          filter(! grepl("YOB", term))
  },
  
  res_together = {
    most_score_names <- setdiff(score_names, 
          c("EA2_noUKB", "whr_combined", "wc_combined", "hip_combined"))
    all_pgs <- paste0(most_score_names, collapse = " + ")
    dep_vars <- c("n_children", "n_sibs")
    names(dep_vars) <- dep_vars
    
    fml <- glue("{{dep_var}} ~ {all_pgs}") # double to quote and use below
    res_all_pgs <- map_dfr(dep_vars,
            ~run_regs_fml(fml, dep_var = .x, famhist = famhist),
            .id = "dep.var"
          )
    famhist <- join_famhist_pcs(famhist, pcs)
    all_pcs <- grep("PC", names(pcs), value = TRUE)
    all_pcs <- paste(all_pcs, collapse = " + ")
    fml_pcs <- glue("{{dep_var}} ~ {all_pgs} + {all_pcs}")
    res_all_pgs_pcs <- map_dfr(dep_vars,
            ~run_regs_fml(fml_pcs, dep_var = .x, famhist = famhist),
            .id = "dep.var"
          )
    
    bind_rows("No" = res_all_pgs, "Yes" = res_all_pgs_pcs, .id = "PCs") %>% 
          filter(term != "(Intercept)", ! grepl("^PC", term))
  },
  
  # crashes out of memory if we run the loop within the plan?
  # res_mnl = run_regs_mnlogit(score_names, fhl_mlogit),
  
  report = rmarkdown::render(
          input       = knitr_in("why-negative-selection.Rmd"), 
          output_file = file_out("why-negative-selection.pdf"), 
          quiet       = TRUE
        )
  
)


drake_config(plan, history = FALSE)
