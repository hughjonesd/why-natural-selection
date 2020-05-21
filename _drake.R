
suppressPackageStartupMessages({
  library(drake)
  library(rlang)
  library(purrr)
  library(tidyr)
  library(readr)
})


source("make-data.R")
source("regressions.R")
source("weight-data.R")

data_dir      <- "../negative-selection-data"
pgs_dir       <- file.path(data_dir, "polygenic_scores/")
ph_file       <- file.path(data_dir, "UKB.EA_pheno.coordinates.QC.david.csv")
pcs_file      <- file.path(data_dir, "UKB.HM3.100PCs.40310.txt")
famhist_file  <- file.path(data_dir, "david.family_history.traits.out.csv")
famhist2_file <- file.path(data_dir, "david.family_history.traits.20042020.out.csv")
famhist3_file <- file.path(data_dir, "david.family_history.traits.05052020.out.csv")
famhist4_file <- file.path(data_dir, "david.family_history.traits.16052020.out.csv")
famhist5_file <- file.path(data_dir, "david.family_history.traits.18052020.out.csv")
rgs_file      <- file.path(data_dir, "EA3_rgs.10052019.rgs.csv")
mf_pairs_file <- file.path(data_dir, "spouse_pair_info/UKB_out.mf_pairs_rebadged.csv")
nomis_file    <- file.path(data_dir, "nomis-2011-statistics.csv")
ghs_file      <- file.path(data_dir, "UKDA-5804-stata8/stata8/Ghs06client.dta")

plan <- drake_plan(
  
  
  score_names  = {
    sub(
      ".*UKB\\.AMC\\.(.*?)\\..*", 
      "\\1", 
      list.files(file_in(!! pgs_dir), pattern = "csv$"), 
      perl = TRUE
    )
  },
  
  
  famhist_raw  = target(
                   make_famhist( 
                      file_in(!! ph_file), 
                      file_in(!! famhist_file),
                      file_in(!! famhist2_file),
                      file_in(!! famhist3_file),
                      file_in(!! famhist4_file),
                      file_in(!! pgs_dir)
                    ), 
                    format = "fst"
                  ),
  
  
  pcs          = read_table2(file_in(!! pcs_file)),
  
  
  famhist      = target(
                  edit_famhist(famhist_raw, score_names), 
                  format = "fst"
                 ),
  
  
  resid_scores = target(
                   make_resid_scores(famhist, pcs, score_names), 
                   format = "fst"
                 ),

  
  ghs_subset   = target(make_ghs_subset(file_in(!! ghs_file)), format = "fst"),
  
  
  ghs_weights  =  weight_by_ghs(
                    ~ factor(sex) + age_at_recruitment + factor(edu_qual), 
                    ghs_subset, 
                    famhist
                  ),
  
  flb_weights = {
    weight_by_ghs(
            ~ age_at_recruitment + factor(flb_cat) + factor(edu_qual), 
            ghs_subset, 
            famhist
          )
  },
  
  mf_pairs     = target(
                   make_mf_pairs(file_in(!! mf_pairs_file), famhist), 
                   format = "fst"
                 ), 
  
  
  rgs          = make_rgs(file_in(!! rgs_file)),
  
  
  res_all      = {
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
  
  
  res_pcs     = {
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
  
  
  res_weighted =  map_dfr(score_names, 
                    run_regs_weighted, 
                    famhist     = famhist, 
                    weight_data = ghs_weights,
                    dep.var     = "n_children"
                  ),
  
  
  res_sibs_weighted = map_dfr(score_names, 
                        run_regs_weighted, 
                        famhist     = famhist, 
                        weight_data = ghs_weights,
                        dep.var     = "n_sibs"
                      ),
  
  
  res_flb_weights = map_dfr(score_names, 
                      run_regs_weighted, 
                      famhist     = famhist, 
                      weight_data = flb_weights,
                      dep.var     = "n_children"
                    ),

  
  res_period   = {
    expand_grid(
            children = c(TRUE, FALSE),
            score_name = score_names
          ) %>% 
          pmap_dfr(run_regs_period, famhist = famhist)
  },
  
  
  res_sex      = {
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
    
    res_sex
  },
  
  
  res_age_flb = {
    res <- map_dfr(setNames(score_names, score_names), 
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
    famhist$age_flb_cat <- santoku::chop_equally(famhist$age_flb, 3, 
          lbl_integer("-"))
    res <- map_dfr(setNames(score_names, score_names), 
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
            score_name = score_names, 
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
            score_name = score_names, 
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
    pars <- expand_grid(score_name = score_names, subset = sexes)
    res <- pmap_dfr(pars,
            ~ run_regs_fml(
              fml        =
                "n_children ~ {score_name}*lo_partners",
              score_name = .x,
              famhist    = famhist,
              subset     = .y
            ),
              .id = "combo"
          )
    pars$combo <- as.character(seq_len(nrow(pars)))
    
    left_join(res, pars, by = "combo") %>%
          mutate(sex = ifelse(subset == "sex == 1", "Male", "Female")) %>%
          select(-combo, -subset)
  }, 
  
  
  res_edu = {
    subsets <- rlang::exprs(
            age_fte_cat == "< 16",
            age_fte_cat == "16-18",
            age_fte_cat == "> 18"
          )
    pars <- expand_grid(score_name = score_names, subset = subsets) 
    res_edu <- pmap_dfr(pars,
            ~ run_regs_fml(
              "n_children ~ {score_name}", 
              score_name = .x, 
              subset     = .y,
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
    pars <- expand_grid(score_name = score_names, subset = subsets) 
    res_income <- pmap_dfr(pars,
            ~ run_regs_fml(
              "n_children ~ {score_name}", 
              score_name = .x, 
              subset     = .y,
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
  
  
  report       =  rmarkdown::render(
                    input       = knitr_in("why-negative-selection.Rmd"), 
                    output_file = file_out("why-negative-selection.pdf"), 
                    quiet       = TRUE
                  )
  
)


drake_config(plan, history = FALSE)
