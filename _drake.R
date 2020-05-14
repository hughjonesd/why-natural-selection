
suppressPackageStartupMessages({
  library(drake)
  library(rlang)
  library(purrr)
  library(tidyr)
})

source("make-data.R")
source("regressions.R")

pgs_dir       <- "~/UKBB data 2019/polygenic_scores/"
ph_file       <- "~/UKBB data 2019/UKB.EA_pheno.coordinates.QC.david.csv"
pcs_file      <- "~/UKBB data 2019/ukb30545.40PCs.csv"
famhist_file  <- "~/Dropbox/assortative mating/biobank-analysis/david.family_history.traits.out.csv"
famhist2_file <- "~/UKBB data 2019/david.family_history.traits.20042020.out.csv"
famhist3_file <- "~/UKBB data 2019/david.family_history.traits.05052020.out.csv"
rgs_file      <- "EA3_rgs.10052019.rgs.csv"
mf_pairs_file <- "~/UKBB data 2019/spouse_pair_info/UKB_out.mf_pairs_rebadged.csv"


plan <- drake_plan(
  
  
  famhist      = edit_famhist(famhist_raw, reverse_code),
  
  
  famhist_raw  = target(
                   make_famhist( 
                      file_in(!! ph_file), 
                      file_in(!! pcs_file), 
                      file_in(!! famhist_file),
                      file_in(!! famhist2_file),
                      file_in(!! famhist3_file),
                      file_in(!! pgs_dir)
                    ), 
                    format = "fst"
                  ),
  
  
  mf_pairs     = target(
                   make_mf_pairs(file_in(!! mf_pairs_file), famhist), 
                   format = "fst"
                 ), 
  
  
  rgs          = make_rgs(file_in(!! rgs_file)),
   
  
  score_names  = {
    sub(
            ".*UKB\\.AMC\\.(.*?)\\..*", 
            "\\1", 
            list.files(file_in(!! pgs_dir), pattern = "csv$"), 
            perl = TRUE
          )
  },
  
  
  reverse_code = c(), # don't reverse code anything
  # {
  #   setdiff(score_names, c("agreeableness", "age_at_menarche",
  #           "age_at_menopauze", "cognitive_ability", 
  #           "conscientiousness", "EA2_noUKB", "EA3_excl_23andMe_UK", 
  #           "extraversion", "height_combined", "openness")
  #         )
  # },
  
  
  res_all      = {
    res_sibs <- map_dfr(score_names, run_regs, 
            dep_var = "n_sibs", 
            famhist = famhist
          )
    res_chn <- map_dfr(score_names, run_regs, 
            dep_var = "n_children", 
            famhist = famhist
          )
    dplyr::bind_rows(
            "N siblings" = res_sibs, 
            "N children" = res_chn, 
            .id = "dep.var"
          )
  },
  
  
  res_pcs     = {
    pcs <- paste0("PC", 1:40)
    res_sibs_pcs <- map_dfr(pcs, run_regs_pcs, 
      dep_var     = "n_sibs", 
      famhist     = famhist
    )
    res_chn_pcs <- map_dfr(pcs, run_regs_pcs, 
      dep_var     = "n_children", 
      famhist     = famhist
    )
    dplyr::bind_rows(
      "N siblings" = res_sibs_pcs, 
      "N children" = res_chn_pcs, 
      .id = "dep.var"
    ) 
  },
  
  
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
         res_sex$sex <- ifelse(res_sex$subset == "sex == 0", 
           "Female", "Male") 
   res_sex
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
    
    all_pcs <- paste0("PC", 1:40, collapse = " + ")
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


drake_config(plan)
