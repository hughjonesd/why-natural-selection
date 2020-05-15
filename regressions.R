
suppressPackageStartupMessages({
  library(glue)
  library(broom)
  library(dplyr)
  library(lmtest)
})


run_regs <- function (dep_var, score_names, famhist, pcs) {
  
  run_reg <- function (score_name) {
    fml <- as.formula(glue("{dep_var} ~ {score_name}"))
    reg_bv <- tidy(lm(fml, famhist), conf.int = TRUE)  
  }
  basic_res <- map_dfr(score_names, run_reg, .id = "score_name")
  
  famhist <- join_famhist_pcs(famhist, pcs)
  pc_names <- grep("PC", names(pcs), value = TRUE)
  run_resid_reg <- function (score_name) {
    fml_resid <- as.formula(glue("{dep_var} ~ {score_name} + {pc_names}"))
    reg_resid <- tidy(lm(fml_resid, famhist), conf.int = TRUE)  
  }
  resid_res <- map_dfr(score_names, run_resid_reg, .id = "score_name")
  
  res <- bind_rows(raw = basic_res, controlled = resid_res, .id = "reg.type")
  res %<>% filter(term != "(Intercept)", ! grepl("^PC", term))
  
  return(res)
}


run_regs_pcs <- function (dep_var, famhist, pcs) {
  pc_names <- grep("PC", names(pcs), value = TRUE)
  
  fh_short <- join_famhist_pcs(famhist[c("f.eid", dep_var)], pcs)
  
  run_reg_pc <- function(pc_name) {
    fml <- as.formula(glue("{dep_var} ~ {pc_name}"))
    reg <- tidy(lm(fml, fh_short), conf.int = TRUE) 
    reg %<>% filter(term != "(Intercept)")
  }
  
  map_dfr(pc_names, run_reg_pc, .id = "pc_name")
}


run_regs_weighted <- function (score_name, famhist, weights) {
  famhist <- left_join(famhist, weights, by = "f.eid")
  
  fml <- as.formula(glue("n_children ~ {score_name}"))
  regs <- tidy(lm(fml, famhist, weights = educ_weight,
        na.action = na.exclude), conf.int = TRUE) 
  regs %<>% filter(term != "(Intercept)")
  
  return(regs)
}


run_regs_period <- function (children, score_name, famhist) {
  dep_var  <- if (children) "n_children" else "n_sibs"
  famhist$year_split <- famhist$YOB >= median(famhist$YOB, na.rm = TRUE)
  famhist$year_split <- factor(famhist$year_split, labels = c("early", "late"))
  fml <- as.formula(glue("{dep_var} ~ {score_name}:year_split"))
  mod <- lm(fml, famhist)
  res <- tidy(mod, conf.int = TRUE)
  res <- filter(res, grepl(score_name, term))
  
  test_spec <- glue("{score_name}:year_splitearly - {score_name}:year_splitlate = 0")
  test <- car::lht(mod, test_spec)
  pval <- test[["Pr(>F)"]]
  res$diff.p.value <- pval[2]
  
  res$children <- children
  
  res
}


run_regs_subset <- function(score_name, subset, famhist) {
  fml <- as.formula(glue::glue("n_children ~ {score_name}"))
  mod <- lm(fml, famhist, subset = eval(subset))
  res <- tidy(mod, conf.int = TRUE)
  res <- filter(res, term == {{score_name}})
  res$subset <- list(subset)

  res
}


run_regs_fml <- function(fml, ..., subset = NULL, famhist) {
  glue_args <- list(...)
  fml <- as.formula(glue::glue_data(fml, .x = glue_args))
  mod <- lm(fml, famhist, subset = eval(subset))
  res <- tidy(mod, conf.int = TRUE)
  
  res
}

