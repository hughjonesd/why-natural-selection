
suppressPackageStartupMessages({
  library(glue)
  library(broom)
  library(dplyr)
  library(lmtest)
})


run_regs <- function (dep_var, score_name, famhist) {
  fml <- as.formula(glue("{dep_var} ~ {score_name}"))
  reg_bv <- tidy(lm(fml, famhist), conf.int = TRUE) 
  
  pcs <- paste0("PC", 1:40, collapse = " + ")
  fml_resid <- as.formula(glue("{dep_var} ~ {score_name} + {pcs}"))
  reg_resid <- tidy(lm(fml_resid, famhist), conf.int = TRUE)
  
  regs <- bind_rows(raw = reg_bv, controlled = reg_resid, .id = "reg.type")
  regs %<>% filter(term != "(Intercept)", ! grepl("^PC", term))
  
  return(regs)
}


run_regs_pcs <- function (dep_var, pc_name, famhist) {
  fml <- as.formula(glue("{dep_var} ~ {pc_name}"))
  regs <- tidy(lm(fml, famhist), conf.int = TRUE) 
  regs %<>% filter(term != "(Intercept)")
  
  return(regs)
}


run_regs_period <- function (children, score_name, famhist) {
  dep_var  <- if (children) "n_children" else "n_sibs"
  year_var <- if (children) "YOB" else "parents_imp_YOB"
  famhist$year_split <- famhist[[year_var]] >= median(famhist[[year_var]], 
        na.rm = TRUE)
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

