
suppressPackageStartupMessages({
  library(glue)
  library(broom)
  library(dplyr)
  library(lmtest)
})


run_regs_basic <- function (dep_var, score_names, famhist) {
  
  run_reg <- function (score_name) {
    fml <- as.formula(glue("{dep_var} ~ {score_name}"))
    reg_bv <- tidy(lm(fml, famhist), conf.int = TRUE)  
  }
  basic_res <- map_dfr(score_names, run_reg, .id = "score_name")
  
  run_resid_reg <- function (score_name) {
    score_resid <- paste0(score_name, "_resid")
    fml_resid <- as.formula(glue("{dep_var} ~ {score_resid}"))
    reg_resid <- tidy(lm(fml_resid, famhist), conf.int = TRUE) 
    reg_resid <- mutate(reg_resid, term = gsub("_resid", "", term))
    reg_resid
  }
  resid_res <- map_dfr(score_names, run_resid_reg, .id = "score_name")
  
  res <- bind_rows(raw = basic_res, controlled = resid_res, .id = "reg.type")
  res %<>% filter(term != "(Intercept)")
  
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


run_regs_weighted <- function (score_name, famhist, weight_data, dep.var) {
  famhist <- inner_join(famhist, weight_data, by = "f.eid")

  fml <- as.formula(glue("{dep.var} ~ {score_name}"))
  mod <- lm(fml, data = famhist, weights = weights, na.action = na.exclude)
  res <- tidy(mod, conf.int = TRUE)
  res %<>% filter(term != "(Intercept)")
  attr(res, "call") <- mod$call
  
  return(res)
}


run_regs_period <- function (children, score_name, famhist) {
  dep_var  <- if (children) "n_children" else "n_sibs"
  famhist$year_split <- famhist$YOB >= 1950
  # evaluates to 1950:
  # median(famhist$YOB, na.rm = TRUE)
  famhist$year_split <- factor(famhist$year_split, labels = c("early", "late"))
  fml <- as.formula(glue("{dep_var} ~ year_split + {score_name}:year_split"))
  mod <- lm(fml, famhist)
  res <- tidy(mod, conf.int = TRUE)
  res <- filter(res, grepl(score_name, term))
  
  test_spec <- glue("year_splitearly:{score_name} - year_splitlate:{score_name} = 0")
  test <- car::lht(mod, test_spec)
  pval <- test[["Pr(>F)"]]
  res$diff.p.value <- pval[2]
  
  res$children <- children
  attr(res, "call") <- mod$call
  
  res
}


run_regs_subset <- function(score_name, subset, famhist) {
  fml <- as.formula(glue::glue("n_children ~ {score_name}"))
  mod <- lm(fml, famhist, subset = eval(subset))
  res <- tidy(mod, conf.int = TRUE)
  res <- filter(res, term == {{score_name}})
  res$subset <- list(subset)
  attr(res, "call") <- mod$call

  res
}


run_regs_fml <- function(fml, ..., subset = NULL, famhist) {
  glue_args <- list(...)
  fml <- as.formula(glue::glue_data(fml, .x = glue_args))
  mod <- lm(fml, famhist, subset = eval(subset))
  res <- tidy(mod, conf.int = TRUE)
  attr(res, "call") <- mod$call
  res
}


run_regs_mnlogit <- function(score_names, fhl_mlogit) {
  loadNamespace("mnlogit")
  run_1_mnlogit <- function(score_name) {
    fml <- as.formula(glue::glue("n_ch_fac ~ 1 | {score_name}"))
    res <- try(mod <- mnlogit::mnlogit(fml, fhl_mlogit, ncores = 2))
    if (inherits(res, "try-error")) return(res)
    mod$model <- NULL # otherwise it stores a copy of the entire data 8-(
    return(mod)
  }
  lapply(score_names, run_1_mnlogit)
}

