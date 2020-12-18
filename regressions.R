
suppressPackageStartupMessages({
  library(glue)
  library(broom)
  library(magrittr)
  library(santoku)
  library(tidyr)
  library(dplyr)
  library(lmtest)
  loadNamespace("tibble")
  loadNamespace("ggplot2")
  loadNamespace("sampleSelection")
  loadNamespace("Formula")
  loadNamespace("fixest")
})

calc_pgs_over_time <- function (famhist, score_names) {
  pgs_over_time <- famhist %>% 
        select(YOB, all_of(score_names)) %>% 
        group_by(
          YOB = santoku::chop(YOB, 
            seq(1940, 1970, 5), 
            labels    = seq(40, 65, 5),
            extend    = FALSE
          ) 
        ) %>% 
        filter(! is.na(YOB)) %>% 
        summarize(across(all_of(score_names), 
          ggplot2::mean_cl_normal, na.rm = TRUE
        )) %>% 
        pivot_longer(-YOB, names_to = "score_name", values_to = "score") %>% 
        mutate(
          score_lo = score$ymin,
          score_hi = score$ymax,
          score    = score$y
        )
  
  run_reg <- function (score_name) {
    fml <- as.formula(glue("{score_name} ~ YOB"))
    tidy(lm(fml, famhist))  
  }
  time_regs <- map_dfr(score_names, run_reg, .id = "score_name")
  time_regs %<>% 
        filter(term != "(Intercept)") %>% 
        select(score_name, p.value, statistic) %>% 
        mutate(
          Change = case_when(
            p.value > 0.05/33 ~ "o", 
            statistic > 0     ~ "+",
            statistic <= 0    ~ "-"
        ))
  
  pgs_over_time %<>% left_join(time_regs, by = "score_name")
  
  pgs_over_time
}


run_regs_basic <- function (dep_var, score_names, famhist, subset = NULL) {
  
  run_reg <- function (score_name) {
    fml <- as.formula(glue("{dep_var} ~ {score_name}_raw"))
    reg_bv <- tidy(lm(fml, famhist, subset = eval(subset)), conf.int = TRUE) 
    reg_bv <- mutate(reg_bv, term = gsub("_raw", "", term))
    reg_bv
  }
  basic_res <- map_dfr(score_names, run_reg, .id = "score_name")
  
  run_resid_reg <- function (score_name) {
    fml_resid <- as.formula(glue("{dep_var} ~ {score_name}"))
    reg_resid <- tidy(lm(fml_resid, famhist, subset = eval(subset)), conf.int = TRUE) 
    reg_resid
  }
  resid_res <- map_dfr(score_names, run_resid_reg, .id = "score_name")
  
  res <- bind_rows(raw = basic_res, controlled = resid_res, .id = "reg.type")
  res %<>% filter(term != "(Intercept)")
  
  return(res)
}


run_regs_pcs <- function (dep_var, famhist, pcs, subset = NULL) {
  pc_names <- grep("PC", names(pcs), value = TRUE)
  
  fh_short <- join_famhist_pcs(famhist[c("f.eid", dep_var, "kids_ss")], pcs)
  
  run_reg_pc <- function(pc_name) {
    fml <- as.formula(glue("{dep_var} ~ {pc_name}"))
    reg <- tidy(lm(fml, fh_short, subset = eval(subset)), conf.int = TRUE) 
    reg %<>% filter(term != "(Intercept)")
  }
  
  map_dfr(pc_names, run_reg_pc, .id = "pc_name")
}


run_regs_weighted <- function (score_name, famhist, weight_data, dep.var, subset = NULL) {
  famhist <- inner_join(famhist, weight_data, by = "f.eid")

  fml <- as.formula(glue("{dep.var} ~ {score_name}"))
  mod <- lm(fml, data = famhist, weights = weights, na.action = na.exclude, 
            subset = eval(subset))
  res <- tidy(mod, conf.int = TRUE)
  res %<>% filter(term != "(Intercept)")
  attr(res, "call") <- mod$call
  
  return(res)
}


run_regs_period <- function (children, score_name, famhist, weight_data) {
  
  famhist <- inner_join(famhist, weight_data, by = "f.eid")
  
  dep_var  <- if (children) "n_children" else "n_sibs"
  subset   <- if (children) quote(kids_ss) else NULL
  famhist$year_split <- famhist$YOB >= 1950
  # evaluates to 1950:
  # median(famhist$YOB, na.rm = TRUE)
  famhist$year_split <- factor(famhist$year_split, labels = c("early", "late"))
  fml <- as.formula(glue("{dep_var} ~ year_split + {score_name}:year_split"))
  mod <- lm(fml, famhist, weights = weights, subset = eval(subset))
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


run_regs_fml <- function(fml, ..., subset = NULL, famhist, weights = NULL) {
  glue_args <- list(...)
  fml <- as.formula(glue::glue_data(fml, .x = glue_args))
  mod <- lm(fml, famhist, subset = eval(subset), weights = eval(weights))
  res <- tidy(mod, conf.int = TRUE)
  attr(res, "call") <- mod$call
  res
}


run_regs_mnlogit <- function (score_names, fhl_mlogit) {
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


run_age_anova <- function (score_name, control, famhist) {
  fml <- as.formula(glue::glue(
    "n_children ~ {score_name}*({control} + 
          age_at_recruitment + I(age_at_recruitment^2))"))
  res <- lm(fml, data = famhist, subset = kids_ss) %>% 
           car::Anova() %>% 
           tidy()
}


run_reg_income_dv <- function (famhist, score_names, ashe_income, age_qual_weights) {
  famhist %<>% add_ashe_income(ashe_income)
  famhist %<>% left_join(age_qual_weights, by = "f.eid")
  
  most_score_names <- setdiff(score_names, c("EA2_noUKB", "whr_combined", 
                                "wc_combined", "hip_combined", "body_fat"))
  
  fml <- paste(most_score_names, collapse = " + ")
  fml <- glue::glue("log(cur_job_pay) ~ age_at_recruitment + ",
                      "I(age_at_recruitment^2) + {fml} | factor(sib_group)")
  fml <- Formula::as.Formula(fml)
  
  res <- fixest::feols(fml, data = famhist, weights = ~ weights)
  
  # predict manually, using causal coefficients from within-family regs
  # and choosing intercept to minimize squared error
  scores <- famhist %>% 
              mutate(agesq = age_at_recruitment^2) %>% 
              select(age_at_recruitment, agesq, all_of(most_score_names)) %>% 
              as.matrix()
  pay_hat <- c(scores %*% coef(res))
  alpha <- mean(log(famhist$cur_job_pay) - pay_hat, na.rm = TRUE)
  pay_hat <- pay_hat + alpha
  
  famhist$pay_hat <- pay_hat
}

run_income_dist <- function (famhist) {
  # predict individual income from PGS
  # 3 densities:
  # - empirically from last_job_pay, weighted by e.g. age_qual
  # - predicted from reg of last_job_pay on PGS, same weighting
  #   - use within-family regressions, then work out the intercept separately
  # - predicted from regression, weighting multiplied by n_children
}
