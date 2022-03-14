
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
  loadNamespace("psych")
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


run_regs_basic <- function (dep_var, score_names, famhist, subset = NULL, 
                              weights = NULL) {
  
  run_reg <- function (score_name) {
    fml <- as.formula(glue("{dep_var} ~ {score_name}_raw"))
    reg_bv <- tidy(lm(fml, famhist, subset = eval(subset), 
                        weights = weights), conf.int = TRUE) 
    reg_bv <- mutate(reg_bv, term = gsub("_raw", "", term))
    reg_bv
  }
  basic_res <- map_dfr(score_names, run_reg, .id = "score_name")
  
  run_resid_reg <- function (score_name) {
    fml_resid <- as.formula(glue("{dep_var} ~ {score_name}"))
    reg_resid <- tidy(lm(fml_resid, famhist, subset = eval(subset), 
                           weights = weights), conf.int = TRUE) 
    reg_resid
  }
  resid_res <- map_dfr(score_names, run_resid_reg, .id = "score_name")
  
  res <- bind_rows(raw = basic_res, controlled = resid_res, .id = "reg.type")
  res %<>% filter(term != "(Intercept)")
  
  return(res)
}


run_regs_pcs <- function (dep_var, famhist, pcs, subset = NULL, weights = NULL) {
  pc_names <- grep("PC", names(pcs), value = TRUE)
  
  fh_short <- join_famhist_pcs(famhist[c("f.eid", dep_var, "kids_ss")], pcs)
  
  run_reg_pc <- function(pc_name) {
    fml <- as.formula(glue("{dep_var} ~ {pc_name}"))
    reg <- tidy(lm(fml, fh_short, subset = eval(subset), 
                     weights = eval(weights)), conf.int = TRUE) 
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
  
  dep_var  <- if (children) "RLRS" else "RLRS_parents"
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
  fml <- as.formula(glue::glue("RLRS ~ {score_name}"))
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


run_regs_age_flb_parents_cross <- function(famhist, score_names, age_birth_var) {
  fml <- sprintf("RLRS_parents ~ %s + {score_name}:%s", age_birth_var, age_birth_var)
  res <- map_dfr(score_names, 
                 ~run_regs_fml(
                   fml        = fml,
                   score_name = .x,
                   famhist    = famhist,
                   subset     = quote(birth_order == 1),
                   weights    = weights
                 ),
                 .id = "score_name"
  )
  res %<>% filter(grepl(":", term))
  
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
    "RLRS ~ {score_name}*({control} + 
          age_at_recruitment + I(age_at_recruitment^2))"))
  res <- lm(fml, data = famhist) %>% 
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


run_reg_fe_fertility <- function (famhist, score_names) {
  # remove scores which correlate highly with others
  # we leave EA3 and bmi_combined
  score_names <- setdiff(score_names, c("EA2_noUKB", "hip_combined", 
                                          "wc_combined", "whr_combined"))
  
  fml_scores <- paste(score_names, collapse = " + ")
  
  fml_raw <- paste("RLRS ~ ", fml_scores, " | sib_group")
  fml_mediators <- paste("RLRS ~ ", fml_scores, " + age_fte_cat | sib_group")
  
  fml_raw <- as.formula(fml_raw)
  fml_mediators <- as.formula(fml_mediators)
  
  mod_raw <- fixest::feols(fml_raw, data = famhist, note = FALSE)
  mod_mediators <- fixest::feols(fml_mediators, data = famhist, note = FALSE)
  
  tidied_raw <- broom::tidy(mod_raw, conf.int = TRUE)
  tidied_raw %<>% dplyr::filter(term %in% score_names)
  
  tidied_mediators <- broom::tidy(mod_mediators, conf.int = TRUE)
  tidied_mediators %<>% dplyr::filter(term %in% score_names)
  
  tidied_all <- bind_rows(Raw = tidied_raw, Controlled = tidied_mediators, 
                            .id = "Regression")
  
  attr(tidied_all, "n_groups") <- length(fixest::fixef(mod_raw)$sib_group)
  attr(tidied_all, "n") <- nobs(mod_raw)
    
  return(tidied_all)
}

run_cor_income <- function (famhist, score_names, age_qual_weights) {
  famhist <- famhist %>% 
               left_join(age_qual_weights, by = "f.eid") %>% 
               filter(
                 ! is.na(RLRS), 
                 ! is.na(income_cat), 
                 ! is.na(weights)
               ) %>%
               mutate(
                 child_weights = weights * n_children
               )
    
  cors_counterfactual <- purrr::map(score_names, ~{
    psych::cor.wt(famhist[c("income_cat", .x)], w = famhist$weights)
  })
  cors_actual <- purrr::map(score_names, ~{
    psych::cor.wt(famhist[c("income_cat", .x)], w = famhist$child_weights)
  })
  
  cors_counterfactual <- purrr::map_dbl(cors_counterfactual, ~ .x$r[2,1])
  cors_actual <- purrr::map_dbl(cors_actual, ~ .x$r[2,1])
  
  res <- tibble(score = score_names, actual = cors_actual, cf = cors_counterfactual)
  res <- res %>% mutate(ratio = actual/cf)
  
  res
}

run_ineq_ea3_calcs <- function (famhist, age_qual_weights, h2) {
  famhist <- left_join(famhist, age_qual_weights, by = "f.eid")
  famhist$child_weights <- famhist$n_children * famhist$weights
  
  res <- list()
  
  reg_edu_ea3 <- lm(age_fulltime_edu ~ EA3_excl_23andMe_UK, 
                    data = famhist, weights = weights)
  
  res$r2_edu_ea3 <- summary(reg_edu_ea3)$r.squared
  res$lambda <- res$r2_edu_ea3/h2
  
  reg_income_ea3 <- lm(income_cat ~ EA3_excl_23andMe_UK, 
                       data = famhist, weights = weights)
  reg_income_ea3_wt <- lm(income_cat ~ EA3_excl_23andMe_UK, 
                         data = famhist, weights = child_weights)
  
  res$r2_income_ea3 <- summary(reg_income_ea3)$r.squared
  res$r2_income_ea3_wt <- summary(reg_income_ea3_wt)$r.squared
  
  res$r2_income_true_psea <- res$r2_income_ea3 / res$lambda
  res$r2_income_true_psea_wt <- res$r2_income_ea3_wt / res$lambda
  
  fh_sibs <- famhist %>% tidyr::drop_na(sib_group, income_cat, 
                                          EA3_excl_23andMe_UK, weights, 
                                          child_weights)
  
  reg_income_ea3_fe <- fixest::feols(income_cat ~ EA3_excl_23andMe_UK | sib_group, 
                                       data = fh_sibs, 
                                       weights = fh_sibs$weights,
                                       note = FALSE
                                     )
  reg_income_ea3_fe_wt <- fixest::feols(income_cat ~ EA3_excl_23andMe_UK | sib_group, 
                                          data = fh_sibs, 
                                          weights = fh_sibs$child_weights,
                                          note = FALSE
                                        )
  
  # can we get std errors for these?
  res$r2_income_ea3_fe <- fixest::r2(reg_income_ea3_fe, type = "wr2")
  res$r2_income_ea3_fe_wt <- fixest::r2(reg_income_ea3_fe_wt, type = "wr2")
  boot_r2s <- replicate(200, {
    fh_boot <- fh_sibs %>% slice_sample(prop = 1, replace = TRUE)
    reg_fe <- fixest::feols(income_cat ~ EA3_excl_23andMe_UK | sib_group, 
                              data = fh_boot, 
                              weights = fh_boot$weights,
                              note = FALSE
                            )
    reg_fe_wt <- fixest::feols(income_cat ~ EA3_excl_23andMe_UK | sib_group, 
                                 data = fh_boot, 
                                 weights = fh_boot$child_weights,
                                 note = FALSE
                               )
    c(fixest::r2(reg_fe, type = "wr2"), fixest::r2(reg_fe_wt, type = "wr2"))
  })
  
  boot_r2_props <- boot_r2s[2,]/boot_r2s[1,] - 1
  res$r2_fe_bounds <- quantile(boot_r2_props, c(0.025, 0.975))
  
  return(res)
}


run_mediation <- function (famhist, res_all) {
  famhist <- famhist %>%
               filter(
                 ! is.na(RLRS),
                 ! is.na(age_fte_cat),
                 ! is.na(sex),
                 ! is.na(age_at_recruitment),
                 ! is.na(height),
                 ! is.na(fluid_iq),
                 ! is.na(f.2040.0.0) # risk attitude
               )
  
  sig_scores <- res_all %>% 
                  filter(dep.var == "RLRS", reg.type == "controlled") %>% 
                  filter(p.value < 0.05/33) %>% 
                  pull(score_name)
  
  run_one_mediation <- function (score_name, famhist = famhist) {
    
    controls <- "age_at_recruitment + sex + fluid_iq + height + f.2040.0.0"
    f_mediator <- as.formula(glue::glue("age_fulltime_edu ~ {score_name} +
                                           {controls}"))
    mod_mediator <- lm(f_mediator, data = famhist)
    f_y <- as.formula(glue::glue("RLRS ~ {score_name} +
                                    age_fulltime_edu +
                                    {controls}"))
    mod_y <- lm(f_y, data = famhist)
    
    # standard Baron/Kenny 1986;
    # dep_var ~ treatment + mediator + covariates
    # mediator ~ treatment + covariates
    # indirect effect = b_mediator * b_med_treat
    b_mediator   <- coef(mod_y)["age_fulltime_edu"]
    b_med_treat  <- coef(mod_mediator)[score_name]
    estimate_ind <- b_mediator * b_med_treat
                         
    # total effect = b_treatment + b_mediator * b_med_treat (from 1)
    b_treat <- coef(mod_y)[score_name]
    estimate_total <- b_treat + estimate_ind
    
    se_b_mediator <- coef(summary(mod_y))["age_fulltime_edu", "Std. Error"]
    se_b_med_treat <- coef(summary(mod_mediator))[score_name, "Std. Error"]
    
    # s.e. of indirect effect = sqrt of:
    # b_mediator^2 * var(b_treatment) + b_treatment^2 * var(b_mediator)
    # = b_mediator^2 * se(b_treatment)^2 + b_treatment^2 * se(b_mediator)^2
    var_ind <- b_mediator^2 * se_b_med_treat^2 + 
                 b_med_treat^2 * se_b_mediator^2 + 
                 se_b_med_treat^2 * se_b_mediator^2
    se_ind <- sqrt(var_ind)
    
    tibble(
      term           = score_name, 
      estimate_total = estimate_total,
      estimate_ind   = estimate_ind,
      prop_ind       = estimate_ind/estimate_total,
      se_ind         = se_ind,
      statistic_ind  = estimate_ind/se_ind,
      p_value_ind    = 2 * pnorm(abs(statistic_ind), mean = 0, 
                                   lower.tail = FALSE)
    )
  }
  
  res <- purrr::map_dfr(sig_scores, run_one_mediation, famhist = famhist)
  
  # bootstraps
  boot_res <- purrr::map_dfr(1:100, function (r) {
                              fh_boot <- famhist %>% 
                                           slice_sample(prop = 1, replace = TRUE)
                              br <- purrr::map_dfr(sig_scores, 
                                                   run_one_mediation, famhist = fh_boot)
                              br <- br[c("term", "estimate_total", "estimate_ind")]
                              br["rep"] <- r
                              br
                            })
  
  # problem: se isn't helpful because ratio not normal
  # but quantiles for 0.05/33 require 660 reps, v slow
  boot_res %<>% 
             group_by(term) %>% 
             mutate(
               prop_ind = estimate_ind/estimate_total
             ) %>%
             summarize(
               prop_ind_conf_low   = quantile(prop_ind, 0.025),
               prop_ind_conf_high  = quantile(prop_ind, 0.975)
             )
  
  res <- left_join(res, boot_res, by = "term")
  
  res
}
