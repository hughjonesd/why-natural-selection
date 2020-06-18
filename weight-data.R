
suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(haven)
  library(forcats)
  library(santoku)
  library(magrittr)
})


#' @param fml A formula to be passed in to `survey::calibrate()`
#'   Possible variables: sex, age_at_recruitment, edu_qual, age_flb...
#'   These could be specified as factors, and should be for e.g. educ_level
weight_by_ghs <- function (fml, ghs_subset, famhist) {

  fml_vars <- all.vars(fml)
  famhist %<>%  filter(across(all_of(fml_vars), Negate(is.na)))
  ghs_subset %<>% filter(across(all_of(fml_vars), Negate(is.na)))
  
  pop_mx <- model.matrix(fml, data = ghs_subset)
  pop_mx <- pop_mx * ghs_subset$weight06
  pop_mx_totals <- colSums(pop_mx)
  
  design <- survey::svydesign(~1, data = famhist)
  # `rake` works on each margin individually. `postStratify` works on multiple 
  # margins, which might be better if we have enough data! `calibrate` uses a 
  # linear model (presumably to predict the population frequencies? or the
  # response probability?)
  calibrate_design <- survey::calibrate(design,
          formula = fml,
          population = pop_mx_totals,
          bounds = c(0, Inf)
        )
  calibrate_weights <- weights(calibrate_design)
  
  data.frame(f.eid = famhist$f.eid, weights = calibrate_weights)
}


weight_parents <- function (famhist, input_weights) {
  famhist <- inner_join(famhist, input_weights, by = "f.eid")
  n_sibs_are_up_to_8 <- ! is.na(famhist$n_sibs) & famhist$n_sibs <= 8
  n_sibs_up_to_8 <- famhist$n_sibs[n_sibs_are_up_to_8]
  wts <- famhist$weights[n_sibs_are_up_to_8]
  # weighted table:
  tbl <- xtabs(wts ~ n_sibs_up_to_8)
  respondent_props <- proportions(tbl)
  respondent_props <- as.vector(respondent_props)
  parent_props <- respondent_props * 1:8
  parent_props <- parent_props/sum(parent_props)
  parent_weights <- rep(NA_real_, nrow(famhist))
  
  parent_weights[n_sibs_are_up_to_8] <- parent_props[n_sibs_up_to_8]
  
  data.frame(f.eid = famhist$f.eid, weights = parent_weights)
}


weight_by_census_age_qual <- function (famhist, census_age_qual) {
  famhist$age_cat <-  santoku::chop(famhist$age_at_recruitment, 
                        c(16, 25, 35, 50, 65), 
                        lbl_discrete(" to ")) %>% 
                      fct_recode("65 and over" = "65 to 73")

  # edu_qual is:
  # 1	College or University degree
  # 2	A levels/AS levels or equivalent
  # 3	O levels/GCSEs or equivalent
  # 4	CSEs or equivalent
  # 5	NVQ or HND or HNC or equivalent
  # 6	Other professional qualifications eg: nursing, teaching
  # -7 None (but I recoded this to 0)
  # 
  # QUAL in the census data is
  # Apprenticeship: only in EW, what it says
  # Level 1: 1-4 O levels/CSE/GCSEs, NVQ level 1, foundation GNVQ; et al.
  # Level 2: 5+ O levels/CSE/GCSEs; 1 A level; or 2-3 AS levels; NVQ level 2; et al
  # Level 3: 2+ A levels, 4+ AS levels; advanced GNVQ; et al
  # Level 4 and above: degree or higher; professional qualification e.g. teaching
  # nursing accountancy; et al.
  # Other: only in EW.
  # 
  # The tricky one is "NVQ or HND or HNC". That could be anything from Level 1-4.
  # We'll assume that the majority are level 1.
  famhist$qual <- dplyr::recode(famhist$edu_qual, 
                    `0` = "No qualifications",
                    `1` = "Level 4 and above",
                    `2` = "Level 3",
                    `3` = "Level 2",
                    `4` = "Level 1",
                    `5` = "Level 1", 
                    `6` = "Level 4 and above"
                  )
  famhist <- famhist[! is.na(famhist$qual), ]
  
  fml <- ~ -1 + age_cat:qual
  pop_totals <- census_age_qual$pop
  names(pop_totals) <- sprintf("age_cat%s:qual%s", 
                         census_age_qual$age_cat, 
                         census_age_qual$qual
                       )
  
  design <- survey::svydesign(~1, probs = ~1, data = famhist)
  calibrate_design <- survey::calibrate(design,
                        formula = fml,
                        population = pop_totals,
                        bounds = c(0, Inf)
                      )
  weights <- weights(calibrate_design)
  
  # passes back only non-NA weights. So use left_join and not cbind
  data.frame(f.eid = famhist$f.eid, weights = weights)
}


weight_by_census_msoa <- function (famhist, census_msoa, famhist_msoa) {
  
  famhist <- inner_join(famhist, famhist_msoa, by = "f.eid")
  
  famhist$age_cat <-  chop(famhist$age_at_recruitment, 
                        c(40, 45, 50, 55, 60, 65, 70, 75), 
                        extend = FALSE,
                        labels = lbl_discrete(" to ")
                      )
  
  famhist$sex <- as.factor(famhist$sex)

  msoa_table <- census_msoa %>% 
        rename(msoa = msoa_code, age_cat = age, Freq = pop) %>%
        mutate(sex = as.factor(sex)) %>% 
        semi_join(famhist, by = c("age_cat", "sex", "msoa")) %>% 
        mutate(
          category = sprintf("%s/%s/%s", age_cat, sex, msoa)
        ) %>% 
        # otherwise we get ONE person whose category isn't in the census data:
        filter(Freq > 0) 
        
  famhist %<>% 
    tidyr::drop_na(age_cat, sex, msoa) %>% 
    semi_join(msoa_table, by = c("age_cat", "sex", "msoa")) %>% 
    mutate(
      category = sprintf("%s/%s/%s", age_cat, sex, msoa)
    )
  
  design <- survey::svydesign(~1, probs = ~1, data = famhist)
  msoa_table %<>% select(category, Freq)
  raked <- survey::rake(design, list(~category), list(msoa_table))

  data.frame(f.eid = famhist$f.eid, weights = weights(raked))
}
