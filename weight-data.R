
suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(haven)
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
