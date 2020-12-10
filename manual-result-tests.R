
# manual tests for each result

library(drake)
library(dplyr)
library(testthat)

loadd(score_names)
loadd(famhist)

lm_coef <- function (fml, target = 2, data = famhist %>% filter(kids_ss)) {
  mod <- lm(fml, data = data)
  coef(mod)[target]
}

estimate_for <- function (data, score_name, term) {
  if (! missing(term)) data$estimate[data$term == term] else 
        data$estimate[data$score_name == score_name]
}

test_that("res_all", {
  loadd(res_all)
  loadd(resid_scores)
  
  expect_equivalent(
    res_all %>% 
      filter(dep.var == "N children", reg.type == "raw") %>% 
      estimate_for("ADHD_2017"),
    lm_coef(n_children ~ ADHD_2017_raw)
  )
  
  expect_equivalent(
    res_all %>% 
      filter(dep.var == "N siblings", reg.type == "raw") %>% 
      estimate_for("alcohol_schumann"),
    lm_coef(n_sibs ~ alcohol_schumann_raw, data = famhist)
  )
  
  
  fh2 <- cbind(famhist, resid_scores[-1])
  expect_equivalent(
    res_all %>% 
      filter(dep.var == "N children", reg.type == "controlled") %>% 
      estimate_for("height_combined"),
    lm_coef(n_children ~ height_combined, data = fh2 %>% filter(kids_ss))
  )
  
  expect_equivalent(
    res_all %>% 
      filter(dep.var == "N siblings", reg.type == "controlled") %>% 
      estimate_for("height_combined"),
    lm_coef(n_sibs ~ height_combined, data = fh2)
  )
})


test_that("res_pcs", {
  loadd(res_pcs)
  loadd(pcs)
  fh <- left_join(famhist, pcs, by = c("f.eid" = "eid"))
  expect_equivalent(
    res_pcs %>% filter(dep.var == "N children", term == "PC33") %>% pull(estimate),
    lm_coef(n_children ~ PC33, data = fh %>% filter(kids_ss))
  )
  
  expect_equivalent(
    res_pcs %>% filter(dep.var == "N siblings", term == "PC45") %>% pull(estimate),
    lm_coef(n_sibs ~ PC45, data = fh)
  )
})


test_that("res_wt_flb_weights", {
  loadd(res_wt_flb_weights)
  loadd(flb_weights)
  fh <- left_join(famhist, flb_weights, by = "f.eid")
  mod <- lm(n_children ~ age_at_menarche, data = fh, weights = weights, subset = kids_ss)
  coef <- coef(mod)[2]
  expect_equivalent(
    coef,
    res_wt_flb_weights %>% estimate_for(term = "age_at_menarche")
  )
})


test_that("res_wt_age_qual_weights", {
  loadd(res_wt_age_qual_weights)
  loadd(age_qual_weights)
  fh <- left_join(famhist, age_qual_weights, by = "f.eid")
  mod <- lm(n_children ~ age_at_menarche, data = fh, weights = weights, subset = kids_ss)
  coef <- coef(mod)[2]
  expect_equivalent(
    coef,
    res_wt_age_qual_weights %>% estimate_for(term = "age_at_menarche")
  )
})


test_that("res_wt_msoa_weights", {
  loadd(res_wt_msoa_weights)
  loadd(msoa_weights)
  fh <- left_join(famhist, msoa_weights, by = "f.eid")
  mod <- lm(n_children ~ bipolar, data = fh, weights = weights, subset = kids_ss)
  coef <- coef(mod)[2]
  expect_equivalent(
    coef,
    res_wt_msoa_weights %>% estimate_for(term = "bipolar")
  )
})


test_that("res_sibs_parent_weights", {
  loadd(res_sibs_parent_weights)
  loadd(parent_weights)
  fh <- left_join(famhist, parent_weights, by = "f.eid")
  mod <- lm(n_sibs ~ age_at_menarche, data = fh, weights = weights)
  coef <- coef(mod)[2]
  expect_equivalent(
    coef,
    res_sibs_parent_weights %>% estimate_for(term = "age_at_menarche")
  )
})


test_that("res_children_comparison", {
  loadd(res_children_comparison)
  loadd(age_qual_weights)
  fh <- left_join(famhist, age_qual_weights, by = "f.eid")
  mod <- lm(n_children ~ age_at_menarche, data = fh, weights = weights, subset = kids_ss & n_children > 0)
  coef <- coef(mod)[2]
  expect_equivalent(
    coef,
    res_children_comparison %>% estimate_for(term = "age_at_menarche")
  )
})


test_that("res_period_children", {
  loadd(res_period_children)
  loadd(age_qual_weights)
  fh <- left_join(famhist, age_qual_weights, by = "f.eid")
  mod <- lm(n_children ~ openness*I(YOB >= 1950), fh, weights = weights, subset = kids_ss)
  coef1 <- coef(mod)["openness"] # early
  coef2 <- coef(mod)[4] + coef1 # late
  expect_equivalent(
    coef1,
    res_period_children %>% estimate_for(term = "year_splitearly:openness")
  )
  expect_equivalent(
    coef2,
    res_period_children %>% estimate_for(term = "year_splitlate:openness")
  )
})



test_that("res_period_parents", {
  loadd(res_period_parents)
  loadd(parent_weights)
  fh <- left_join(famhist, parent_weights, by = "f.eid")
  mod <- lm(n_sibs ~ openness*I(YOB >= 1950), fh, weights = weights)
  coef1 <- coef(mod)["openness"] # early
  coef2 <- coef(mod)[4] + coef1 # late
  expect_equivalent(
    coef1,
    res_period_parents %>% estimate_for(term = "year_splitearly:openness")
  )
  expect_equivalent(
    coef2,
    res_period_parents %>% estimate_for(term = "year_splitlate:openness")
  )
})


test_that("res_sex", {
  loadd(res_sex)
  coef1 <- lm_coef(n_children ~ cannabis*factor(sex), "cannabis")
  coef2 <- coef1 + lm_coef(n_children ~ cannabis*factor(sex), "cannabis:factor(sex)1")
  expect_equivalent(
    coef1,
    res_sex %>% filter(term == "cannabis", sex == "Female") %>% pull(estimate)
  )
  expect_equivalent(
    coef2,
    res_sex %>% filter(term == "cannabis", sex == "Male") %>% pull(estimate)
  )
})


test_that("res_age_flb", {
  loadd(res_age_flb)
  expect_equivalent(
    lm_coef(n_children ~ openness + age_flb),
    res_age_flb %>% filter(score_name == "openness", term == "openness") %>% pull(estimate)
  )
})


test_that("res_age_flb_cross", {
  loadd(res_age_flb_cross)
  coefs <- lm_coef(n_children ~ age_flb_cat + neuroticism:age_flb_cat, 4:6)
  expect_equivalent(
    coefs,
    res_age_flb_cross %>% estimate_for("neuroticism")
  )
})


test_that("res_age_flb_dv", {
  loadd(res_age_flb_dv)
  coef <- lm_coef(age_flb ~ hip_combined)
  expect_equivalent(
    coef,
    res_age_flb_dv %>% estimate_for(term = "hip_combined")
  )
})


test_that("res_age_birth_parents", {
  loadd(res_age_birth_parents)
  coef <- lm_coef(n_sibs ~ hip_combined + fath_age_birth, 
                  data = famhist %>% filter(birth_order == 1))
  expect_equivalent(
    coef,
    res_age_birth_parents %>% 
      filter(control == "fath_age_birth", term == "hip_combined") %>% 
      pull(estimate)
  )
  
  coef <- lm_coef(n_sibs ~ hip_combined + moth_age_birth, 
                  data = famhist %>% filter(birth_order == 1))
  expect_equivalent(
    coef,
    res_age_birth_parents %>% 
      filter(control == "moth_age_birth", term == "hip_combined") %>% 
      pull(estimate)
  )
})


test_that("res_age_birth_parents_dv", {
  loadd(res_age_birth_parents_dv)
  coef <- lm_coef(fath_age_birth ~ age_at_menopauze, 
                  data = famhist %>% filter(birth_order == 1))
  expect_equivalent(
    coef, 
    res_age_birth_parents_dv %>% 
      filter(dep.var == "fath_age_birth", term == "age_at_menopauze") %>% 
      pull(estimate)
  )
})


test_that("res_partners", {
  loadd(res_partners)
  coef <- lm_coef(n_children ~ age_at_menopauze*lo_partners, c(2, 4),
                  data = famhist %>% filter(kids_ss, sex == 0))
  expect_equivalent(
    coef[1],
    res_partners %>% 
      filter(sex == "Female") %>% 
      estimate_for(term = "lo_partnersFALSE:age_at_menopauze")
  )
  expect_equivalent(
    coef[1] + coef[2],
    res_partners %>% 
      filter(sex == "Female") %>% 
      estimate_for(term = "lo_partnersTRUE:age_at_menopauze")
  )
  
  coef <- lm_coef(n_children ~ age_at_menopauze*lo_partners, c(2, 4),
                  data = famhist %>% filter(kids_ss, sex == 1))
  expect_equivalent(
    coef[1],
    res_partners %>% 
      filter(sex == "Male") %>% 
      estimate_for(term = "lo_partnersFALSE:age_at_menopauze")
  )
  expect_equivalent(
    coef[1] + coef[2],
    res_partners %>% 
      filter(sex == "Male") %>% 
      estimate_for(term = "lo_partnersTRUE:age_at_menopauze")
  )
})



test_that("res_with_partner_sex", {
  loadd(res_with_partner_sex)
  coef <- lm_coef(n_children ~ age_at_menopauze*with_partner, c(2, 4),
                  data = famhist %>% filter(kids_ss, sex == 0))
  expect_equivalent(
    coef[1],
    res_with_partner_sex %>% 
      filter(sex == "Female") %>% 
      estimate_for(term = "with_partnerFALSE:age_at_menopauze")
  )
  expect_equivalent(
    coef[1] + coef[2],
    res_with_partner_sex %>% 
      filter(sex == "Female") %>% 
      estimate_for(term = "with_partnerTRUE:age_at_menopauze")
  )
  
  coef <- lm_coef(n_children ~ age_at_menopauze*with_partner, c(2, 4),
                  data = famhist %>% filter(kids_ss, sex == 1))
  expect_equivalent(
    coef[1],
    res_with_partner_sex %>% 
      filter(sex == "Male") %>% 
      estimate_for(term = "with_partnerFALSE:age_at_menopauze")
  )
  expect_equivalent(
    coef[1] + coef[2],
    res_with_partner_sex %>% 
      filter(sex == "Male") %>% 
      estimate_for(term = "with_partnerTRUE:age_at_menopauze")
  )
})


test_that("res_edu", {
  loadd(res_edu)
  coef <- lm_coef(n_children ~ age_fte_cat + sc_substance_use:age_fte_cat, 4:6)
  expect_equivalent(
    coef,
    res_edu %>% estimate_for(term = "sc_substance_use")
  )
})


test_that("res_income", {
  loadd(res_income)
  coef <- lm_coef(n_children ~ factor(income_cat) + alzheimer:factor(income_cat), 6:10)
  expect_equivalent(
    coef,
    res_income %>% estimate_for(term = "alzheimer")
  )
})



test_that("res_income_controlled", {
  loadd(res_income_controlled)
  av <- car::Anova(lm(n_children ~ cpd_substance_use*(factor(income_cat) + 
                     age_at_recruitment + I(age_at_recruitment^2)), 
                     data = famhist,
                     subset = kids_ss
                  ))
  
  expect_equivalent(
    av %>% broom::tidy() %>% filter(term == "cpd_substance_use:factor(income_cat)") %>% pull(p.value),
    res_income_controlled %>% filter(term == "cpd_substance_use:factor(income_cat)") %>% pull(p.value)
  )
})


test_that("res_edu_controlled", {
  loadd(res_edu_controlled)
  av <- car::Anova(lm(n_children ~ cpd_substance_use*(factor(age_fte_cat) + 
                     age_at_recruitment + I(age_at_recruitment^2)), 
                     data = famhist,
                     subset = kids_ss
                  ))
  
  expect_equivalent(
    av %>% broom::tidy() %>% filter(term == "cpd_substance_use:factor(age_fte_cat)") %>% pull(p.value),
    res_edu_controlled %>% filter(term == "cpd_substance_use:factor(age_fte_cat)") %>% pull(p.value)
  )
  
})


test_that("res_partners_controlled", {
  loadd(res_partners_controlled)
  coefs <- lm_coef(n_children ~ age_at_recruitment + I(age_at_recruitment^2) + lo_partners +
              SCZ2:lo_partners + SCZ2:(age_at_recruitment + I(age_at_recruitment^2)), 
              target = c("lo_partnersFALSE:SCZ2", "lo_partnersTRUE:SCZ2"),
              data   = famhist %>% filter(kids_ss & sex == 0)
            )
  expect_equivalent(coefs[1],
                    res_partners_controlled %>% 
                      filter(sex == "Female") %>% 
                      estimate_for(term = "lo_partnersFALSE:SCZ2")
                   )
  expect_equivalent(coefs[2],
                    res_partners_controlled %>% 
                      filter(sex == "Female") %>% 
                      estimate_for(term = "lo_partnersTRUE:SCZ2")
                   )
  
  coefs <- lm_coef(n_children ~ age_at_recruitment + I(age_at_recruitment^2) + lo_partners +
              SCZ2:lo_partners + SCZ2:(age_at_recruitment + I(age_at_recruitment^2)), 
              target = c("lo_partnersFALSE:SCZ2", "lo_partnersTRUE:SCZ2"),
              data   = famhist %>% filter(kids_ss & sex == 1)
            )
  expect_equivalent(coefs[1],
                    res_partners_controlled %>% 
                      filter(sex == "Male") %>% 
                      estimate_for(term = "lo_partnersFALSE:SCZ2")
                   )
  expect_equivalent(coefs[2],
                    res_partners_controlled %>% 
                      filter(sex == "Male") %>% 
                      estimate_for(term = "lo_partnersTRUE:SCZ2")
                   )
})

test_that("res_ee_control", {
  loadd(res_ee_control)
  loadd(ashe_income)
  source("~/import-ukbb-data/import-ukbb-data.R")
  fh <- add_ashe_income(famhist, ashe_income)
  fh %<>% filter(kids_ss)
  coef <- lm_coef(n_children ~ dpw_substance_use + age_fte_cat + first_job_pay, 
                    data = fh)
         
  expect_equivalent(
    coef, 
    res_ee_control %>% estimate_for(term = "dpw_substance_use")
  )          
})

