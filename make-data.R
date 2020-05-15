
suppressPackageStartupMessages({
  library(readr)
  library(magrittr)
  library(dplyr)
  library(rlang)
  library(santoku)
})


make_pcs <- function (pcs_file) {
  pcs <- read_table(pcs_file)
  pc_names <- grep("PC", names(pcs), value = TRUE)
  pcs[pc_names] <- scale(pcs[pc_names])
  pcs
}


join_famhist_pcs <- function (famhist, pcs) {
  left_join(famhist, pcs, by = c("f.eid" = "IID"))
}


join_famhist_resid_scores <- function (famhist, resid_scores) {
  cbind(famhist, resid_scores)
}


make_resid_scores <- function (famhist, pcs, score_names) {
  famhist <- join_famhist_pcs(famhist, pcs)
  resid_scores <- data.frame(dummy = numeric(nrow(famhist))) 
  
  for (score_name in score_names) {
    resid_fml <- paste(score_name, "~", paste0("PC", 1:40, collapse = " + "))
    resid_score <- resid(lm(as.formula(resid_fml), famhist, 
          na.action = na.exclude))
    resid_scores[[paste0(score_name, "_resid")]] <- resid_score
  }
  
  resid_score$dummy <- NULL
  resid_scores
}


make_famhist <- function (
        ph_file, 
        famhist_file,
        famhist2_file,
        famhist3_file,
        pgs_dir
      ) {
  ph <- read_csv(ph_file, col_types = cols(
    geno_measurement_plate = col_skip(),
    geno_measurement_well = col_skip(),
    .default = col_double()
  ))
  
  famhist <- read_csv(famhist_file, col_types = strrep("d", 40))
  famhist2 <- read_csv(famhist2_file, col_types = strrep("d", 17))
  famhist3 <- read_csv(famhist3_file)
  
  names(famhist) <- paste0("f.", names(famhist))
  names(famhist) <- gsub("\\-", ".", names(famhist))
  names(famhist2) <- paste0("f.", names(famhist2))
  names(famhist2) <- gsub("\\-", ".", names(famhist2))
  names(famhist3) <- paste0("f.", names(famhist3))
  names(famhist3) <- gsub("\\-", ".", names(famhist3))
  
  famhist %<>% left_join(ph, by = c("f.eid" = "eid"))
  famhist %<>% left_join(famhist2, by = "f.eid")
  famhist %<>% left_join(famhist3, by = "f.eid")
  
  # only "genetic" whites
  famhist %<>% filter(! is.na(genetic_ethnic_grouping))
  
  for (pgs_file in list.files(pgs_dir, pattern = "csv$", full.names = TRUE)) {
    score_name <- sub(".*UKB\\.AMC\\.(.*?)\\..*", "\\1", pgs_file, perl = TRUE)
    pgs <- read_delim(pgs_file, delim = " ", col_types = "dd")
    pgs %<>% filter(FID > 0) 
    names(pgs)[2] <- score_name # instead of "SCORE"
    famhist %<>% left_join(pgs, by = c("f.eid" = "FID"))
  }
  
  return(famhist)
}


edit_famhist <- function (famhist, score_names) {
  # we get very few extra cases from adding f.2946.1.0 etc, and it makes calculating
  # father's year of birth more complex
  
  # roughly speaking, these are ages in 2007-10
  famhist$fath_age <- famhist$f.2946.0.0
  famhist$fath_age[famhist$fath_age < 0] <- NA
  famhist$moth_age <- famhist$f.1845.0.0
  famhist$moth_age[famhist$moth_age < 0] <- NA
  
  # full brothers and sisters
  famhist$nbro <- pmax(famhist$f.1873.0.0, famhist$f.1873.1.0, 
    famhist$f.1873.2.0, na.rm = TRUE)
  famhist$nbro[famhist$nbro < 0] <- NA
  famhist$nsis <- pmax(famhist$f.1883.0.0, famhist$f.1883.1.0, 
    famhist$f.1883.2.0, na.rm = TRUE)
  famhist$nsis[famhist$nsis < 0] <- NA
  famhist$n_sibs <- famhist$nbro + famhist$nsis + 1
  
  famhist$n_partners <- pmax(famhist$f.2149.0.0, famhist$f.2149.1.0, 
    famhist$f.2149.2.0, na.rm = TRUE)
  famhist$n_partners[famhist$n_partners < 0] <- NA
  famhist$lo_partners <- famhist$n_partners <= 3
  
  famhist$n_children <- pmax(famhist$f.2405.0.0, famhist$f.2405.1.0,
    famhist$f.2405.2.0, famhist$f.2734.0.0, famhist$f.2734.1.0, 
    famhist$f.2734.2.0, 
    na.rm = TRUE
  )
  famhist$n_children[famhist$n_children < 0] <- NA
  
  famhist$age_fulltime_edu[famhist$age_fulltime_edu < 0] <- NA
  famhist$age_fte_cat <- santoku::chop(famhist$age_fulltime_edu, 
    c(16, 18), 
    c("< 16", "16-18", "> 18"))
  child_vars <- rlang::quos(f.2754.0.0, f.2754.1.0, f.2754.2.0)
  famhist %<>% mutate_at(vars(!!!child_vars), ~ ifelse(.x < 0, NA, .x))
  famhist$age_flb <- rowMeans(famhist %>% select(!!!child_vars), na.rm = TRUE)
  famhist$year_flb <- famhist$YOB + famhist$age_flb
  
  famhist$income_cat <- famhist$f.738.0.0
  famhist$income_cat[famhist$income_cat < 0] <- NA
  
  famhist[score_names] <- scale(famhist[score_names])
  
  return(famhist)
}


make_rgs <- function (rgs_file) {
  rgs <- read_csv(rgs_file)
  
  rgs$p2 <- rgs$p2 %>% dplyr::recode(
    EA3_excl_23andMe_and_allUK = "EA3_excl_23andMe_UK",
    cognitve_ability.noUKB     = "cognitive_ability",
    autism_2017.ipsych.pgc     = "autism_2017",
  ) %>% {sub("\\.GPC\\.23andme$", "", .)}
  
  rgs
}


make_mf_pairs <- function (mf_pairs_file, famhist) {
  mf_pairs <- read_csv(mf_pairs_file, col_types = "dddccccc")
  famhist_tmp <- famhist %>% dplyr::select(! (starts_with("f.") | 
      starts_with("PC") | starts_with("home") |
      starts_with("assessment") | starts_with("birth_")), f.eid)
  mf_pairs %<>% 
    left_join(famhist_tmp, by = c("ID.m" = "f.eid")) %>% 
    left_join(famhist_tmp, by = c("ID.f" = "f.eid"), suffix = c(".m", ".f"))
  
  mf_pairs
}


weight_by_ghs <- function (ghs_file, famhist) {
  library(haven)
  
  ghs <- haven::read_dta(ghs_file)
  # variables of interest: 
  # sex (1 male)
  # age
  # ethnic = 1 = white nbritish
  # EdAge - age left FTE
  # hiqual - highest qualification 
  # edlev00 - education level
  # edlev10 - ditto, fewer categories
  # weight06 - "weight you should use to weight the data"
  # 
  # grfam1h - gross weekly income of family (harmonized)
  # grhhold - gross weekly household income
  # for comparison, the UKBB question was "what is the average total
  # income before tax received by your household" (p.a.)
  # 
  # ten1 - tenure
  # llord - landlord
  
  ghs <- ghs %>% 
    zap_labels() %>% 
    filter(ethnic == 1, sex %in% 1:2, edage > 0, edage < age) %>% 
    select(sex, age, edage, weight06, ten1, llord) 
  
  ghs_props <- ghs %>% 
    group_by(sex, age, edage) %>% 
    summarize(
      weighted_n = sum(weight06)
    )
  ghs_prop_totals <- ghs %>% 
    group_by(sex, age) %>% 
    summarize(
      weighted_n_sex_age = sum(weight06)
    )
  ghs_props <- ghs_props %>% 
    left_join(ghs_prop_totals, by = c("sex", "age")) %>% 
    mutate(
      prop_age_fte = weighted_n/weighted_n_sex_age,
      sex = recode(sex, "1" = "male", "2" = "female")
    ) %>% 
    select(-weighted_n, -weighted_n_sex_age)
  
  ukbb_props <- famhist %>% 
    group_by(sex, age_at_reqruitment, age_fulltime_edu) %>% 
    summarize(n = n())
  ukbb_prop_totals <- famhist %>% 
    group_by(sex, age_at_reqruitment) %>% 
    summarize(n_sex_age = n())
  ukbb_props <- ukbb_props %>% 
    left_join(ukbb_prop_totals, by = c("sex", "age_at_reqruitment")) %>% 
    mutate(
      prop_age_fte = n/n_sex_age,
      sex = recode(sex, "0" = "female", "1" = "male")
    ) %>% 
    select(-n, -n_sex_age)
  
  ghs_weight_df <- left_join(ukbb_props, ghs_props, 
    by = c("sex", "age_at_reqruitment" = "age", "age_fulltime_edu" = "edage"),
    suffix = c(".ukbb", ".ghs")) %>% 
    mutate(
      educ_weight = prop_age_fte.ghs/prop_age_fte.ukbb
    ) %>% 
    select(-starts_with("prop")) %>% 
    filter(
      ! is.na(educ_weight)
    ) %>% mutate(
      sex = ifelse(sex == "female", 0, 1)
    )
  
  famhist %>% left_join(ghs_weight_df) %>% select(f.eid, educ_weight)
}
