
suppressPackageStartupMessages({
  library(readr)
  library(magrittr)
  library(dplyr)
  library(rlang)
  library(santoku)
  library(haven)
  loadNamespace("car") # very annoying if it overwrites recode
  loadNamespace("matrixStats")
  loadNamespace("mlogit")
})

# utility function:
negative_to_na <- function (x) {
  x[x < 0] <- NA
  x
}

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
  
  resid_scores$dummy <- NULL
  resid_scores
}


make_famhist <- function (
        ph_file, 
        famhist_file,
        famhist2_file,
        famhist3_file,
        famhist4_file,
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
  famhist4 <- read_csv(famhist4_file, col_types = strrep("d", 33))
  famhist5 <- read_csv(famhist5_file, col_types = strrep("d", 4))
  names(famhist) <- paste0("f.", names(famhist))
  names(famhist) <- gsub("\\-", ".", names(famhist))
  names(famhist2) <- paste0("f.", names(famhist2))
  names(famhist2) <- gsub("\\-", ".", names(famhist2))
  names(famhist3) <- paste0("f.", names(famhist3))
  names(famhist3) <- gsub("\\-", ".", names(famhist3))
  names(famhist4) <- paste0("f.", names(famhist4))
  names(famhist4) <- gsub("\\-", ".", names(famhist4))
  names(famhist5) <- paste0("f.", names(famhist5))
  names(famhist5) <- gsub("\\-", ".", names(famhist5))
  
  famhist %<>% left_join(ph, by = c("f.eid" = "eid"))
  famhist %<>% left_join(famhist2, by = "f.eid")
  famhist %<>% left_join(famhist3, by = "f.eid")
  famhist %<>% left_join(famhist4, by = "f.eid")
  famhist %<>% left_join(famhist5, by = "f.eid")
  
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
  
  names(famhist) <- sub("age_at_reqruitment", "age_at_recruitment", 
        names(famhist))
  
  # remove negatives
  famhist %<>% mutate(across(
      c(age_fulltime_edu, starts_with(c(
        "f.2946", "f.1845", "f.2754", "f.738",  "f.2764", "f.2405", "f.2734",
        "f.2149", "f.1873", "f.1883", "f.2784", "f.2794", "f.709",  "f.3872",
        "f.5057"
      ))), 
      negative_to_na
    )
  )
  
  famhist$income_cat <- famhist$f.738.0.0
  
  # roughly speaking, these are ages in 2007-10
  famhist$fath_age <- famhist$f.2946.0.0
  famhist$moth_age <- famhist$f.1845.0.0
  famhist$fath_age_birth <- famhist$fath_age - famhist$age_at_recruitment
  famhist$moth_age_birth <- famhist$moth_age - famhist$age_at_recruitment
  
  # full brothers and sisters
  famhist$nbro <- pmax(famhist$f.1873.0.0, famhist$f.1873.1.0, 
        famhist$f.1873.2.0, na.rm = TRUE)
  famhist$nsis <- pmax(famhist$f.1883.0.0, famhist$f.1883.1.0, 
        famhist$f.1883.2.0, na.rm = TRUE)
  famhist$n_sibs <- famhist$nbro + famhist$nsis + 1
  # a few people give varying answers, we assume median is fine.
  # including later answers picks up c. 10K extra people
  famhist$n_older_sibs <- matrixStats::rowMedians(
          as.matrix(famhist[, c("f.5057.0.0", "f.5057.1.0", "f.5057.2.0")]),
          na.rm = TRUE
        )
  
  famhist$n_partners <- pmax(famhist$f.2149.0.0, famhist$f.2149.1.0, 
    famhist$f.2149.2.0, na.rm = TRUE)
  famhist$lo_partners <- famhist$n_partners <= 3
  
  famhist$n_children <- pmax(famhist$f.2405.0.0, famhist$f.2405.1.0,
    famhist$f.2405.2.0, famhist$f.2734.0.0, famhist$f.2734.1.0, 
    famhist$f.2734.2.0, 
    na.rm = TRUE
  )
  
  
  famhist$n_in_household <- famhist$f.709.0.0
  
  famhist$with_partner   <- famhist$f.6141.0.0 == 1
  # Many NAs, almost all from people living alone i.e. f.709 == 1
  famhist$with_partner[famhist$n_in_household == 1] <- FALSE
  famhist$with_partner[famhist$f.6141.0.0 == -3] <- NA
  
  famhist$age_fte_cat <- santoku::chop(famhist$age_fulltime_edu, 
    c(16, 18), 
    c("< 16", "16-18", "> 18"))
  # -7 means never went to school. We recode to 0 for simpliciy
  famhist$edu_qual[famhist$edu_qual == -7] <- 0
  famhist$edu_qual[famhist$edu_qual == -3] <- NA
  
  # we use pmax, assuming that people *can* have given birth for the first
  # time in between surveys.
  famhist$age_flb <- pmax(
          famhist$f.3872.0.0, famhist$f.3872.1.0, famhist$f.3872.2.0,
          famhist$f.2754.0.0, famhist$f.2754.1.0, famhist$f.2754.2.0,
          na.rm = TRUE
        )
  famhist$age_llb <- pmax(
          famhist$f.2764.0.0, famhist$f.2764.1.0, famhist$f.2764.2.0,
          na.rm = TRUE
        )
  
  famhist$year_flb <- famhist$YOB + famhist$age_flb
  famhist$year_llb <- famhist$YOB + famhist$age_llb
  
  famhist$flb_cat <- santoku::fillet(famhist$age_flb, c(13, 20, 23, 26, 30, 47))
  famhist$flb_cat %<>% forcats::fct_expand("No children")
  famhist$flb_cat[famhist$sex == 0 &  famhist$n_children == 0] <- "No children"
  
  famhist$urbrur <- car::recode(famhist$f.20118.0.0,
          "c(1, 5, 11, 12) = 'urban';
          c(2, 6, 13, 14, 15)  = 'town';
          c(3, 4, 7, 8, 16, 17, 18) = 'rural';
          9 = NA_character_
          "
        )
  
  famhist[score_names] <- scale(famhist[score_names])
  
  return(famhist)
}


make_famhist_long_mlogit <- function (famhist, score_names) {
  famhist$n_ch_fac <- santoku::chop(famhist$n_children, 
        brk_left(0:5, close_end = FALSE), lbl_integer())
  fh_subset <- famhist %>% select(n_ch_fac, all_of(score_names))
  mlogit::mlogit.data(fh_subset, choice = "n_ch_fac", shape = "wide", 
        alt.levels = levels(fh_subset$n_ch_fac))
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


make_ghs_subset <- function(ghs_file) {
  # calculate population totals from ghs file for variables of interest
  # using only 40-70 year olds I guess!
  
  ghs <- haven::read_dta(ghs_file)
  # variables of interest: 
  # sex (1 male)
  # age
  # ethnic = 1 = white nbritish
  # EdAge - age left FTE
  # hiqual - highest qualification 
  # edlev00 - education level
  # edlev10 - ditto, fewer categories 
  #   1 = a level or above; 2 = o level; 3 other ; 4 none. -
  #   -9 does not apply - including if age > 59
  # weight06 - "weight you should use to weight the data"
  # 
  # grfam1h - gross weekly income of family (harmonized)
  # grhhold - gross weekly household income
  # for comparison, the UKBB question was "what is the average total
  # income before tax received by your household" (p.a.)
  # 
  # ten1 - tenure
  # llord - landlord
  
  ghs_subset <- ghs %>% 
    zap_labels() %>% 
    filter(
      ethnic == 1, sex %in% 1:2, edage < age, age %in% c(40:71),
      chbnbm1 < age
    ) %>% 
    select(sex, age, edage, weight06, ten1, llord, chbnbm1, edlev00, chnbrnt) 
  
  ghs_subset$edlev00[ghs_subset$edlev00 == -9] <- 0 # "never went to school"
  ghs_subset %<>% mutate(across(
    c(edage, chbnbm1, ten1, llord, chbnbm1, edlev00, chnbrnt), 
    negative_to_na
  ))
  
  # create variables with the same meaning as in famhist
  # ghs: 1 male, 2 female; famhist: 0 female, 1 male
  ghs_subset$sex                <- 2 - ghs_subset$sex
  ghs_subset$age_at_recruitment <- ghs_subset$age
  ghs_subset$age_fulltime_edu   <- ghs_subset$edage
  ghs_subset$age_flb            <- ghs_subset$chbnbm1

  # YearsEdu is 7, 10, 13, 15, 19 or 20. It maps from
  # edu_qual, which is:
  # 1	College or University degree
  # 2	A levels/AS levels or equivalent
  # 3	O levels/GCSEs or equivalent
  # 4	CSEs or equivalent
  # 5	NVQ or HND or HNC or equivalent
  # 6	Other professional qualifications eg: nursing, teaching
  # -7 None (but I recoded this to 0)
  # 
  # ghs$edlev00 is:
  # -9    Never attended school
  # -8    NA
  # -6    CHILD/OUT AGE/NO INT
  # 1     Higher Degree
  # 2    First Degree
  # 3    Teaching qualification
  # 4    Other higher qualification
  # 5    Nursing qualification
  # 6    GCE A level in two or more subjects
  # 7    GCE A level in one subject
  # 8    GCSE/Olevel, standard grades, 5+
  # 9    GCSE/Olevel 1-4
  # 10    CSE below grade 1, GCSE below grade C
  # 11    Apprenticeship
  # 12    Other qualification
  # 13    no qualification
  ghs_subset$edu_qual     <- dplyr::recode(ghs_subset$edlev00,
    "0" = 0,
    "1"  = 1,
    "2"  = 1,
    "3"  = 6,
    "4"  = 5,
    "5"  = 6,
    "6"  = 2,
    "7"  = 2,
    "8"  = 3,
    "9"  = 3,
    "10"  = 4,
    "11"  = 5,
    "12"  = 0, # mostly "started an apprenticeship, not yet finished"
    "13"  = 0
  )
  
  # our calibration model for women will be: 
  # age_at_recruitment, edu_qual, and categories including:
  # never had a child + age_flb quantiles
  # values below are 
  ghs_subset$flb_cat <- santoku::fillet(ghs_subset$age_flb, c(13, 20, 23, 26, 30, 47))
  ghs_subset$flb_cat %<>% forcats::fct_expand("No children")
  ghs_subset$flb_cat[ghs_subset$sex == 0 & ghs_subset$chnbrnt == 0] <- "No children"
  
  ghs_subset
}

