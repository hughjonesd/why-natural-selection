
suppressPackageStartupMessages({
  library(readr)
  library(magrittr)
  library(dplyr)
  library(rlang)
  library(haven)
  library(sf)
  loadNamespace("santoku")
  loadNamespace("ggplot2")
  loadNamespace("car") # very annoying if it overwrites recode
  loadNamespace("matrixStats")
  loadNamespace("mlogit")
  loadNamespace("readxl")
})


join_famhist_pcs <- function (famhist, pcs) {
  left_join(famhist, pcs, by = c("f.eid" = "eid"))
}


join_famhist_resid_scores <- function (famhist, resid_scores) {
  cbind(famhist, resid_scores %>% select(-f.eid))
}


# edit_famhist <- function (famhist, score_names, ashe_income) {
#   # we get very few extra cases from adding f.2946.1.0 etc, and it makes calculating
#   # father's year of birth more complex
#   
#   names(famhist) <- sub("age_at_reqruitment", "age_at_recruitment", 
#         names(famhist))
#   
#   # remove negatives
#   famhist %<>% mutate(across(
#       c(age_fulltime_edu, starts_with(c(
#         "f.2946", "f.1845", "f.2754", "f.738",  "f.2764", "f.2405", "f.2734",
#         "f.2149", "f.1873", "f.1883", "f.2784", "f.2794", "f.709",  "f.3872",
#         "f.5057", "f.6138"
#       ))), 
#       negative_to_na
#     )
#   )
# 
#   famhist$female <- famhist$f.31.0.0 == 0
#   
#   # "Field 845 was collected from all participants except those who indicated 
#   # they have a College or University degree, as defined by their answers to 
#   # Field 6138". So, we impute this to be 21.
#   famhist$age_fulltime_edu[is.na(famhist$age_fulltime_edu) & famhist$edu_qual == 1] <- 21
#   
#   famhist$income_cat <- famhist$f.738.0.0
#   
#   # roughly speaking, these are ages in 2007-10
#   famhist$fath_age <- famhist$f.2946.0.0
#   famhist$moth_age <- famhist$f.1845.0.0
#   famhist$fath_age_birth <- famhist$fath_age - famhist$age_at_recruitment
#   famhist$moth_age_birth <- famhist$moth_age - famhist$age_at_recruitment
#   
#   # full brothers and sisters
#   famhist$nbro <- pmax(famhist$f.1873.0.0, famhist$f.1873.1.0, 
#         famhist$f.1873.2.0, na.rm = TRUE)
#   famhist$nsis <- pmax(famhist$f.1883.0.0, famhist$f.1883.1.0, 
#         famhist$f.1883.2.0, na.rm = TRUE)
#   famhist$n_sibs <- famhist$nbro + famhist$nsis + 1
#   # a few people give varying answers, we assume median is fine.
#   # including later answers picks up c. 10K extra people
#   famhist$n_older_sibs <- matrixStats::rowMedians(
#           as.matrix(famhist[, c("f.5057.0.0", "f.5057.1.0", "f.5057.2.0")]),
#           na.rm = TRUE
#         )
#   
#   famhist$n_partners <- pmax(famhist$f.2149.0.0, famhist$f.2149.1.0, 
#     famhist$f.2149.2.0, na.rm = TRUE)
#   # f.2139 is age at first sexual intercourse. -2 means "never had sex";
#   # the question about number of partners was then not asked.
#   famhist$n_partners[famhist$f.2139 == -2] <- 0
#   famhist$lo_partners <- famhist$n_partners <= 3
#   
#   famhist$n_children <- pmax(famhist$f.2405.0.0, famhist$f.2405.1.0,
#     famhist$f.2405.2.0, famhist$f.2734.0.0, famhist$f.2734.1.0, 
#     famhist$f.2734.2.0, 
#     na.rm = TRUE
#   )
#   
#   
#   famhist$n_in_household <- famhist$f.709.0.0
#   
#   famhist$with_partner   <- famhist$f.6141.0.0 == 1
#   # Many NAs, almost all from people living alone i.e. f.709 == 1
#   famhist$with_partner[famhist$n_in_household == 1] <- FALSE
#   famhist$with_partner[famhist$f.6141.0.0 == -3] <- NA
#   
#   famhist$age_fte_cat <- santoku::chop(famhist$age_fulltime_edu, 
#     c(16, 18), 
#     c("< 16", "16-18", "> 18"))
#   
#   # -7 means never went to school. We recode to 0 for simpliciy
#   famhist$edu_qual[famhist$edu_qual == -7] <- 0
#   famhist$edu_qual[famhist$edu_qual == -3] <- NA
#   
#   # we use pmax, assuming that people *can* have given birth for the first
#   # time in between surveys.
#   famhist$age_flb <- pmax(
#           famhist$f.3872.0.0, famhist$f.3872.1.0, famhist$f.3872.2.0,
#           famhist$f.2754.0.0, famhist$f.2754.1.0, famhist$f.2754.2.0,
#           na.rm = TRUE
#         )
#   famhist$age_flb_cat <- santoku::chop_equally(famhist$age_flb, 3, 
#                                                labels = lbl_discrete("-"))
#   famhist$age_llb <- pmax(
#           famhist$f.2764.0.0, famhist$f.2764.1.0, famhist$f.2764.2.0,
#           na.rm = TRUE
#         )
#   
#   famhist$year_flb <- famhist$YOB + famhist$age_flb
#   famhist$year_llb <- famhist$YOB + famhist$age_llb
#   
#   famhist$flb_cat <- santoku::fillet(famhist$age_flb, c(13, 20, 23, 26, 30, 47))
#   famhist$flb_cat %<>% forcats::fct_expand("No children")
#   famhist$flb_cat[famhist$sex == 0 &  famhist$n_children == 0] <- "No children"
#   
#   famhist$urbrur <- car::recode(famhist$f.20118.0.0,
#           "c(1, 5, 11, 12) = 'urban';
#           c(2, 6, 13, 14, 15)  = 'town';
#           c(3, 4, 7, 8, 16, 17, 18) = 'rural';
#           9 = NA_character_
#           "
#         )
#   
#   famhist[score_names] <- scale(famhist[score_names])
#   
#   famhist %<>% 
#         mutate(f.22617.0.0 = as.character(f.22617.0.0)) %>% 
#         left_join(ashe_income, by = c("f.22617.0.0" = "Code")) %>% 
#         select(-Description, -mean_pay) %>% 
#         rename(first_job_pay = median_pay) %>% 
#         mutate(first_job_pay = first_job_pay/1000)
#   
#   famhist$kids_ss <- famhist$age_at_recruitment >= 45
#   
#   return(famhist)
# }


make_famhist_long_mlogit <- function (famhist, score_names) {
  famhist$n_ch_fac <- santoku::chop(famhist$n_children, 0:5, 
                                      labels = santoku::lbl_discrete())
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
  
  reverse_coded <- c("ai_substance_use", "dpw_substance_use", 
                       "cpd_substance_use", "si_substance_use")
  rgs$rg[rgs$p2 %in% reverse_coded] <- rgs$rg[rgs$p2 %in% reverse_coded] * -1
  
  rgs
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


add_deprivation_data <- function (famhist, dep_data_dir) {
  lsoa <- sf::read_sf(file.path(dep_data_dir, "GIS"))
  sf::st_crs(lsoa) <- "EPSG:27700"
  
  famhist %<>% 
            filter(! is.na(home_lon.0.0), ! is.na(home_lat.0.0)) %>% 
            sf::st_as_sf(
              coords = c("home_lon.0.0", "home_lat.0.0"), 
              crs = "EPSG:27700" # OSGB: 1936, see 
              # https://epsg.io/27700
            )
  lsoa_codes <- sf::st_intersects(famhist, lsoa) %>% 
                  map_dbl(~ if (length(.x) == 0) NA_real_ else .x[[1]] )
  
  famhist %<>% sf::st_drop_geometry()
  famhist$lsoa_dz_cd <- lsoa$CD2011[lsoa_codes]
  
  lsoa_deprivation <- readr::read_csv(
                        file.path(dep_data_dir, "1971-to-2011-dep-density.csv"),
                        col_types = cols(
                          "Quin71" = col_factor(levels = as.character(1:5))
                        )
                      )
  lsoa_deprivation %<>% 
        select(lsoa_dz_cd, town71, Quin71) 
  famhist %<>% left_join(lsoa_deprivation, by = "lsoa_dz_cd")
  
  famhist
}
