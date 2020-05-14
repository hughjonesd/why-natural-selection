
library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(santoku)


add_educ_weights <- function (famhist_nowt, nomis_stats) {
  nomis_stats %<>% 
        slice(1) %>%  
        pivot_longer(
          cols      = c(-date,-geography, -`geography code`), 
          names_to  = c("age", "qual", "ethn", "measure"), 
          names_sep = ";"
        ) %>% 
        mutate(
          across(c(age, qual, ethn), ~sub(".*: ", "", .x))
        ) %>% 
        filter(ethn == "English/Welsh/Scottish/Northern Irish/British") %>% 
        select(-starts_with("geog"), -ethn, -measure)
  
  census_props <- nomis_stats %>% 
        filter(
          qual != "Highest level of qualification", 
          age != "Age 16 and over"
        ) %>% group_by(age) %>% 
        mutate(percent = value/sum(value)) %>% 
        select(-date, -value)
  
  
  famhist_nowt %<>% mutate(
    age = santoku::chop(
        age_at_reqruitment, 
        breaks = c(35, 50, 65, 75, 100),
        labels = c("Age 35 to 49", "Age 50 to 64", "Age 65 to 74", "Age 75 and over")
      ),
    qual = dplyr::recode(edu_qual,
        "1" = "Level 4 qualifications and above",
        "2" = "Level 3 qualifications",
        "3" = "Level 2 qualifications",
        "4" = "Level 1 qualifications", # "CSEs or equivalent"
        "5" = "Level 1 qualifications",  # "NVQ or HND or equivalent"
        "6" = "Level 4 qualifications and above", # "Other prof qualifications eg: nursing, teaching"
       "-7" = "No qualifications"
      ) 
    ) 
  
  ukbb_props <- famhist_nowt %>% 
    filter(! is.na(qual)) %>% 
    group_by(age, qual) %>% 
    summarize(n = n()) %>% 
    mutate(age_n = sum(n), percent = n/age_n) %>% 
    select(-n, -age_n)
  
  
  prop_table <- ukbb_props %>% left_join(census_props, by = c("age", "qual"), 
        suffix = c(".ukbb", ".census"))
  prop_table$educ_weight <- prop_table$percent.census/prop_table$percent.ukbb
  
  famhist <- famhist_nowt %>% left_join(
          prop_table %>% select(age, qual, educ_weight), 
          by = c("age", "qual")
        ) %>% 
        select(-age, -qual)
  
  return(famhist)
}
