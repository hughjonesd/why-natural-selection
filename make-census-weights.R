
suppressPackageStartupMessages({
  loadNamespace("readr")
  loadNamespace("tidyr")
  library(dplyr)
  library(sf)
})

make_census_age_qual <- function (DC1108EW_file, DC5202SC_file) {
  ew_age_qual <- readr::read_csv(DC1108EW_file, col_names = FALSE)
  # cols 1-3 are date, geography and geography code. We only want 
  # row 1: "England and Wales"
  ew_age_qual <- t(ew_age_qual)
  ew_age_qual %<>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    select(name = 1, value = 2) %>% 
    slice(-(1:3)) %>% 
    filter(
      ! grepl("All categories", name, fixed = TRUE)
    ) %>% 
    tidyr::separate(name, 
                    into = c("age_cat", "qual", "ethnic", "measure"), 
                    sep  = ";"
    ) %>% 
    select(-measure) %>% 
    mutate(
      age_cat = sub("Age: Age ", "", age_cat),
      qual    = sub(".*: ", "", qual),
      qual    = sub("(Level.*) qualifications", "\\1", qual),
      ethnic  = sub("Ethnic Group: *", "", ethnic)
    ) %>% 
    filter(ethnic == " White: Total") %>% 
    select(-ethnic) %>% 
    mutate(
      value = as.numeric(value)
    )
  
  # alternatively, "Ethnic Group: White: English" would
  # capture English/Welsh/Scottish/Northern Irish/British only
  
  # England splits up age categories to 65-74, 75+. Scotland doesn't.
  ew_age_qual %<>% 
    tidyr::pivot_wider(
      id_cols     = qual,
      names_from  = age_cat, 
      values_from = value
    ) %>% 
    mutate(
      `65 and over` = `65 to 74` + `75 and over`
    ) %>% 
    select(- `65 to 74`, - `75 and over`) %>% 
    tidyr::pivot_longer(
      -qual,
      names_to  = "age_cat",
      values_to = "value"
    ) %>% 
    mutate(
      place = "England and Wales"
    )
  
  sc_age_qual <- readr::read_csv(DC5202SC_file, skip = 12)
  sc_age_qual %<>% 
    select(
      age_cat = `X1`,              # new names reflect values in the columns
      qual    = `Ethnicity (Flat)`, 
      value   = `White: Total`
    ) %>% 
    tidyr::fill(age_cat) %>% 
    slice(-1) %>%                  # removes a second row of headings
    filter(
      ! grepl("All people", age_cat),
      ! grepl("Total", qual),
      ! is.na(qual)                # removes some notes at the bottom
    ) %>% 
    mutate(
      age_cat = sub("Aged ", "", age_cat, fixed = TRUE),
      age_cat = sub(":", "", age_cat, fixed = TRUE),
      place = "Scotland"
    )
  
  census_age_qual <- bind_rows(sc_age_qual, ew_age_qual)
  # since 4.0.0 could use `proportions()`
  
  # merge the countries to form Great Britain
  # we exclude categories that don't appear in the famhist data
  # i.e. young people and those with "Other" or "Apprenticeship" qualifications
  census_age_qual %<>% 
        filter(age_cat %in% c("35 to 49", "50 to 64", "65 and over")) %>% 
        filter(! qual %in% c("Other qualifications", "Apprenticeship")) %>% 
        group_by(age_cat, qual) %>% 
        summarize(pop = sum(value)) %>% 
        ungroup() %>% 
        mutate(percent = pop/sum(pop))
  
  return(census_age_qual)
}


find_containing_msoas <- function (famhist, msoa_shapefile) {
  msoa_shape <- st_read(msoa_shapefile)
  
  famhist %<>% 
    filter(! is.na(home_lon.0.0), ! is.na(home_lat.0.0)) %>% 
    st_as_sf(
      coords = c("home_lon.0.0", "home_lat.0.0"), 
      crs = "EPSG:27700" # OSGB: 1936, see 
      # https://epsg.io/27700
    )
  
  msoa <- st_intersects(famhist, msoa_shape)
  # a few famhist rows intersect 2 msoas. These were always neighbours 
  # (checked manually); we take the first and don't worry about it.
  # a 1000 or so don't intersect any. Here, rounding seems to have
  # put people out at sea :-) again, we don't worry about it
  # (using `st_is_within_distance()` would take forever)
  msoa <- map_dbl(msoa, ~ if (length(.x) == 0) NA_real_ else .x[[1]] )
  msoa <- msoa_shape$geo_code[msoa]
  
  data.frame(f.eid = famhist$f.eid, msoa = msoa)
}


make_census_msoa <- function (census_msoa_file) {
  census_msoa <-  readr::read_csv(census_msoa_file) 
  census_msoa %<>% 
        dplyr::select(
          -matches("Total"), 
          -matches("All categories"), 
          -matches("All persons")
        ) %>% 
        tidyr::pivot_longer(
          -c(date, geography, `geography code`),
          names_to = c("age", "sex", "hh_type", NA),
          names_sep = "; "
        ) %>% 
        rename(
          msoa_code = `geography code`
        ) %>% 
        mutate(
          age     = gsub("Age: Age ", "", age, fixed = TRUE),
          sex     = case_when(
                      (sex == "Sex: Males") ~ 1,
                      (sex == "Sex: Females") ~ 0
                    ),
          hh_type = gsub("living arrangements: ", "", hh_type, fixed = TRUE),
          with_partner = grepl("Living in a couple", hh_type)
        ) %>% 
        filter(
          as.numeric(gsub("^(\\d+).*", "\\1", age)) >= 40
        )
  
  census_msoa %<>% 
        group_by(msoa_code, age, sex, with_partner) %>% 
        summarize(pop = sum(value), .groups = "drop") %>% 
        ungroup() %>% 
        mutate(percent = pop/sum(pop))
  
  return(census_msoa)
}
