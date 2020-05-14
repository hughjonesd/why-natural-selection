
library(dplyr)
library(glue)
library(readr)
loadd(famhist)
loadd(score_names)

grf <- readr::read_table2("~/UKBB data 2019/relatedness_file.txt")
grf <- grf %>% filter(Kinship > 0.45)
sibs <- left_join(grf, famhist, by = c("ID1" = "f.eid"))
sibs <- left_join(sibs, famhist, by = c("ID2" = "f.eid"))


diff_scores <- c(
        score_names,
        paste0(score_names, "_resid"),
        "age_fulltime_edu", 
        "YearsEdu.ISCED"
      )

for (ds in diff_scores){
  dsx <- paste0(ds, ".x")
  dsy <- paste0(ds, ".y")
  sibs[paste0("diff_", ds)] <- sibs[dsx] - sibs[dsy]
}

target_phenotypes <- list(
  EA3_excl_23andMe_UK = list("age_fulltime_edu", "YearsEdu.ISCED")
)

models <- list()
for (sn in score_names) {
  phenos <- target_phenotypes[[sn]]
  if (is.null(phenos)) next
  
  models[[sn]] <- list()
  if (! is.list(phenos)) phenos <- list(phenos)
  for (pheno in phenos) {
    fml_raw <- as.formula(glue::glue("diff_{pheno} ~ diff_{sn}"))
    fml_resid <- as.formula(glue::glue("diff_{pheno} ~ diff_{sn}_resid"))
    models[[sn]][[pheno]] <- list(
      raw   = lm(fml_raw, sibs),
      resid = lm(fml_resid, sibs)
    )
  }
}

summary(lm(diff_YearsEdu.ISCED ~ diff_EA3_excl_23andMe_UK, sibs))
summary(lm(diff_YearsEdu.ISCED ~ diff_EA3_excl_23andMe_UK_resid, sibs))
summary(lm(diff_age_fulltime_edu ~ diff_EA3_excl_23andMe_UK, sibs))
summary(lm(diff_age_fulltime_edu ~ diff_EA3_excl_23andMe_UK_resid, sibs))
