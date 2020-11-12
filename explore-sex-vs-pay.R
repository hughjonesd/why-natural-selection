
map_dfr(score_names, 
        ~ (paste("n_children ~ sex * fj_q + sex * fj_q:", .x) %>% 
             as.formula %>% 
             lm(data = famhist) %>% 
             tidy()
           ), .id = "score_name") -> tmp
tmp %>% 
  mutate(female = grepl("sex", term), pay = sub(".*?(\\d+%).*", "\\1", term)) %>% 
  mutate(term_bit = gsub(".*:", "", term)) %>% 
  filter(term_bit == score_name) %>% 
  select(-term_bit) %>% 
  ggplot(aes(pay, estimate, color = female, group = female)) + 
    geom_line(position = "dodge") + 
    facet_wrap(vars(score_name), ) +
    geom_hline(yintercept = 0)

map_dfr(score_names, 
        ~ (paste("n_children ~ sex * ", .x) %>% 
             as.formula %>% 
             lm(data = famhist) %>% 
             tidy()), .id = "score_name") -> tmp2
map_dfr(score_names, 
        ~ (paste("first_job_pay ~ factor(sex):", .x) %>% 
             as.formula %>% lm(data = famhist) %>% 
             tidy()), .id = "score_name") -> tmp3

tmp2 <- left_join(tmp2, tmp3, by = c("score_name", "term"), suffix = c(".kids", ".pay"))
tmp2 <- tmp2 %>% filter(term != "(Intercept)") %>% mutate(female = grepl("sex.0", term))
ggplot(tmp2, aes(estimate.kids, estimate.pay, color = female)) + 
  geom_point() +
  geom_text(aes(label = score_name))


# explore controlling for pay and/or education
fml <- paste("n_children ~ ", paste(score_names, collapse = " + "))
fml <- as.formula(fml)
fml <- update(fml, .~.+sex+age_at_recruitment+I(age_at_recruitment^2))
fml <- update(fml, . ~ . - EA2_noUKB - wc_combined - whr_combined - hip_combined - autism_2017 - body_fat)

tmp <- lm(fml, data=famhist)
tmp2 <- lm(update(fml, .~.+log(first_job_pay) + factor(edu_qual)), data= famhist)

tmp <- tidy(tmp)
tmp2 <- tidy(tmp2)
tmp3 <- inner_join(tmp, tmp2, by = "term")

tmp3 %<>% mutate(prop_change = estimate.y/estimate.x)
tmp3 %>% 
  filter(p.value.x < 0.05, term %in% score_names) %>% 
  mutate(
    `P value` = santoku::chop(p.value.y, c(0.05/33, 0.01, 0.05), c("< 0.05/33", "< 0.01", "< 0.05", "n.s."))
  ) %>% 
  arrange(desc(prop_change)) %>% 
  ggplot(aes(abs(estimate.x), abs(estimate.y))) + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
    geom_point(aes(color = `P value`)) + 
    scale_color_viridis_d(direction = -1) +
    geom_text(aes(label = term), size = 2, nudge_y = -.0005) +
    xlim(0, .03) + 
    ylim(0, .03)

# couples, effect of male and female PGS when combined
# economic theory might suggest more negative selection for females
mf_pairs %<>% filter(
  n_children.f == n_children.m,
  with_partner.f, 
  with_partner.m
  )
couples <- map_dfr(
    score_names, 
    ~ {
      fm <- as.formula(paste0("n_children.f ~ ", .x, ".f + ", .x, ".m + age_at_recruitment.f"))
      tidy(lm(fm, mf_pairs))
    }, 
    .id = "score_name"
  ) 
couples %<>% 
  filter(term!="(Intercept)", term != "age_at_recruitment.f") %>% 
  mutate(
    female = grepl(".f$", term),
    score_name = forcats::fct_reorder(score_name, estimate, function (x) x[2])
  ) %>% 
  arrange(score_name, female)
ggplot(couples, aes(estimate, score_name, color = female)) + 
  geom_point() + 
  geom_vline(xintercept = 0)

# positive at 0.28 but insignificant
cor.test(couples$estimate[couples$female], couples$estimate[! couples$female])
