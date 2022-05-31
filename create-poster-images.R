
library(dplyr)
library(ggplot2)

# ==== setup ====

library(drake)
library(dplyr)
library(magrittr)
library(forcats)
library(ggplot2)
library(tidyr)
loadNamespace("scales")
loadNamespace("shades")


drake::loadd(rgs)
drake::loadd(score_names)

theme_set(theme_minimal())
theme_update(text = element_text(family = "Abadi MT Condensed Light"))

my_hline <- geom_hline(yintercept = 0, colour = "grey20", linetype = "dotted")
my_vline <- geom_vline(xintercept = 0, colour = "grey20", linetype = "dotted")


standard_ggplot <- function (
    dfr, 
    fill_col, 
    n_regs, 
    ...,
    score_col = quo(term), 
    order_idx = 1, 
    fill_direction = 1,
    n_cats    = NULL,
    conf_int  = NULL
) {
  if (! missing(score_col)) score_col <- enquo(score_col)
  
  mfc <- missing(fill_col)
  fill_col <- if (mfc) quo(NULL) else enquo(fill_col)
  if (missing(n_cats)) {
    n_cats <- if (mfc) 1 else length(unique(pull(dfr, {{fill_col}})))
  }
  
  if (is.null(conf_int)) conf_int <- "conf.low" %in% names(dfr)
  
  if (conf_int) {
    confint_geom_segment <- geom_segment(
      aes(
        x        = conf.low, 
        xend     = conf.high,
        linetype = "95% c.i. uncorrected"
      ), 
      color   = "grey45", 
      alpha   = 0.2, 
      size    = 1.1, 
      lineend = "round"
    )
  } else {
    confint_geom_segment <- NULL
  }
  
  my_shape_fill_guide <- guide_legend(order = 1, title = NULL)
  
  n_regs <- as.double(n_regs)
  shape_values <- paste(c("circle", "square", "triangle", "diamond", 
                          "triangle down"), "filled")
  shape_values <- shape_values[1:n_cats]
  if (fill_direction == -1) shape_values <- rev(shape_values)
  dfr %>% mutate(
    !! score_col := pretty_names(!! score_col),
    !! score_col := fct_reorder(!! score_col, estimate, order_abs(order_idx))
  ) %>% 
    rename(p = p.value) %>% 
    ggplot(aes(estimate, {{score_col}}, yend = {{score_col}}, 
               shape = {{fill_col}}, ...)) +
    my_vline + 
    confint_geom_segment +
    geom_point(size = 1.8, alpha = 0.8, aes(color = {{fill_col}}, 
                                            fill = stage(p < 0.05/{{n_regs}}, after_scale = ifelse(
                                              fill == "white", fill, color)))) +
    my_fill_scale(aesthetics = "color", n = n_cats, direction = fill_direction, 
                  guide = my_shape_fill_guide) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"), 
                      guide = guide_legend(
                        title = sprintf("p < 0.05/%s", n_regs),
                        override.aes = list(shape = "circle filled"),
                        order = 2
                      )) + 
    scale_shape_manual(values = shape_values, guide = my_shape_fill_guide) +
    labs(x = "β RLRS", y = "") +
    theme(
      legend.justification = "top", 
      panel.grid.major.y   = element_blank(),
      plot.title.position = "plot",
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) + guides(
      linetype = "none"
    )
}


my_fill_scale <- function (n, aesthetics = "fill", direction = 1, ...) {
  
  fill_colors <- c("steelblue4", "darkred")
  
  scale <- if (n <= 2) {
    dots <- list(...)
    dots$direction <- NULL
    do.call(scale_colour_manual, c(list(values = fill_colors[1:n],
                                        aesthetics = aesthetics), 
                                   dots))
  } else {
    # scale <- scale_color_viridis_d(aesthetics = aesthetics, option = "D",
    #                                  direction = direction, end = 0.6, ...)
    fill_colors <- c(fill_colors[1], "#008b8b", fill_colors[2])
    if (direction != -1) fill_colors <- rev(fill_colors)
    scale_color_manual(aesthetics = aesthetics, 
                       values = shades::gradient(fill_colors,
                                                 steps = n, space = "Lab"), 
                       ...) 
  }
  
  scale
}


order_abs <- function (n = 1) {
  function (x) x[n] 
}


pretty_names <- function (names) {
  pretty <- c(
    ADHD_2017               = "ADHD",
    age_at_menarche         = "Age at menarche",
    age_at_menopauze        = "Age at menopause",
    agreeableness           = "Agreeableness",
    ai_substance_use        = "Age at smoking initiation",
    alcohol_schumann        = "Alcohol use",
    alzheimer               = "Alzheimer",
    autism_2017             = "Autism",
    bipolar                 = "Bipolar",
    bmi_combined            = "BMI",
    body_fat                = "Body Fat",
    caffeine                = "Caffeine",
    cannabis                = "Cannabis (ever vs. never)",
    cognitive_ability       = "Cognitive Ability",
    conscientiousness       = "Conscientiousness",
    coronary_artery_disease = "Coronary Artery Disease",
    cpd_substance_use       = "Cigarettes per day",
    diagram_T2D             = "Type 2 Diabetes",
    dpw_substance_use       = "Drinks per week",
    EA2_noUKB               = "Educ. attainment 2 (no UKBB)",
    EA3_excl_23andMe_UK     = "Educ. attainment 3 (no UK)",
    eating_disorder         = "Eating disorder",
    extraversion            = "Extraversion",
    height_combined         = "Height",
    hip_combined            = "Hip circumference",
    MDD_PGC2_noUKB          = "Major Depressive Disorder",
    neuroticism             = "Neuroticism",
    openness                = "Openness",
    sc_substance_use        = "Smoking cessation",
    SCZ2                    = "Schizophrenia",
    si_substance_use        = "Smoking initiation",
    wc_combined             = "Waist circumference",
    whr_combined            = "Waist-hip ratio"
  )
  
  pretty[names]
}



# ==== Main plot ====


drake::loadd(res_wt_van_alten_weights)
drake::loadd(res_unweighted)

res_wt_combined <- bind_rows(
  None                = res_unweighted %>% 
    select(-score_name) %>% 
    filter(term != "(Intercept)"),
  Main                = res_wt_van_alten_weights,
  .id = "Weights"
)

n_regs <- as.double(length(score_names))

res_wt_combined %>% 
  filter(Weights %in% c("Main", "None")) %>% 
  mutate(
    Weights  = recode(Weights, "Main" = "Weighted", "None" = "Unweighted"),
  ) %>% 
  standard_ggplot(fill_col = Weights, n_regs = n_regs, fill_direction = -1,
                  order_idx = 2)

ggsave("poster-images/main.jpeg", width = 13, height = 8, bg = "white",
         units = "cm")

# ==== Correlations ====
# 
drake::loadd(res_all)
drake::loadd(res_cor_educ)
drake::loadd(res_cor_income)

effect_size <-  res_all %>% 
  filter(dep.var == "RLRS", reg.type == "controlled") %>% 
  pull(estimate)


dfr <- data.frame(
  PGS          = rep(score_names, 2),
  Correlation  = rep(c("Earnings", "Education"), each = 33),
  cor          = c(res_cor_income[,1], res_cor_educ[,1]),
  effect_size  = rep(effect_size, 2)
)

dfr %>% 
  mutate(
    PGS = ifelse(abs(effect_size) > 0.027, PGS, ""),
    PGS = sub("(^\\w*?)_.*", "\\1", PGS)
  ) %>% 
  ggplot(aes(effect_size, cor)) + 
  my_hline +
  my_vline +
  geom_point(aes(color = Correlation), size = 2) +
  geom_text(aes(label = PGS), size = 3, colour = "black", 
            check_overlap = TRUE, nudge_x = 0.008) +
  facet_wrap( ~ Correlation, 
              scales = "free") + 
  labs(x = "β RLRS", y = "Corr.") +
  coord_cartesian(clip = "off") +
  theme(
    text = element_text(size = 18),
    legend.position = "none",
    panel.spacing = unit(20, "pt"),
    panel.grid = element_blank()
  ) + 
  scale_color_manual(values = c("steelblue", "darkred")) +
  scale_x_continuous(breaks = c(-0.04, 0, 0.04), 
                     limits = c(-.04, 0.04))

ggsave("poster-images/correlations.jpeg", width = 13, height = 8,
       units = "cm", bg = "white")


# ==== Income ====
# 

drake::loadd(res_income)

n_regs <- as.double(nrow(res_income))
res_income %>% 
  mutate(
    Income = factor(income_cat, 
                    labels = c("< £18K", "£18-30K", "£31-51K", "£52-100K", "> £100K")
    )
  ) %>% 
  standard_ggplot(n_regs = n_regs, fill_col = `Income`)

ggsave("poster-images/income.jpeg", width = 13, height = 8, bg = "white",
       units = "cm")

# ==== Education ====
# 
drake::loadd(res_edu)
n_regs <- as.double(nrow(res_edu))

res_edu %>% 
  mutate(
    "Age left FTE" = fct_relevel(age_fte_cat, "< 16", "16-18", "> 18"),
  ) %>% 
  standard_ggplot(n_regs = n_regs, fill_col = `Age left FTE`)

ggsave("poster-images/education.jpeg", width = 13, height = 8, bg = "white",
       units = "cm")

# ==== Single parenthood ====
# 

drake::loadd(res_with_partner)
res_with_partner %<>% filter(grepl(":", term)) 

n_regs <- as.double(nrow(res_with_partner))

res_with_partner %>% 
  mutate(
    Household = ifelse(grepl("TRUE", term), "With partner", 
                       "Without partner")
  ) %>% 
  standard_ggplot(fill_col = Household, n_regs = n_regs, 
                  score_col = score_name, fill_direction = -1)

ggsave("poster-images/partner.jpeg", width = 13, height = 8, bg = "white",
       units = "cm")


# ==== AFLB ====

drake::loadd(res_age_flb_cross)

# this counts each cross term as a separate test:
n_regs <- as.double(nrow(res_age_flb_cross)) 

res_age_flb_cross %>% 
  mutate(
    `Age at first live birth` = gsub("age_flb_cat(.*):.*", "\\1", term),
  ) %>% 
  standard_ggplot(score_col = score_name, 
                  fill_col = `Age at first live birth`, n_regs = n_regs)

ggsave("poster-images/AFLB.jpeg", width = 13, height = 8, bg = "white",
       units = "cm")

# ==== Mediation ====
# 

drake::loadd(res_mediation)

res_mediation %<>% mutate(
  # 0.025 not 0.05 because two-sided
  sig     = abs(statistic_ind) > qnorm(1 - 0.025/
                                         nrow(res_mediation)),
  pos_sig = sig & prop_ind > 0
)

res_mediation %>% 
  mutate(
    term = pretty_names(term),
    term = fct_reorder(term, prop_ind, order_abs(1)),
    ci_contains_0 = prop_ind_conf_low < 0 & prop_ind_conf_high> 0,
    prop_ind_conf_low = ifelse(ci_contains_0, NA, prop_ind_conf_low),
    prop_ind_conf_high = ifelse(ci_contains_0, NA, prop_ind_conf_high)
  ) %>% 
  ggplot(aes(y = term)) + 
  geom_col(aes(x = prop_ind), width = .5, fill = "steelblue4") + 
  geom_errorbar(aes(
      xmin = prop_ind_conf_low, 
      xmax = prop_ind_conf_high,
      linetype = "Bootstrap 95% c.i. uncorrected",
    ),
    width = 0.2,
    color = "black"
  ) +
  my_vline +
  scale_x_continuous(labels = scales::percent, n.breaks = 7) +
  coord_cartesian(xlim = c(-0.5, 1)) +
  labs(x = "Proportion mediated", y = "") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14)
  ) +
  guides(
    linetype = "none"
  )

ggsave("poster-images/mediation.jpeg", width = 10, height = 8, bg = "white",
       units = "cm")

# ==== Income-EA3 barplot ====

drake::loadd(famhist)

income_EA3 <- famhist %>% 
  filter(! is.na(n_children), ! is.na(income_cat)) %>% 
  mutate(
    `Household income` = factor(income_cat, 
                                labels = c("< £18K", "£18-30K", "£31-51K", "£52-100K", 
                                           "> £100K"))
  ) %>% 
  group_by(`Household income`) %>% 
  summarize(
    `Without selection` = mean(EA3_excl_23andMe_UK, na.rm = TRUE),
    Actual              = weighted.mean(EA3_excl_23andMe_UK, n_children, 
                                        na.rm = TRUE)
  ) %>% 
  tidyr::pivot_longer(-`Household income`) %>% 
  rename(
    `Children's mean EA3` = value
  )

income_EA3 %>% 
  ggplot(aes(`Household income`, `Children's mean EA3`, fill = name)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("steelblue4", "grey65")) +
  theme(
    legend.background = element_rect(
      fill = "white", 
      colour = "black",
      size = rel(0.5)
    ),
    legend.margin = margin(t = 0, b = 5.5, l = 5.5, 
                           r = 5.5),
    legend.title = element_blank(), 
    legend.position = c(.25, .85),
    panel.grid = element_blank()
  )

ggsave("poster-images/income-EA3.jpeg", width = 13, height = 8, bg = "white",
       units = "cm")
