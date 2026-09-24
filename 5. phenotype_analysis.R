library(dplyr)
library(purrr)
library(tibble)
library(survival)
library(survey)

saved <- readRDS("data/HF_continuum_imputed_cohorts.rds")
imputed_cohorts <- saved$imputed_cohorts
cohort <- saved$cohort
# cox_formula <- saved$cox_formula
# cox_formula_m3 <- saved$cox_formula_m3
n_imputations <- saved$n_imputations
message("Loaded imputed_cohorts (n = ", n_imputations, ")")

## ============================================================
## Mortality trajectories across phenotypes
## ============================================================

## ---- 9.1a Unadjusted survey-weighted Kaplan-Meier curves ------------
km_design <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~wtmec6yr,
                       nest = TRUE, data = cohort)
km_fit <- svykm(Surv(permth_exm, mortstat) ~ hf_continuum_label, design = km_design) #se=TRUE

message("\n---- 9.1 Survey-weighted Kaplan-Meier: unadjusted survival at 60/120 months ----")
km_summary <- purrr::imap_dfr(km_fit, function(fit, label) {
  tibble(
    hf_continuum_label = label,
    surv_60mo = suppressWarnings(approx(fit$time, fit$surv, xout = 60, method = "constant", rule = 2)$y),
    surv_120mo = suppressWarnings(approx(fit$time, fit$surv, xout = 120, method = "constant", rule = 2)$y)
  )
})
print(km_summary)


## ---- 9.1a Unadjusted Kaplan-Meier curves (light) -----------
if (exists("imputed_cohorts")) rm(imputed_cohorts)
if (exists("boot_results_list")) rm(boot_results_list)
gc()

km_fit_light <- survfit(Surv(permth_exm, mortstat) ~ hf_continuum_label,
                        data = cohort, weights = wtmec6yr)

km_summary_wide <- summary(km_fit_light, times = c(60, 120))
km_summary <- tibble(
  hf_continuum_label = gsub("hf_continuum_label=", "", km_summary_wide$strata),
  time = km_summary_wide$time,
  surv = km_summary_wide$surv
) %>%
  tidyr::pivot_wider(names_from = time, values_from = surv, names_prefix = "surv_") %>%
  rename(surv_60mo = surv_60, surv_120mo = surv_120)

message("\n---- 9.1 Case-weighted Kaplan-Meier: unadjusted survival at 60/120 months ----")
print(km_summary)

km_curve_data <- tibble(
  hf_continuum_label = rep(names(km_fit_light$strata), km_fit_light$strata),
  time = km_fit_light$time,
  surv = km_fit_light$surv,
  lower = km_fit_light$lower,
  upper = km_fit_light$upper
) %>%
  mutate(hf_continuum_label = gsub("hf_continuum_label=", "", hf_continuum_label))


km_plot <- ggplot2::ggplot(km_curve_data, ggplot2::aes(x = time, y = surv, color = hf_continuum_label,
                                                       fill = hf_continuum_label)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  ggplot2::geom_step(linewidth = 0.8) +
  ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  ggplot2::scale_x_continuous(limits = c(0, 240), breaks = seq(0, 240, by = 60), expand = c(0.01, 0.01)) +
  ggplot2::labs(
    x = "Follow-up time (months)", y = "Survival probability",
    color = "HF-continuum phenotype", fill = "HF-continuum phenotype",
    title = "Weighted Kaplan-Meier survival by HF phenotype"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom",
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave("data/HF_continuum_KM_curves.png", km_plot, width = 12, height = 6, dpi = 300)
message("Saved: data/HF_continuum_KM_curves.png (case-weighted, with approximate 95% CI)")

mortality_rate_table <- cohort %>%
  group_by(hf_continuum_label) %>%
  summarise(
    n = n(),
    deaths_unweighted = sum(mortstat == 1, na.rm = TRUE),
    person_years_unweighted = round(sum(permth_exm, na.rm = TRUE) / 12, 1),
    deaths_weighted = round(sum(wtmec6yr * (mortstat == 1), na.rm = TRUE), 0),
    person_years_weighted = round(sum(wtmec6yr * permth_exm, na.rm = TRUE) / 12, 0),
    .groups = "drop"
  ) %>%
  mutate(
    unweighted_rate_per_1000py = round(1000 * deaths_unweighted / person_years_unweighted, 2),
    weighted_rate_per_1000py = round(1000 * deaths_weighted / person_years_weighted, 2)
  )

message("\n---- 9.1b Deaths / person-years / mortality rate by phenotype ----")
print(mortality_rate_table)

## ---- 9.2 Adjusted survey-weighted Cox models, pooled across MI -----
##   Model 1: age, sex, race/ethnicity
##   Model 2: Model 1 + education, poverty-income ratio, insurance
##   Model 3: Model 2 + hypertension, diabetes, obesity, CKD, smoking, CHD/MI/stroke

cox_formula_m0 <- Surv(permth_exm, mortstat) ~ hf_continuum_label
cox_formula_m1 <- Surv(permth_exm, mortstat) ~ hf_continuum_label +
  ridageyr + riagendr + ridreth1
cox_formula_m2 <- update(cox_formula_m1, . ~ . + dmdeduc2 + indfmpir + hiq011)
cox_formula_m3 <- update(cox_formula_m2, . ~ . + flag_hypertension + flag_diabetes +
                           flag_obesity + flag_ckd + flag_smoking + flag_chd_mi + flag_stroke)

cox_formula <- cox_formula_m3

fit_pooled_cox <- function(formula, label) {
  message("\n---- Fitting ", label, " across ", n_imputations, " imputed datasets ----")
  models <- purrr::map(imputed_cohorts, function(d) {
    d <- d %>% mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                               flag_hypertension, flag_diabetes, flag_obesity,
                               flag_ckd, flag_smoking, flag_chd_mi, flag_stroke,
                               hiq011), as.factor))
    design_i <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~wtmec6yr,
                          nest = TRUE, data = d)
    svycoxph(formula, design = design_i)
  })
  pooled <- mitools::MIcombine(models)
  
  summ <- summary(pooled) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    mutate(model = label, HR = exp(results), HR_lower = exp(`(lower`), HR_upper = exp(`upper)`)) %>%
    select(model, term, HR, HR_lower, HR_upper, everything())
  list(models = models, pooled = pooled, summary = summ)
}

cox_m0 <- fit_pooled_cox(cox_formula_m0, "Model 0 (Unadjusted)")
cox_m1 <- fit_pooled_cox(cox_formula_m1, "Model 1 (age/sex/race)")
cox_m2 <- fit_pooled_cox(cox_formula_m2, "Model 2 (+ education/PIR/insurance)")
cox_m3 <- fit_pooled_cox(cox_formula_m3, "Model 3 (+ HTN/DM/obesity/CKD/smoking/CHD-MI/stroke)")

cox_models <- cox_m3$models
cox_pooled <- cox_m3$pooled
cox_pooled_summary <- cox_m3$summary

message("\n---- Adjusted Cox models (Models 1-3), pooled across ", n_imputations, " imputations ----")
cox_all_models_summary <- bind_rows(cox_m0$summary, cox_m1$summary, cox_m2$summary, cox_m3$summary)
print(cox_all_models_summary %>% filter(grepl("hf_continuum_label", term)))

## ---- Forest Plots: ggforest-style, ALL covariates, one figure per model ----
term_labels <- c(
  "hf_continuum_labelStage A (at risk)" = "Stage A (at risk)",
  "hf_continuum_labelStage B (pre-HF, biomarker)" = "Stage B (pre-HF, biomarker)",
  "hf_continuum_labelStage C (self-reported HF)" = "Stage C (self-reported HF)",
  "ridageyr" = "Age, per year",
  "riagendr2" = "Female (vs Male)",
  "ridreth12" = "Other Hispanic (vs Mexican American)",
  "ridreth13" = "Non-Hispanic White (vs Mexican American)",
  "ridreth14" = "Non-Hispanic Black (vs Mexican American)",
  "ridreth15" = "Other Race (vs Mexican American)",
  "dmdeduc22" = "9-11th Grade (vs <9th Grade)",
  "dmdeduc23" = "High School Grad/GED (vs <9th Grade)",
  "dmdeduc24" = "Some College/AA (vs <9th Grade)",
  "dmdeduc25" = "College Graduate+ (vs <9th Grade)",
  "indfmpir" = "Poverty-Income Ratio, per unit",
  "hiq0112" = "Uninsured (vs Insured)",
  "flag_hypertension1" = "Hypertension (Yes vs No)",
  "flag_diabetes1" = "Diabetes (Yes vs No)",
  "flag_obesity1" = "Obesity (Yes vs No)",
  "flag_ckd1" = "CKD, self-report (Yes vs No)",
  "flag_smoking1" = "Current Smoking (Yes vs No)",
  "flag_chd_mi1" = "CHD/MI History (Yes vs No)",
  "flag_stroke1" = "Stroke History (Yes vs No)"
)

term_order <- rev(names(term_labels))

forest_data <- cox_all_models_summary %>%
  filter(term %in% names(term_labels)) %>%
  mutate(
    term_label = factor(term_labels[term], levels = unname(term_labels[term_order])),
    label_text = sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper)
  )

hr_range <- range(c(forest_data$HR_lower, forest_data$HR_upper), na.rm = TRUE)
hr_limits_shared <- c(hr_range[1] * 0.9, hr_range[2] * 1.1)

build_classic_forest_plot <- function(model_data, model_title, x_limits = NULL) {
  model_data <- model_data %>% mutate(term_label = droplevels(term_label))
  n_rows <- nlevels(model_data$term_label)
  
  panel_labels <- ggplot2::ggplot(model_data, ggplot2::aes(x = 0, y = term_label)) +
    ggplot2::geom_text(ggplot2::aes(label = term_label), hjust = 0, size = 3.4) +
    ggplot2::xlim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(r = 2, l = 4))
  
  ## log
  x_limits_log <- c(0.5, 10)

  panel_forest <- ggplot2::ggplot(model_data, ggplot2::aes(x = HR, y = term_label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = pmax(HR_lower, x_limits_log[1]),
                                          xmax = pmin(HR_upper, x_limits_log[2])),
                             color = "#2C3E50", size = 0.55, linewidth = 0.7) +
    ggplot2::geom_segment(
      data = model_data %>% dplyr::filter(HR_upper > x_limits_log[2]),
      ggplot2::aes(x = x_limits_log[2] * 0.85, xend = x_limits_log[2] * 0.98,
                   y = term_label, yend = term_label),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.08, "inches")),
      color = "#2C3E50", linewidth = 0.7
    ) +
    ggplot2::scale_x_log10(
      limits = x_limits_log,
      breaks = c(0.5, 1, 2, 4, 8, 16),
      labels = scales::label_number(accuracy = 0.1, drop0trailing = TRUE)
    ) +
    ggplot2::labs(x = "Hazard Ratio (log scale)", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  
  panel_hr_text <- ggplot2::ggplot(model_data, ggplot2::aes(x = 0, y = term_label)) +
    ggplot2::geom_text(ggplot2::aes(label = label_text), hjust = 0.5, size = 3.4) +
    ggplot2::xlim(-0.5, 0.5) +
    ggplot2::labs(title = "HR (95% CI)") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(l = 2, r = 2),
                   plot.title = ggplot2::element_text(size = 9, hjust = 0.5))
  
  combined <- tryCatch({
    p <- patchwork::wrap_plots(panel_labels, panel_forest, panel_hr_text,
                               widths = c(2.2, 2.3, 1.2))
    p + patchwork::plot_annotation(
      title = model_title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5))
    )
  }, error = function(e) {
    message("  patchwork combination failed (", conditionMessage(e), ") -- falling back to",
            " gridExtra::grid.arrange().")
    if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
    gridExtra::arrangeGrob(
      panel_labels, panel_forest, panel_hr_text,
      ncol = 3, widths = c(2.2, 2.3, 1.2),
      top = grid::textGrob(model_title, gp = grid::gpar(fontface = "bold", fontsize = 12))
    )
  })
  
  list(plot = combined, n_rows = n_rows)
}

model_captions <- c(
  "Model 0 (Unadjusted)" =
    "Fig. Forest plot for the unadjusted Cox proportional hazards model. HR: hazard ratio; CI: confidence interval. No covariates adjusted.",
  "Model 1 (age/sex/race)" =
    "Fig. Forest plot for Model 1. HR: hazard ratio; CI: confidence interval. Model 1 included age, sex, and race/ethnicity.",
  "Model 2 (+ education/PIR/insurance)" =
    "Fig. Forest plot for Model 2. HR: hazard ratio; CI: confidence interval. Model 2 was adjusted as for Model 1 and additionally included education, poverty-income ratio, and health insurance status.",
  "Model 3 (+ HTN/DM/obesity/CKD/smoking/CHD-MI/stroke)" =
    "Fig. Forest plot for Model 3. HR: hazard ratio; CI: confidence interval. Model 3 was adjusted as for Model 2 and additionally included hypertension, diabetes, obesity, chronic kidney disease, current smoking, and history of coronary heart disease/myocardial infarction or stroke."
)

model_file_suffix <- c(
  "Model 0 (Unadjusted)" = "Model0_Unadjusted",
  "Model 1 (age/sex/race)" = "Model1",
  "Model 2 (+ education/PIR/insurance)" = "Model2",
  "Model 3 (+ HTN/DM/obesity/CKD/smoking/CHD-MI/stroke)" = "Model3"
)

built_plots <- list()
for (m in unique(forest_data$model)) {
  built <- build_classic_forest_plot(forest_data %>% filter(model == m), m, x_limits = hr_limits_shared)
  built_plots[[m]] <- built
  fname <- paste0("data/HF_continuum_ForestPlot_", model_file_suffix[m], ".png")
  fig_height <- max(3, 0.45 * built$n_rows + 1.2)
  ggplot2::ggsave(fname, built$plot, width = 9.5, height = fig_height, dpi = 300, limitsize = FALSE)
  message("Saved: ", fname, " (", built$n_rows, " terms)")
  message("  Suggested caption: ", model_captions[m])
}


## ============================================================
## Model 3 Pairwise comparison：Stage B vs Stage A, Stage C vs Stage B, Stage C vs Stage A
## ============================================================
fit_model3_releveled <- function(ref_level, imputed_cohorts, cox_formula_m3) {
  models <- purrr::map(imputed_cohorts, function(d) {
    d <- d %>%
      mutate(
        hf_continuum_label = relevel(as.factor(hf_continuum_label), ref = ref_level),
        across(c(riagendr, ridreth1, dmdeduc2, flag_hypertension, flag_diabetes,
                 flag_obesity, flag_ckd, flag_smoking, flag_chd_mi, flag_stroke,
                 hiq011), as.factor)
      )
    design_i <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~wtmec6yr,
                          nest = TRUE, data = d)
    svycoxph(cox_formula_m3, design = design_i)
  })
  pooled <- mitools::MIcombine(models)
  summary(pooled) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    mutate(
      reference = ref_level,
      HR = exp(results), HR_lower = exp(`(lower`), HR_upper = exp(`upper)`),
      z = results / se,
      p_value = 2 * pnorm(-abs(z))
    ) %>%
    select(reference, term, HR, HR_lower, HR_upper, p_value, everything())
}

## ---- 1. Stage A as ref: Stage B vs A & Stage C vs A -----
message("\n---- Stage A as ref ----")
fit_ref_A <- fit_model3_releveled("Stage A (at risk)", imputed_cohorts, cox_formula_m3)
print(fit_ref_A %>% filter(grepl("^hf_continuum_label", term)) %>%
        select(term, HR, HR_lower, HR_upper, p_value))

## ---- 2. Stage B as ref：Stage C vs B ----
message("\n---- Stage B as ref ----")
fit_ref_B <- fit_model3_releveled("Stage B (pre-HF, biomarker)", imputed_cohorts, cox_formula_m3)
print(fit_ref_B %>% filter(grepl("^hf_continuum_label", term)) %>%
        select(term, HR, HR_lower, HR_upper, p_value))

## ----------- 3. Combine -------------------------
pairwise_comparisons <- bind_rows(
  ## vs No apparent risk（from cox_pooled_summary, existing）
  cox_pooled_summary %>%
    filter(grepl("^hf_continuum_label", term)) %>%
    transmute(
      comparison = paste0(sub("^hf_continuum_label", "", term), " vs No apparent HF risk"),
      HR, HR_lower, HR_upper
    ),
  ## New：Stage A
  fit_ref_A %>%
    filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
    transmute(comparison = "Stage B vs Stage A", HR, HR_lower, HR_upper),
  fit_ref_A %>%
    filter(term == "hf_continuum_labelStage C (self-reported HF)") %>%
    transmute(comparison = "Stage C vs Stage A", HR, HR_lower, HR_upper),
  ## New：Stage B
  fit_ref_B %>%
    filter(term == "hf_continuum_labelStage C (self-reported HF)") %>%
    transmute(comparison = "Stage C vs Stage B", HR, HR_lower, HR_upper)
) %>%
  mutate(HR_CI = sprintf("%.2f (%.2f\u2013%.2f)", HR, HR_lower, HR_upper))

p_values_lookup <- bind_rows(
  fit_ref_A %>% filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
    transmute(comparison = "Stage B vs Stage A", p_value),
  fit_ref_A %>% filter(term == "hf_continuum_labelStage C (self-reported HF)") %>%
    transmute(comparison = "Stage C vs Stage A", p_value),
  fit_ref_B %>% filter(term == "hf_continuum_labelStage C (self-reported HF)") %>%
    transmute(comparison = "Stage C vs Stage B", p_value)
)

pairwise_comparisons <- pairwise_comparisons %>%
  left_join(p_values_lookup, by = "comparison")

print(pairwise_comparisons %>% select(comparison, HR_CI, p_value))

writexl::write_xlsx(pairwise_comparisons, "data/HF_continuum_pairwise_stage_comparisons.xlsx")
message("\nSaved: data/HF_continuum_pairwise_stage_comparisons.xlsx")


## ---- 9.3 Adjusted survival curves standardized to the cohort -------
rep_data <- imputed_cohorts[[1]] %>%
  mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                  flag_hypertension, flag_diabetes, flag_obesity,
                  flag_ckd, flag_smoking, flag_chd_mi, flag_stroke,
                  hiq011), as.factor))

time_grid_93 <- seq(0, 120, by = 2)

## Case-weighted Cox fit 
fit_93 <- do.call("coxph", list(
  formula = cox_formula, data = rep_data, weights = rep_data$wtmec6yr,
  x = TRUE, na.action = na.exclude
))

X_full_93 <- model.matrix(fit_93)
beta_hat_93 <- coef(fit_93)
X_full_93 <- X_full_93[, names(beta_hat_93), drop = FALSE]

used_rows_93 <- as.integer(rownames(X_full_93))
weights_vec_93 <- rep_data$wtmec6yr[used_rows_93]
stopifnot(nrow(X_full_93) == length(weights_vec_93))

bh_93 <- basehaz(fit_93, centered = FALSE)
H0_grid_93 <- suppressWarnings(approx(bh_93$time, bh_93$hazard, xout = time_grid_93,
                                      method = "constant", rule = 2)$y)

pheno_cols_93 <- grep("^hf_continuum_label", names(beta_hat_93), value = TRUE)
levels_93 <- levels(rep_data$hf_continuum_label)

standardized_curves <- purrr::map_dfr(levels_93, function(lvl) {
  X_lvl <- X_full_93
  X_lvl[, pheno_cols_93] <- 0
  col_name <- paste0("hf_continuum_label", lvl)
  if (col_name %in% pheno_cols_93) X_lvl[, col_name] <- 1
  
  lp <- as.numeric(X_lvl %*% beta_hat_93)
  exp_lp <- exp(lp)
  
  ## Correctly weighted at every time point -- this is the fix.
  mean_surv <- sapply(H0_grid_93, function(h0) weighted.mean(exp(-h0 * exp_lp), w = weights_vec_93))
  tibble(hf_continuum_label = lvl, time = time_grid_93, mean_surv = mean_surv)
})

message("\n---- 9.3 Standardized adjusted survival at 60/120 months (weighted marginal standardization, CORRECTED) ----")
standardized_summary <- standardized_curves %>%
  group_by(hf_continuum_label) %>%
  summarise(
    surv_60mo = suppressWarnings(approx(time, mean_surv, xout = 60, method = "constant", rule = 2)$y),
    surv_120mo = suppressWarnings(approx(time, mean_surv, xout = 120, method = "constant", rule = 2)$y),
    .groups = "drop"
  )
print(standardized_summary)

## ---- 9.3b Bootstrap CI for the corrected weighted point estimate -----
cluster_bootstrap_ids_93 <- function(data, psu_var = "sdmvpsu", strata_var = "sdmvstra") {
  strata_ids <- unique(data[[strata_var]])
  purrr::map_dfr(strata_ids, function(s) {
    psus <- unique(data[[psu_var]][data[[strata_var]] == s])
    sampled <- sample(psus, length(psus), replace = TRUE)
    purrr::map_dfr(sampled, function(p) {
      data %>% dplyr::filter(.data[[strata_var]] == s, .data[[psu_var]] == p) %>% dplyr::select(seqn)
    })
  })
}

run_one_marginal_boot_93 <- function(boot_data) {
  fit_b <- tryCatch(
    do.call("coxph", list(formula = cox_formula, data = boot_data,
                          weights = boot_data$wtmec6yr, x = TRUE, na.action = na.exclude)),
    error = function(e) {
      message("    coxph() failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit_b)) return(NULL)
  
  X_full_b <- model.matrix(fit_b)
  beta_hat_b <- coef(fit_b)
  X_full_b <- X_full_b[, names(beta_hat_b), drop = FALSE]
  
  used_rows_b <- as.integer(rownames(X_full_b))
  weights_vec_b <- boot_data$wtmec6yr[used_rows_b]
  if (nrow(X_full_b) != length(weights_vec_b)) return(NULL)
  
  bh_b <- basehaz(fit_b, centered = FALSE)
  H0_grid_b <- suppressWarnings(approx(bh_b$time, bh_b$hazard, xout = time_grid_93,
                                       method = "constant", rule = 2)$y)
  
  pheno_cols_b <- grep("^hf_continuum_label", names(beta_hat_b), value = TRUE)
  lvls_b <- levels(boot_data$hf_continuum_label)
  result_list <- vector("list", length(lvls_b))
  
  for (i in seq_along(lvls_b)) {
    lvl <- lvls_b[i]
    X_lvl <- X_full_b
    X_lvl[, pheno_cols_b] <- 0
    col_name <- paste0("hf_continuum_label", lvl)
    if (col_name %in% pheno_cols_b) X_lvl[, col_name] <- 1
    
    lp <- as.numeric(X_lvl %*% beta_hat_b)
    exp_lp <- exp(lp)
    surv_vec <- sapply(H0_grid_b, function(h0) weighted.mean(exp(-h0 * exp_lp), w = weights_vec_b))
    result_list[[i]] <- tibble(hf_continuum_label = lvl, time = time_grid_93, surv = surv_vec)
  }
  bind_rows(result_list)
}

n_boot_93 <- 200

message("\n---- 9.3b Bootstrap CI (B = ", n_boot_93, ", weighted matrix computation) ----")
t0_93 <- Sys.time()
boot_curves_list_93 <- vector("list", n_boot_93)
for (b in seq_len(n_boot_93)) {
  if (b %% 20 == 0) message("  Bootstrap ", b, " of ", n_boot_93)
  boot_ids_93 <- cluster_bootstrap_ids_93(rep_data %>% dplyr::select(seqn, sdmvpsu, sdmvstra))
  boot_data_b_93 <- boot_ids_93 %>% dplyr::left_join(rep_data, by = "seqn", relationship = "many-to-many")
  boot_curves_list_93[[b]] <- tryCatch(run_one_marginal_boot_93(boot_data_b_93), error = function(e) {
    message("    Replicate ", b, " failed: ", conditionMessage(e))
    NULL
  })
}
message("Bootstrap runtime: ", round(as.numeric(difftime(Sys.time(), t0_93, units = "mins")), 1), " min")

boot_curves_93 <- bind_rows(boot_curves_list_93, .id = "replicate")
message("Successful replicates: ", length(unique(boot_curves_93$replicate)), " / ", n_boot_93)
if (nrow(boot_curves_93) == 0) stop("boot_curves_93 is empty -- check messages above.")

standardized_curves_ci_93 <- boot_curves_93 %>%
  group_by(hf_continuum_label, time) %>%
  summarise(
    surv_lower = quantile(surv, 0.025, na.rm = TRUE),
    surv_upper = quantile(surv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## ---- Figure 4: adjusted cumulative mortality through 10 years -------
## Point estimate and CI now both computed via the identical weighted matrix algorithm
fig4_data <- standardized_curves %>%
  filter(time <= 120) %>%
  left_join(standardized_curves_ci_93, by = c("hf_continuum_label", "time")) %>%
  mutate(
    cum_mortality = 1 - mean_surv,
    cum_mortality_lower = 1 - surv_upper,
    cum_mortality_upper = 1 - surv_lower,
    hf_continuum_label = factor(hf_continuum_label,
                                levels = c("No apparent HF risk", "Stage A (at risk)",
                                           "Stage B (pre-HF, biomarker)", "Stage C (self-reported HF)"))
  )

message("\nCI non-NA row check:")
print(fig4_data %>% group_by(hf_continuum_label) %>%
        summarise(non_na = sum(!is.na(cum_mortality_lower)), total = n(), .groups = "drop"))

fig4_plot <- ggplot2::ggplot(fig4_data, ggplot2::aes(x = time, y = cum_mortality,
                                                     color = hf_continuum_label, fill = hf_continuum_label)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = cum_mortality_lower, ymax = cum_mortality_upper),
                       alpha = 0.15, color = NA, na.rm = TRUE) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 24)) +
  ggplot2::labs(
    x = "Follow-up time (months)", y = "Adjusted cumulative mortality",
    color = "HF-continuum phenotype", fill = "HF-continuum phenotype",
    title = "Adjusted cumulative mortality through 10 years by HF-continuum phenotype",
    subtitle = paste0("Weighted marginal standardization; 95% CI via cluster bootstrap (B=", n_boot_93, ")")
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave("data/HF_continuum_Figure4_adjusted_cumulative_mortality.png", fig4_plot,
                width = 12, height = 6, dpi = 300)
message("Saved: data/HF_continuum_Figure4_adjusted_cumulative_mortality.png (weighted marginal standardization, with bootstrap CI)")


## ============================================================
## 9.3c Estimand documentation and covariate-overlap diagnostic
## ============================================================
message("\n---- 9.3c Covariate overlap check (age, PIR) across phenotype groups ----")
overlap_summary <- rep_data %>%
  group_by(hf_continuum_label) %>%
  summarise(
    age_min = min(ridageyr, na.rm = TRUE), age_median = median(ridageyr, na.rm = TRUE),
    age_max = max(ridageyr, na.rm = TRUE),
    pir_min = min(indfmpir, na.rm = TRUE), pir_median = median(indfmpir, na.rm = TRUE),
    pir_max = max(indfmpir, na.rm = TRUE),
    .groups = "drop"
  )
print(overlap_summary)

ref_age_range <- overlap_summary %>% filter(hf_continuum_label == "No apparent HF risk") %>%
  summarise(lo = age_min, hi = age_max)
overlap_pct_in_ref_age_range <- rep_data %>%
  group_by(hf_continuum_label) %>%
  summarise(
    pct_within_ref_age_range = round(100 * mean(ridageyr >= ref_age_range$lo & ridageyr <= ref_age_range$hi,
                                                na.rm = TRUE), 1),
    .groups = "drop"
  )
message("\nPercent of each group falling within the no-apparent-risk group's observed age range:")
print(overlap_pct_in_ref_age_range)
message("A low percentage for Stage B/C indicates the standardized estimate for those groups relies on",
        " extrapolation beyond the age range actually observed in the no-apparent-risk reference group --",
        " report this explicitly as a limitation on the standardized estimates' reliability.")


## ============================================================
## Conventional (3-level, no-biomarker) phenotype construction
## ============================================================
imputed_cohorts_conventional <- purrr::map(imputed_cohorts, function(d) {
  d <- d %>%
    mutate(stage_a_model3_flags = as.numeric(
      flag_hypertension == 1 | flag_diabetes == 1 | flag_obesity == 1 |
        flag_ckd == 1 | flag_smoking == 1 | flag_chd_mi == 1 | flag_stroke == 1
    ))
  d %>% mutate(
    conventional_label = factor(
      case_when(
        mcq160b == 1 ~ "Stage C (self-reported HF)",
        stage_a_model3_flags == 1 ~ "Stage A (at risk)",
        TRUE ~ "No apparent HF risk"
      ),
      levels = c("No apparent HF risk", "Stage A (at risk)", "Stage C (self-reported HF)")
    )
  )
})

message("\n---- Verification: conventional_label missingness (should be 0) ----")
print(sapply(imputed_cohorts_conventional[1:min(3, length(imputed_cohorts_conventional))],
             function(d) sum(is.na(d$conventional_label))))

## ---- reassignment breakdown -----
reassignment_check <- imputed_cohorts_conventional[[1]] %>%
  filter(hf_continuum_label == "Stage B (pre-HF, biomarker)") %>%
  count(conventional_label, name = "n") %>%
  mutate(pct = round(100 * n / sum(n), 1))

message("\n---- Reassignment of original Stage B participants under the conventional",
        " (no-NT-proBNP) label ----")
message("Total original Stage B participants: ",
        sum(imputed_cohorts_conventional[[1]]$hf_continuum_label == "Stage B (pre-HF, biomarker)"))
print(reassignment_check)

## ---- Same-N verification -----------------------
same_n_check <- purrr::map2_dfr(imputed_cohorts, imputed_cohorts_conventional,
                                function(d_full, d_conv) {
                                  tibble(n_full = nrow(d_full), n_conv = nrow(d_conv), same_n = nrow(d_full) == nrow(d_conv))
                                })
message("\n---- Same-participants verification across all ", n_imputations, " imputations ----")
print(same_n_check)
if (!all(same_n_check$same_n)) {
  stop("Full and conventional models are NOT evaluated on the same participants",
       " in at least one imputation -- investigate before proceeding.")
}

conventional_formula <- update(cox_formula_m3, . ~ . - hf_continuum_label + conventional_label)
message("\n---- Model formulas (identical Model 3 covariate set, phenotype term differs) ----")
message("Full model (with Stage B): ", deparse(cox_formula_m3))
message("Conventional model (no NT-proBNP): ", deparse(conventional_formula))


## ============================================================
## 9.4 Restricted mean survival time (RMST) -- unadjusted and adjusted
## ============================================================
rmst <- function(time, surv, tau) {
  ord <- order(time)
  time <- time[ord]; surv <- surv[ord]
  keep <- time <= tau
  knot_times <- c(0, time[keep], tau)
  last_surv_before_tau <- if (any(keep)) tail(surv[keep], 1) else 1
  knot_surv <- c(1, surv[keep], last_surv_before_tau)
  sum(diff(knot_times) * head(knot_surv, -1))
}

## ---- Unadjusted RMST (from survey-weighted KM curves, 10.1) ----------
rmst_results <- purrr::imap_dfr(km_fit, function(fit, label) {
  tibble(
    hf_continuum_label = label,
    rmst_60mo = rmst(fit$time, fit$surv, 60),
    rmst_120mo = rmst(fit$time, fit$surv, 120)
  )
})

ref_rmst_120 <- rmst_results$rmst_120mo[rmst_results$hf_continuum_label == "No apparent HF risk"]
rmst_diff_results <- rmst_results %>%
  mutate(rmst_120mo_diff_vs_ref = round(rmst_120mo - ref_rmst_120, 2))

message("\n---- 9.4a Unadjusted RMST (months) and 10y difference vs reference ----")
print(rmst_diff_results %>% select(hf_continuum_label, rmst_60mo, rmst_120mo, rmst_120mo_diff_vs_ref))

## ---- Adjusted RMST (from §9.3's corrected, weighted standardized_curves) ----
rmst_results_adjusted <- standardized_curves %>%
  group_by(hf_continuum_label) %>%
  summarise(
    rmst_60mo_adj = rmst(time, mean_surv, 60),
    rmst_120mo_adj = rmst(time, mean_surv, 120),
    .groups = "drop"
  )

ref_rmst_120_adj <- rmst_results_adjusted$rmst_120mo_adj[
  rmst_results_adjusted$hf_continuum_label == "No apparent HF risk"
]
rmst_diff_results_adjusted <- rmst_results_adjusted %>%
  mutate(rmst_120mo_diff_vs_ref_adj = round(rmst_120mo_adj - ref_rmst_120_adj, 2))

message("\n---- 9.4b Adjusted RMST (months, from corrected weighted marginal standardization)",
        " and 10y difference vs reference ----")
print(rmst_diff_results_adjusted)


## ============================================================
## Table 3: Mortality risk by HF-continuum phenotype
## ============================================================
## Columns: Phenotype, Deaths/person-years, Weighted mortality rate,
## Adjusted 10-year mortality risk (from corrected §9.3), Adjusted RMST
## difference.
table3_final <- mortality_rate_table %>%
  left_join(
    standardized_summary %>%
      mutate(adjusted_10y_mortality_risk_pct = round(100 * (1 - surv_120mo), 1)) %>%
      select(hf_continuum_label, adjusted_10y_mortality_risk_pct),
    by = "hf_continuum_label"
  ) %>%
  left_join(
    rmst_diff_results_adjusted %>% select(hf_continuum_label, rmst_120mo_diff_vs_ref_adj),
    by = "hf_continuum_label"
  ) %>%
  rename(rmst_diff_10y_months_adjusted = rmst_120mo_diff_vs_ref_adj)

message("\n---- Table 3 (final): Mortality risk by HF-continuum phenotype ----")
message("(Unadjusted mortality rate; Adjusted 10-year risk and RMST difference both from the",
        " corrected, weighted marginal-standardization pipeline in 9.3/9.4b; HR columns are in",
        " the Forest Plot figure, not repeated here.)")
print(table3_final)


## ============================================================
## 9.5 C-statistic: with Stage B vs conventional model
## ============================================================
c_stat_list <- vector("list", length(imputed_cohorts))
for (i in seq_along(imputed_cohorts)) {
  message("  C-statistic: imputation ", i, " of ", length(imputed_cohorts))
  
  d_full <- imputed_cohorts[[i]] %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
  d_conv <- imputed_cohorts_conventional[[i]] %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, conventional_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
  
  fit_full <- coxph(cox_formula_m3, data = d_full, weights = d_full$wtmec6yr)
  fit_conv <- coxph(conventional_formula, data = d_conv, weights = d_conv$wtmec6yr)
  
  c_stat_list[[i]] <- tibble(
    imputation = i,
    c_with_stageB = concordance(fit_full)$concordance,
    c_conventional = concordance(fit_conv)$concordance
  )
}
c_stat_results <- bind_rows(c_stat_list)

c_stat_summary <- c_stat_results %>%
  summarise(
    mean_c_with_stageB = mean(c_with_stageB),
    mean_c_conventional = mean(c_conventional),
    mean_c_improvement = mean(c_with_stageB - c_conventional)
  )

message("\n---- 9.5 C-statistic: with Stage B vs conventional (Stage A/C only) model ----")
print(c_stat_summary)


## ============================================================
## 9.6 Calibration slope, IDI, continuous NRI (10-year horizon)
## ============================================================
rep_data_full <- imputed_cohorts[[1]] %>%
  mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                  flag_hypertension, flag_diabetes, flag_obesity,
                  flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
rep_data_conv <- imputed_cohorts_conventional[[1]] %>%
  mutate(across(c(riagendr, ridreth1, dmdeduc2, conventional_label,
                  flag_hypertension, flag_diabetes, flag_obesity,
                  flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))

compute_point_estimates_one_imputation <- function(imp_idx) {
  rep_data_full_i <- imputed_cohorts[[imp_idx]] %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
  rep_data_conv_i <- imputed_cohorts_conventional[[imp_idx]] %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, conventional_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
  
  fit_full_i <- do.call("coxph", list(formula = cox_formula_m3, data = rep_data_full_i,
                                      weights = rep_data_full_i$wtmec6yr, na.action = na.exclude))
  fit_conv_i <- do.call("coxph", list(formula = conventional_formula, data = rep_data_conv_i,
                                      weights = rep_data_conv_i$wtmec6yr, na.action = na.exclude))
  
  lp_full_i <- predict(fit_full_i, type = "lp")
  calib_full_i <- unname(coef(coxph(Surv(rep_data_full_i$permth_exm, rep_data_full_i$mortstat) ~ lp_full_i,
                                    na.action = na.exclude))[1])
  lp_conv_i <- predict(fit_conv_i, type = "lp")
  calib_conv_i <- unname(coef(coxph(Surv(rep_data_conv_i$permth_exm, rep_data_conv_i$mortstat) ~ lp_conv_i,
                                    na.action = na.exclude))[1])
  
  bh_full_i <- basehaz(fit_full_i, centered = FALSE)
  H0_full_i <- suppressWarnings(approx(bh_full_i$time, bh_full_i$hazard, xout = 120,
                                       method = "constant", rule = 2)$y)
  pred_full_i <- 1 - exp(-H0_full_i * exp(lp_full_i))
  
  bh_conv_i <- basehaz(fit_conv_i, centered = FALSE)
  H0_conv_i <- suppressWarnings(approx(bh_conv_i$time, bh_conv_i$hazard, xout = 120,
                                       method = "constant", rule = 2)$y)
  pred_conv_i <- 1 - exp(-H0_conv_i * exp(lp_conv_i))
  
  dat_i <- tibble(
    event_10y = case_when(
      rep_data_full_i$mortstat == 1 & rep_data_full_i$permth_exm <= 120 ~ 1,
      rep_data_full_i$permth_exm >= 120 ~ 0,
      TRUE ~ NA_real_
    ),
    pred_full = pred_full_i, pred_conv = pred_conv_i
  ) %>% filter(!is.na(event_10y) & !is.na(pred_full) & !is.na(pred_conv))
  
  ev_i <- dat_i %>% filter(event_10y == 1)
  ne_i <- dat_i %>% filter(event_10y == 0)
  
  idi_i <- (mean(ev_i$pred_full) - mean(ev_i$pred_conv)) - (mean(ne_i$pred_full) - mean(ne_i$pred_conv))
  nri_i <- (mean(ev_i$pred_full > ev_i$pred_conv) - mean(ev_i$pred_full < ev_i$pred_conv)) +
    (mean(ne_i$pred_full < ne_i$pred_conv) - mean(ne_i$pred_full > ne_i$pred_conv))
  
  tibble(imputation = imp_idx, calibration_slope_full = calib_full_i,
         calibration_slope_conv = calib_conv_i, idi = idi_i, nri = nri_i, analytic_n = nrow(dat_i))
}

message("\n---- 9.6 Calibration slope, IDI, continuous NRI across all ", n_imputations, " imputed datasets ----")
mi_point_estimates <- purrr::map_dfr(seq_len(n_imputations), compute_point_estimates_one_imputation)
print(mi_point_estimates)

calibration_slope_full <- mean(mi_point_estimates$calibration_slope_full, na.rm = TRUE)
calibration_slope_conv <- mean(mi_point_estimates$calibration_slope_conv, na.rm = TRUE)
idi <- mean(mi_point_estimates$idi, na.rm = TRUE)
nri <- mean(mi_point_estimates$nri, na.rm = TRUE)

discrimination_summary <- tibble(
  metric = c("Calibration slope (with Stage B)", "Calibration slope (conventional)",
             "IDI (with Stage B vs conventional)", "Continuous NRI (with Stage B vs conventional)"),
  value = round(c(calibration_slope_full, calibration_slope_conv, idi, nri), 4)
)
print(discrimination_summary)


## ============================================================
## 9.7 Bootstrap 95% CIs (cluster bootstrap, PSU within stratum)
## ============================================================
n_boot <- 200

cluster_bootstrap_ids <- function(data, psu_var = "sdmvpsu", strata_var = "sdmvstra") {
  strata_ids <- unique(data[[strata_var]])
  purrr::map_dfr(strata_ids, function(s) {
    psus_in_stratum <- unique(data[[psu_var]][data[[strata_var]] == s])
    sampled_psus <- sample(psus_in_stratum, length(psus_in_stratum), replace = TRUE)
    purrr::map_dfr(sampled_psus, function(p) {
      data %>% filter(.data[[strata_var]] == s, .data[[psu_var]] == p) %>% select(seqn)
    })
  })
}

run_one_bootstrap <- function(boot_full, boot_conv) {
  fit_full_b <- tryCatch(coxph(cox_formula_m3, data = boot_full, weights = boot_full$wtmec6yr,
                               na.action = na.exclude), error = function(e) NULL)
  fit_conv_b <- tryCatch(coxph(conventional_formula, data = boot_conv, weights = boot_conv$wtmec6yr,
                               na.action = na.exclude), error = function(e) NULL)
  if (is.null(fit_full_b) || is.null(fit_conv_b)) return(NULL)
  
  c_full_b <- tryCatch(concordance(fit_full_b)$concordance, error = function(e) NA_real_)
  c_conv_b <- tryCatch(concordance(fit_conv_b)$concordance, error = function(e) NA_real_)
  
  lp_full_b <- predict(fit_full_b, type = "lp")
  calib_full_b <- tryCatch(unname(coef(coxph(Surv(boot_full$permth_exm, boot_full$mortstat) ~ lp_full_b,
                                             na.action = na.exclude))[1]), error = function(e) NA_real_)
  lp_conv_b <- predict(fit_conv_b, type = "lp")
  calib_conv_b <- tryCatch(unname(coef(coxph(Surv(boot_conv$permth_exm, boot_conv$mortstat) ~ lp_conv_b,
                                             na.action = na.exclude))[1]), error = function(e) NA_real_)
  
  bh_full_b <- basehaz(fit_full_b, centered = FALSE)
  H0_full_b <- suppressWarnings(approx(bh_full_b$time, bh_full_b$hazard, xout = 120,
                                       method = "constant", rule = 2)$y)
  pred_full_b <- 1 - exp(-H0_full_b * exp(lp_full_b))
  
  bh_conv_b <- basehaz(fit_conv_b, centered = FALSE)
  H0_conv_b <- suppressWarnings(approx(bh_conv_b$time, bh_conv_b$hazard, xout = 120,
                                       method = "constant", rule = 2)$y)
  pred_conv_b <- 1 - exp(-H0_conv_b * exp(lp_conv_b))
  
  idi_data_b <- tibble(
    event_10y = case_when(
      boot_full$mortstat == 1 & boot_full$permth_exm <= 120 ~ 1,
      boot_full$permth_exm >= 120 ~ 0,
      TRUE ~ NA_real_
    ),
    pred_full = pred_full_b, pred_conv = pred_conv_b
  ) %>% filter(!is.na(event_10y) & !is.na(pred_full) & !is.na(pred_conv))
  
  if (nrow(idi_data_b) < 10 || sum(idi_data_b$event_10y == 1) < 2) {
    idi_b <- NA_real_; nri_b <- NA_real_
  } else {
    ev_b <- idi_data_b %>% filter(event_10y == 1)
    ne_b <- idi_data_b %>% filter(event_10y == 0)
    idi_b <- (mean(ev_b$pred_full) - mean(ev_b$pred_conv)) - (mean(ne_b$pred_full) - mean(ne_b$pred_conv))
    nri_b <- (mean(ev_b$pred_full > ev_b$pred_conv) - mean(ev_b$pred_full < ev_b$pred_conv)) +
      (mean(ne_b$pred_full < ne_b$pred_conv) - mean(ne_b$pred_full > ne_b$pred_conv))
  }
  
  tibble(c_with_stageB = c_full_b, c_conventional = c_conv_b,
         calibration_slope_full = calib_full_b, calibration_slope_conv = calib_conv_b,
         idi = idi_b, nri = nri_b)
}

message("\n---- 9.7 Bootstrap 95% CIs (B = ", n_boot, ") ----")
boot_results_list <- vector("list", n_boot)
for (b in seq_len(n_boot)) {
  if (b %% 20 == 0) message("  Bootstrap replicate ", b, " of ", n_boot)
  boot_ids <- cluster_bootstrap_ids(rep_data_full %>% select(seqn, sdmvpsu, sdmvstra))
  boot_full <- boot_ids %>% left_join(rep_data_full, by = "seqn", relationship = "many-to-many")
  boot_conv <- boot_ids %>% left_join(rep_data_conv, by = "seqn", relationship = "many-to-many")
  
  boot_results_list[[b]] <- tryCatch(run_one_bootstrap(boot_full, boot_conv), error = function(e) {
    message("    Bootstrap replicate ", b, " failed (", conditionMessage(e), ") -- skipped.")
    NULL
  })
}
boot_results <- bind_rows(boot_results_list)
message("Successful bootstrap replicates: ", nrow(boot_results), " of ", n_boot)

boot_ci <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(0.025, 0.975)))
}

rubins_pool <- function(within_var, between_estimates) {
  point <- mean(between_estimates, na.rm = TRUE)
  between_var <- var(between_estimates, na.rm = TRUE)
  m <- sum(!is.na(between_estimates))
  total_var <- within_var + (1 + 1 / m) * between_var
  se <- sqrt(total_var)
  tibble(point_estimate = point, ci_lower = point - 1.96 * se, ci_upper = point + 1.96 * se)
}

pooled_c_full <- rubins_pool(var(boot_results$c_with_stageB, na.rm = TRUE), c_stat_results$c_with_stageB)
pooled_c_conv <- rubins_pool(var(boot_results$c_conventional, na.rm = TRUE), c_stat_results$c_conventional)
pooled_calib_full <- rubins_pool(var(boot_results$calibration_slope_full, na.rm = TRUE),
                                 mi_point_estimates$calibration_slope_full)
pooled_calib_conv <- rubins_pool(var(boot_results$calibration_slope_conv, na.rm = TRUE),
                                 mi_point_estimates$calibration_slope_conv)
pooled_idi <- rubins_pool(var(boot_results$idi, na.rm = TRUE), mi_point_estimates$idi)
pooled_nri <- rubins_pool(var(boot_results$nri, na.rm = TRUE), mi_point_estimates$nri)

discrimination_summary_ci <- tibble(
  metric = c("C-statistic (with Stage B)", "C-statistic (conventional)",
             "Calibration slope (with Stage B)", "Calibration slope (conventional)",
             "IDI (with Stage B vs conventional)", "Continuous NRI (with Stage B vs conventional)"),
  point_estimate = round(c(pooled_c_full$point_estimate, pooled_c_conv$point_estimate,
                           pooled_calib_full$point_estimate, pooled_calib_conv$point_estimate,
                           pooled_idi$point_estimate, pooled_nri$point_estimate), 4),
  boot_ci_lower = round(c(pooled_c_full$ci_lower, pooled_c_conv$ci_lower,
                          pooled_calib_full$ci_lower, pooled_calib_conv$ci_lower,
                          pooled_idi$ci_lower, pooled_nri$ci_lower), 4),
  boot_ci_upper = round(c(pooled_c_full$ci_upper, pooled_c_conv$ci_upper,
                          pooled_calib_full$ci_upper, pooled_calib_conv$ci_upper,
                          pooled_idi$ci_upper, pooled_nri$ci_upper), 4)
)
message("\n---- Discrimination/calibration summary (Rubin's-rules-pooled) ----")
print(discrimination_summary_ci)


## ============================================================
## Save Results outputs
## ============================================================
writexl::write_xlsx(
  list(
    "Table3 Final" = table3_final,
    "Mortality Rate Table" = mortality_rate_table,
    "KM Unadjusted Summary" = km_summary,
    "Cox Models 0-3 (MI)" = cox_all_models_summary,
    "Standardized Adj Survival" = standardized_summary,
    "RMST Unadjusted" = rmst_diff_results,
    "RMST Adjusted" = rmst_diff_results_adjusted,
    "Stage B Reassignment" = reassignment_check,
    "C-statistic Stage B" = c_stat_summary,
    "Calibration IDI NRI" = discrimination_summary_ci,
    "Fig3 Stratified Prevalence" = fig3_data
  ),
  "data/HF_continuum_trajectory_and_sensitivity.xlsx"
)
message("\nSaved: data/HF_continuum_trajectory_and_sensitivity.xlsx")


## ============================================================
## Table 4: Incremental prognostic value of Stage B/pre-HF
## ============================================================
fmt_ci <- function(point, lower, upper, digits = 3) {
  if (is.na(lower) || is.na(upper)) return(sprintf(paste0("%.", digits, "f"), point))
  sprintf(paste0("%.", digits, "f", " (%.", digits, "f", "-%.", digits, "f", ")"), point, lower, upper)
}
get_ci <- function(metric_name) {
  row <- discrimination_summary_ci %>% filter(metric == metric_name)
  c(row$point_estimate, row$boot_ci_lower, row$boot_ci_upper)
}

c_conv  <- get_ci("C-statistic (conventional)")
c_full  <- get_ci("C-statistic (with Stage B)")
calib_conv <- get_ci("Calibration slope (conventional)")
calib_full <- get_ci("Calibration slope (with Stage B)")
idi_vals <- get_ci("IDI (with Stage B vs conventional)")
nri_vals <- get_ci("Continuous NRI (with Stage B vs conventional)")
c_improvement_ci <- boot_ci(boot_results$c_with_stageB - boot_results$c_conventional)

table4 <- tibble(
  `Reference Model` = c(
    "Model 3 covariates + 3-level phenotype (Stage B participants reassigned to Stage A or no apparent risk using remaining criteria)",
    NA
  ),
  `Model + Stage B` = c(NA, "Same covariates + 4-level phenotype restoring biomarker-defined Stage B"),
  Metric = c("C-statistic (10-year, 95% CI)", "C-statistic"),
  Reference = c(fmt_ci(c_conv[1], c_conv[2], c_conv[3]), NA),
  `With Stage B` = c(NA, fmt_ci(c_full[1], c_full[2], c_full[3]))
)

table4_wide <- tibble(
  Metric = c("N (participants)", "C-statistic (95% CI)", "Delta C-statistic (95% CI)",
             "Calibration slope (95% CI)", "IDI (95% CI)", "Continuous NRI (95% CI)"),
  `Reference model (no NT-proBNP)` = c(
    as.character(nrow(rep_data_conv)),
    fmt_ci(c_conv[1], c_conv[2], c_conv[3]),
    "Reference",
    fmt_ci(calib_conv[1], calib_conv[2], calib_conv[3]),
    "\u2014", "\u2014"
  ),
  `Model + Stage B (with NT-proBNP)` = c(
    as.character(nrow(rep_data_full)),
    fmt_ci(c_full[1], c_full[2], c_full[3]),
    fmt_ci(c_stat_summary$mean_c_improvement, c_improvement_ci[1], c_improvement_ci[2]),
    fmt_ci(calib_full[1], calib_full[2], calib_full[3]),
    fmt_ci(idi_vals[1], idi_vals[2], idi_vals[3]),
    fmt_ci(nri_vals[1], nri_vals[2], nri_vals[3])
  )
)

message("\n---- Table 4: Incremental prognostic value of Stage B/pre-HF ----")
message("Both models fit on N=", nrow(rep_data_full), " identical participants (verified above).",
        " Stage B participants in the reference model are reassigned to Stage A or no apparent",
        " risk (see 'Stage B Reassignment' sheet), not excluded.",
        " Prediction horizon: 10 years (120 months). Pooled across ", n_imputations,
        " imputed datasets; 95% CIs from cluster bootstrap (B=", n_boot,
        ") combined with between-imputation variance via Rubin's rules.")
print(table4_wide)

writexl::write_xlsx(
  list("Table4" = table4_wide, "Stage B Reassignment" = reassignment_check,
       "Same-N Check" = same_n_check),
  "data/HF_continuum_Table4_incremental_prognostic_value.xlsx"
)
message("Saved: data/HF_continuum_Table4_incremental_prognostic_value.xlsx")