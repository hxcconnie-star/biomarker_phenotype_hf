## ============================================================
## Sensitivity Analyses
## ============================================================
##   1. Alternative NT-proBNP threshold
##   2. Excluding advanced CKD
##   3. Excluding early deaths
##   4. Obesity-stratified analysis
##   5. Complete-case analysis
##   6. Age-as-time-scale Cox model

library(dplyr)
library(purrr)
library(tibble)
library(survival)
library(survey)
library(mitools)
library(writexl)

stopifnot(exists("cohort"), exists("imputed_cohorts"), exists("cox_formula_m3"),
          exists("cox_pooled_summary"))

fit_pooled_cox_sens <- function(formula, data_list, label) {
  models <- purrr::map(data_list, function(d) {
    design_i <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~wtmec6yr,
                          nest = TRUE, data = d)
    svycoxph(formula, design = design_i)
  })
  pooled <- mitools::MIcombine(models)
  summary(pooled) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    mutate(sensitivity = label, HR = exp(results),
           HR_lower = exp(`(lower`), HR_upper = exp(`upper)`)) %>%
    select(sensitivity, term, HR, HR_lower, HR_upper, everything())
}

check_cell_counts <- function(data, phenotype_var = "hf_continuum_label",
                              label = "", min_events = 10) {
  tab <- data %>%
    group_by(.data[[phenotype_var]]) %>%
    summarise(n = n(), events = sum(mortstat == 1, na.rm = TRUE), .groups = "drop")
  message("  Cell counts for '", label, "':")
  print(tab)
  thin_cells <- tab %>% filter(events < min_events)
  if (nrow(thin_cells) > 0) {
    message("  *** WARNING: ", nrow(thin_cells), " phenotype level(s) have <", min_events,
            " events -- HRs referencing or estimating these levels are at real risk of complete",
            " separation (absurdly large HR, CI spanning orders of magnitude). Consider collapsing",
            " categories before trusting the fit. ***")
  }
  invisible(tab)
}

## Reference (primary) Stage B HR, carried through every comparison below.
primary_n <- nrow(imputed_cohorts[[1]])
primary_stageB <- cox_pooled_summary %>%
  filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
  transmute(sensitivity = "Primary analysis (Model 3, MI-pooled)", n = primary_n, HR, HR_lower, HR_upper)

sensitivity_results <- list("primary" = primary_stageB)

## ============================================================
## 1. Alternative NT-proBNP threshold
## ============================================================
## Primary Stage B definition: age/sex-specific 90th percentile among
## the non-HF reference population. Alternative here: a FIXED clinical
## screening threshold (125 pg/mL)
alt_threshold <- 125  # pg/mL

cohort_alt_threshold <- cohort %>%
  mutate(
    nt_probnp_elevated_alt = as.numeric(ssbnp > alt_threshold),
    hf_continuum_level_alt = case_when(
      mcq160b == 1 ~ 3L,
      is.na(mcq160b) ~ NA_integer_,
      nt_probnp_elevated_alt == 1 ~ 2L,
      is.na(nt_probnp_elevated_alt) & is.na(stage_a_risk_factor) ~ NA_integer_,
      stage_a_risk_factor == 1 ~ 1L,
      TRUE ~ 0L
    ),
    hf_continuum_label_alt = factor(
      hf_continuum_level_alt, levels = 0:3,
      labels = c("No apparent HF risk", "Stage A (at risk)",
                 "Stage B (pre-HF, biomarker)", "Stage C (self-reported HF)")
    )
  ) %>%
  filter(!is.na(hf_continuum_level_alt))

message("\n---- 1. Alternative NT-proBNP threshold (>", alt_threshold, " pg/mL) ----")
message("N = ", nrow(cohort_alt_threshold), " (primary N = ", nrow(cohort), ")")

imputed_cohorts_alt_threshold <- purrr::map(imputed_cohorts, function(d) {
  d %>%
    select(-any_of(c("hf_continuum_level", "hf_continuum_label"))) %>%
    inner_join(
      cohort_alt_threshold %>% select(seqn, hf_continuum_level_alt, hf_continuum_label_alt),
      by = "seqn"
    ) %>%
    rename(hf_continuum_level = hf_continuum_level_alt, hf_continuum_label = hf_continuum_label_alt) %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
})

check_cell_counts(imputed_cohorts_alt_threshold[[1]], label = "Alt. NT-proBNP threshold")

alt_threshold_fit <- fit_pooled_cox_sens(
  cox_formula_m3, imputed_cohorts_alt_threshold,
  paste0("Alt. NT-proBNP threshold (>", alt_threshold, " pg/mL, fixed)")
)
message("  Full covariate table (all terms, not just Stage B):")
print(alt_threshold_fit %>% select(term, HR, HR_lower, HR_upper))

sensitivity_results[["alt_threshold"]] <- alt_threshold_fit %>%
  filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
  transmute(sensitivity, n = nrow(imputed_cohorts_alt_threshold[[1]]), HR, HR_lower, HR_upper)
print(sensitivity_results[["alt_threshold"]])

## ============================================================
## 2. Excluding advanced CKD
## ============================================================
## recompute the NT-proBNP percentile threshold EXCLUDING advanced CKD (eGFR <30, KDIGO stage 4-5)
advanced_ckd_egfr_cutoff <- 30

cohort_excl_ckd <- cohort %>%
  filter(is.na(egfr) | egfr >= advanced_ckd_egfr_cutoff) %>%
  group_by(age_band, riagendr) %>%
  mutate(
    nt_probnp_pctile_ckdexcl = {
      if (is.na(age_band[1])) {
        rep(NA_real_, n())
      } else {
        ref <- ssbnp[!is.na(mcq160b) & mcq160b == 2 & !is.na(ssbnp)]
        if (length(ref) < 10) rep(NA_real_, n()) else ecdf(ref)(ssbnp) * 100
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    nt_probnp_elevated_ckdexcl = as.numeric(nt_probnp_pctile_ckdexcl > 90),
    hf_continuum_level_ckdexcl = case_when(
      mcq160b == 1 ~ 3L,
      is.na(mcq160b) ~ NA_integer_,
      nt_probnp_elevated_ckdexcl == 1 ~ 2L,
      is.na(nt_probnp_elevated_ckdexcl) & is.na(stage_a_risk_factor) ~ NA_integer_,
      stage_a_risk_factor == 1 ~ 1L,
      TRUE ~ 0L
    ),
    hf_continuum_label_ckdexcl = factor(
      hf_continuum_level_ckdexcl, levels = 0:3,
      labels = c("No apparent HF risk", "Stage A (at risk)",
                 "Stage B (pre-HF, biomarker)", "Stage C (self-reported HF)")
    )
  ) %>%
  filter(!is.na(hf_continuum_level_ckdexcl))

message("\n---- 2. Excluding advanced CKD (eGFR <", advanced_ckd_egfr_cutoff, ") ----")
message("N = ", nrow(cohort_excl_ckd), " (primary N = ", nrow(cohort), "; ",
        sum(cohort$egfr < advanced_ckd_egfr_cutoff, na.rm = TRUE), " excluded for advanced CKD)")

imputed_cohorts_excl_ckd <- purrr::map(imputed_cohorts, function(d) {
  d %>%
    select(-any_of(c("hf_continuum_level", "hf_continuum_label"))) %>%
    inner_join(
      cohort_excl_ckd %>% select(seqn, hf_continuum_level_ckdexcl, hf_continuum_label_ckdexcl),
      by = "seqn"
    ) %>%
    rename(hf_continuum_level = hf_continuum_level_ckdexcl, hf_continuum_label = hf_continuum_label_ckdexcl) %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
})

check_cell_counts(imputed_cohorts_excl_ckd[[1]], label = "Excluding advanced CKD")

excl_ckd_fit <- fit_pooled_cox_sens(
  cox_formula_m3, imputed_cohorts_excl_ckd,
  paste0("Excluding advanced CKD (eGFR<", advanced_ckd_egfr_cutoff, ")")
)
message("  Full covariate table (all terms, not just Stage B):")
print(excl_ckd_fit %>% select(term, HR, HR_lower, HR_upper))

sensitivity_results[["excl_ckd"]] <- excl_ckd_fit %>%
  filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
  transmute(sensitivity, n = nrow(imputed_cohorts_excl_ckd[[1]]), HR, HR_lower, HR_upper)
print(sensitivity_results[["excl_ckd"]])

## ============================================================
## 3. Excluding early deaths
## ============================================================
## Excludes participants who died within the first 6 months of follow-up
early_death_window_months <- 6
n_early_deaths <- sum(cohort$mortstat == 1 & cohort$permth_exm < early_death_window_months, na.rm = TRUE)

message("\n---- 3. Excluding early deaths (<", early_death_window_months, " months) ----")
message(n_early_deaths, " early deaths excluded")

imputed_cohorts_excl_early <- purrr::map(imputed_cohorts, function(d) {
  d %>%
    filter(!(mortstat == 1 & permth_exm < early_death_window_months)) %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
})

check_cell_counts(imputed_cohorts_excl_early[[1]], label = "Excluding early deaths")

excl_early_fit <- fit_pooled_cox_sens(
  cox_formula_m3, imputed_cohorts_excl_early,
  paste0("Excluding early deaths (<", early_death_window_months, "mo)")
)
message("  Full covariate table (all terms, not just Stage B):")
print(excl_early_fit %>% select(term, HR, HR_lower, HR_upper))

sensitivity_results[["excl_early"]] <- excl_early_fit %>%
  filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
  transmute(sensitivity, n = nrow(imputed_cohorts_excl_early[[1]]), HR, HR_lower, HR_upper)
print(sensitivity_results[["excl_early"]])

## ============================================================
## 4. Obesity-stratified analysis
## ============================================================
## Fits Model 3 separately within obese (BMI>=30) and non-obese strata.

cox_formula_m3_noobesity <- update(cox_formula_m3, . ~ . - flag_obesity)

collapse_ref_for_obesity <- function(lbl) {
  factor(
    ifelse(as.character(lbl) %in% c("No apparent HF risk", "Stage A (at risk)"),
           "No apparent risk / Stage A", as.character(lbl)),
    levels = c("No apparent risk / Stage A", "Stage B (pre-HF, biomarker)", "Stage C (self-reported HF)")
  )
}

for (obesity_stratum in c(0, 1)) {
  stratum_label <- if (obesity_stratum == 1) "Obesity-stratified: Obese (BMI>=30)" else "Obesity-stratified: Non-obese (BMI<30)"
  
  imputed_cohorts_stratum <- purrr::map(imputed_cohorts, function(d) {
    d %>%
      filter(as.character(flag_obesity) == as.character(obesity_stratum)) %>%
      mutate(hf_continuum_label = collapse_ref_for_obesity(hf_continuum_label)) %>%
      mutate(across(c(riagendr, ridreth1, dmdeduc2,
                      flag_hypertension, flag_diabetes,
                      flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
  })
  
  message("\n---- 4. ", stratum_label, " (N = ", nrow(imputed_cohorts_stratum[[1]]), ") ----")
  message("  NOTE: reference category collapsed to 'No apparent risk / Stage A' -- see comment above",
          " this section for why.")
  check_cell_counts(imputed_cohorts_stratum[[1]], label = stratum_label, min_events = 10)
  
  stratum_fit <- fit_pooled_cox_sens(cox_formula_m3_noobesity, imputed_cohorts_stratum, stratum_label)
  message("  Full covariate table (all terms, not just Stage B):")
  print(stratum_fit %>% select(term, HR, HR_lower, HR_upper))
  
  key <- paste0("obesity_", obesity_stratum)
  sensitivity_results[[key]] <- stratum_fit %>%
    filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
    transmute(sensitivity, n = nrow(imputed_cohorts_stratum[[1]]), HR, HR_lower, HR_upper)
  print(sensitivity_results[[key]])
}

## ============================================================
## 5. Complete-case analysis
## ============================================================
cox_m3_covariates <- c("ridageyr", "riagendr", "ridreth1", "dmdeduc2", "indfmpir", "hiq011",
                       "flag_hypertension", "flag_diabetes", "flag_obesity", "flag_ckd",
                       "flag_smoking", "flag_chd_mi", "flag_stroke")
complete_case_data <- cohort %>%
  filter(if_all(all_of(cox_m3_covariates), ~ !is.na(.x))) %>%
  mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                  flag_hypertension, flag_diabetes, flag_obesity,
                  flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))

message("\n---- 5. Complete-case N = ", nrow(complete_case_data),
        " (vs full cohort N = ", nrow(cohort), ") ----")
check_cell_counts(complete_case_data, label = "Complete-case")

cc_design <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~wtmec6yr,
                       nest = TRUE, data = complete_case_data)
cc_cox <- svycoxph(cox_formula_m3, design = cc_design)
cc_cox_summary <- summary(cc_cox)$coefficients %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  mutate(HR = exp(coef), HR_lower = exp(coef - 1.96 * `se(coef)`),
         HR_upper = exp(coef + 1.96 * `se(coef)`))
message("  Full covariate table (all terms, not just Stage B):")
print(cc_cox_summary %>% select(term, HR, HR_lower, HR_upper))

comparison <- cox_pooled_summary %>%
  select(term, HR_mi = HR, HR_mi_lower = HR_lower, HR_mi_upper = HR_upper) %>%
  left_join(
    cc_cox_summary %>% select(term, HR_cc = HR, HR_cc_lower = HR_lower, HR_cc_upper = HR_upper),
    by = "term"
  )
message("\n---- MI vs complete-case Cox model comparison (all terms) ----")
print(comparison)

sensitivity_results[["complete_case"]] <- cc_cox_summary %>%
  filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
  transmute(sensitivity = "Complete-case (no MI)", n = nrow(complete_case_data), HR, HR_lower, HR_upper)
print(sensitivity_results[["complete_case"]])

## ============================================================
## 6. Age-as-time-scale Cox model
## ============================================================
cox_formula_age_timescale <- update(
  cox_formula_m3, Surv(age_entry, age_exit, mortstat) ~ . - ridageyr
)

imputed_cohorts_age_ts <- purrr::map(imputed_cohorts, function(d) {
  d %>%
    mutate(age_entry = ridageyr, age_exit = ridageyr + permth_exm / 12) %>%
    mutate(across(c(riagendr, ridreth1, dmdeduc2, hf_continuum_label,
                    flag_hypertension, flag_diabetes, flag_obesity,
                    flag_ckd, flag_smoking, flag_chd_mi, flag_stroke, hiq011), as.factor))
})

message("\n---- 6. Age-as-time-scale Cox model (N = ", nrow(imputed_cohorts_age_ts[[1]]), ") ----")
check_cell_counts(imputed_cohorts_age_ts[[1]], label = "Age-as-time-scale")

age_ts_fit <- fit_pooled_cox_sens(cox_formula_age_timescale, imputed_cohorts_age_ts, "Age-as-time-scale")
message("  Full covariate table (all terms, not just Stage B):")
print(age_ts_fit %>% select(term, HR, HR_lower, HR_upper))

sensitivity_results[["age_timescale"]] <- age_ts_fit %>%
  filter(term == "hf_continuum_labelStage B (pre-HF, biomarker)") %>%
  transmute(sensitivity, n = nrow(imputed_cohorts_age_ts[[1]]), HR, HR_lower, HR_upper)
print(sensitivity_results[["age_timescale"]])


## ============================================================
## Table 5: Sensitivity Analyses
## ============================================================
row_lookup <- tibble::tribble(
  ~sensitivity,                                                                  ~row_no, ~analysis,                          ~specification,
  "Primary analysis (Model 3, MI-pooled)",                                       "Ref.",  "Primary analysis",                 "Model 3, MI-pooled (m=20)",
  paste0("Alt. NT-proBNP threshold (>", alt_threshold, " pg/mL, fixed)"),        "1",     "Alternative NT-proBNP threshold",  paste0("Fixed threshold, >", alt_threshold, " pg/mL"),
  paste0("Excluding advanced CKD (eGFR<", advanced_ckd_egfr_cutoff, ")"),        "2",     "Excluding advanced CKD",           paste0("eGFR <", advanced_ckd_egfr_cutoff, " mL/min/1.73m\u00b2 excluded"),
  paste0("Excluding early deaths (<", early_death_window_months, "mo)"),         "3",     "Excluding early deaths",           paste0("Deaths <", early_death_window_months, " months excluded"),
  "Obesity-stratified: Non-obese (BMI<30)",                                      "4a",    "Obesity-stratified",               "Non-obese (BMI <30 kg/m\u00b2)",
  "Obesity-stratified: Obese (BMI>=30)",                                         "4b",    "Obesity-stratified",               "Obese (BMI \u226530 kg/m\u00b2)",
  "Complete-case (no MI)",                                                       "5",     "Complete-case analysis",           "No multiple imputation",
  "Age-as-time-scale",                                                           "6",     "Age-as-time-scale Cox model",      "Attained age as the Cox time scale"
)

table5 <- bind_rows(sensitivity_results, .id = NULL) %>%
  mutate(HR_CI = sprintf("%.2f (%.2f\u2013%.2f)", HR, HR_lower, HR_upper)) %>%
  left_join(row_lookup, by = "sensitivity") %>%
  arrange(match(row_no, c("Ref.", "1", "2", "3", "4a", "4b", "5", "6"))) %>%
  select(row_no, analysis, specification, n, HR_CI, HR, HR_lower, HR_upper)


message("\n============================================================")
message("Table 5. Sensitivity Analyses -- Stage B/pre-HF HR across specifications")
message("(Row 7, broader 1999-2018 replication, NOT included -- not implemented, see script comments)")
message("============================================================")
print(table5 %>% select(row_no, analysis, specification, n, HR_CI))

## ---- Excel output
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "Table 5")

display_tbl <- table5 %>%
  transmute(`#` = row_no, Analysis = analysis, Specification = specification,
            N = n, `HR (95% CI)` = HR_CI)

writeData(wb, "Table 5", "Table 5. Sensitivity Analyses", startRow = 1, startCol = 1)
addStyle(wb, "Table 5", createStyle(fontSize = 13, textDecoration = "bold"), rows = 1, cols = 1)

writeData(wb, "Table 5", display_tbl, startRow = 3, startCol = 1, headerStyle = createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", border = "TopBottomLeftRight"
))

body_style <- createStyle(border = "TopBottomLeftRight", borderColour = "#B7B7B7", valign = "center")
addStyle(wb, "Table 5", body_style, rows = 4:(3 + nrow(display_tbl)),
         cols = 1:ncol(display_tbl), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Table 5", createStyle(halign = "center"), rows = 4:(3 + nrow(display_tbl)),
         cols = c(1, 4, 5), gridExpand = TRUE, stack = TRUE)

addStyle(wb, "Table 5", createStyle(fgFill = "#D9EAD3", textDecoration = "bold"),
         rows = 4, cols = 1:ncol(display_tbl), gridExpand = TRUE, stack = TRUE)

setColWidths(wb, "Table 5", cols = 1:ncol(display_tbl), widths = c(6, 26, 34, 10, 20))

addWorksheet(wb, "Full Numeric Detail")
writeData(wb, "Full Numeric Detail", table5 %>% select(row_no, analysis, specification, n, HR, HR_lower, HR_upper))
addWorksheet(wb, "MI vs Complete-Case (all terms)")
writeData(wb, "MI vs Complete-Case (all terms)", comparison)

saveWorkbook(wb, "data/HF_continuum_Table5_sensitivity_analyses.xlsx", overwrite = TRUE)
message("\nSaved: data/HF_continuum_Table5_sensitivity_analyses.xlsx (styled, with footnotes)")