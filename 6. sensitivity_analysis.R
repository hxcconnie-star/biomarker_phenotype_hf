## ============================================================
## Biomarker-Informed Phenotyping of the HF Continuum
## Sensitivity Analyses (Results §7, Table 5)
## ============================================================
## PREREQUISITE: run hf_continuum_phenotype_analysis.R FIRST, in the
## SAME R session. This script reuses (without rebuilding):
##   cohort, imputed_cohorts, imputed_cohorts_conventional,
##   cox_formula_m3, cox_pooled_summary (= cox_m3$summary, the primary
##   MI-pooled Model 3 result used as the reference row below)
## None of the checks below touch the imputation itself (mice()), so
## there's no need to re-run MI -- only the cohort/phenotype/model
## specification changes per row.
##
## Covers protocol Table 5 rows:
##   1. Alternative NT-proBNP threshold
##   2. Excluding advanced CKD
##   3. Excluding early deaths
##   4. Obesity-stratified analysis
##   5. Complete-case analysis (moved here from the main script's old §10)
##   6. Age-as-time-scale Cox model
##   7. Broader NHANES 1999-2018 replication -- NOT IMPLEMENTED, see the
##      note near the bottom of this script for exactly why, and what a
##      partial version would require.

library(dplyr)
library(purrr)
library(tibble)
library(survival)
library(survey)
library(mitools)
library(writexl)

stopifnot(exists("cohort"), exists("imputed_cohorts"), exists("cox_formula_m3"),
          exists("cox_pooled_summary"))

## ---- Shared helper 1: fit Model 3 (or a variant) pooled across MI ----
## Self-contained here (doesn't depend on fit_pooled_cox() from the main
## script still being in memory) -- mirrors its logic exactly.
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

## ---- Shared helper 2: cell-count / event-count diagnostic ------------
## Checks the number of PEOPLE and EVENTS (deaths) in every phenotype
## level of a given dataset. A phenotype level with very few events (or
## zero) acting as the reference category is exactly what produces
## complete/quasi-complete separation -- absurd HRs like 148,000 with a
## CI spanning several orders of magnitude, a numerical artifact, not a
## real effect. Run this BEFORE fitting any stratified/subset analysis
## below, not just after something looks wrong.
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
## screening threshold (125 pg/mL), not age/sex-stratified -- commonly
## cited (e.g. ESC guidelines) for ruling out chronic HF in a non-acute/
## outpatient setting.
## ASSUMPTION: 125 pg/mL is my choice for "a clinically used threshold"
## since the protocol doesn't specify an exact number -- confirm this is
## what you intend, or tell me the value you'd rather use.
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
## Protocol's Stage B secondary definition #3: recompute the NT-proBNP
## percentile threshold EXCLUDING advanced CKD (eGFR <30, KDIGO stage
## 4-5) from the reference population -- NT-proBNP is renally cleared
## and can be elevated by reduced kidney function independent of cardiac
## status. Advanced-CKD participants are also excluded from the analytic
## sample for this specific check (testing whether the phenotype's
## prognostic value holds outside this confounded group).
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
## Excludes participants who died within the first 6 months of
## follow-up -- addresses reverse causation (someone already terminally
## ill at baseline could have elevated NT-proBNP or a recent HF
## diagnosis BECAUSE they were dying, not as a genuine prospective risk
## marker).
## ASSUMPTION: 6-month window is my choice; the protocol doesn't specify
## an exact cutoff -- 6 and 12 months are both common in the literature.
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
## Obesity is known to LOWER circulating NT-proBNP independent of
## cardiac status (adipose tissue clearance effect), so stratifying
## checks whether the phenotype's prognostic value is consistent across
## this potential effect-modifier. flag_obesity is dropped from the
## formula WITHIN each stratum (it's constant there, not a variable).
##
## CONFIRMED ISSUE, FIXED: flag_obesity is itself one of the 10
## components of the Stage A risk-factor composite -- so within the
## OBESE stratum specifically, "No apparent HF risk" (which by
## definition requires ZERO risk factors) is a near-empty reference
## category (confirmed: 5 people, ZERO deaths in this cohort's obese
## subgroup). A reference category with 0 events produces complete
## separation -- the fitted Stage B HR came out as ~148,000 with a CI
## spanning multiple orders of magnitude, a numerical artifact, not a
## real effect estimate. NOT an algorithm problem; a genuinely near-empty
## cell created by this subgroup's structure.
##
## Fix: within the obesity-stratified analysis ONLY, collapse "No
## apparent HF risk" and "Stage A (at risk)" into one combined reference
## category ("No apparent risk / Stage A"). This is defensible to report
## as-is: within the obese stratum, "at risk" is close to the default
## state anyway (since obesity itself is a risk-factor trigger), so the
## meaningful contrast becomes "no overt/pre-clinical disease" (combined
## reference) vs "Stage B" vs "Stage C" -- both strata now use the SAME
## 3-level structure for a clean side-by-side comparison, though the
## reference category's meaning differs slightly from the primary
## analysis's 4-level structure. Disclose this explicitly in the
## manuscript wherever this row is reported (see notes printed below).
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
## Drop anyone missing ANY Cox-model covariate, then fit the SAME model
## once (no MI). Compare against the primary MI-pooled result -- a large
## discrepancy would suggest the missing covariates aren't missing-at-
## random, and the MI/complete-case choice matters for the conclusions.
##
## Filters completeness on ALL 13 Model 3 covariates (not a subset) --
## confirmed necessary earlier in this project: filtering on a subset
## leaves residual NA that coxph/svycoxph silently drop during fitting,
## while the survey design object (built from the pre-drop data) still
## "thinks" it has the full N, causing a dimension mismatch inside
## svycoxph's internal variance calculation.
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
## Uses attained AGE (not follow-up time since baseline) as the
## underlying time scale -- a well-established alternative Cox
## parameterization that can better control confounding by age when age
## is a strong risk factor (as it is here; Korn, Graubard & Midthune
## 1997). Requires the "counting process" Surv(start, stop, event)
## format: start = age at baseline exam, stop = age at event/censoring.
## Age itself is DROPPED from the covariate list (it's the time scale
## now, not a covariate).
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
## 7. Broader NHANES 1999-2018 replication -- NOT IMPLEMENTED
## ============================================================
## NT-proBNP (SSBNP_A) was measured ONLY in the 1999-2004 surplus-sera
## sub-study -- no NHANES cycle after 2004 re-measured it. A true
## replication of the 4-level HF-continuum phenotype (including Stage B)
## is therefore not possible in 2005-2018 data under any
## operationalization.
##
## A PARTIAL replication IS feasible: extending the cohort through 2018
## using only the 3-level "conventional" phenotype (No apparent risk /
## Stage A / Stage C -- matching imputed_cohorts_conventional's
## construction in the main script) to check whether the conventional-
## risk-factor mortality gradient replicates in the broader window. This
## requires a SEPARATE extraction pipeline: nhanes_hf_continuum_extraction.R
## currently only pulls 1999-2004 cycles, and NHANES variable names/
## availability shift across cycles (the same kind of quirks documented
## in HF_Continuum_Codebook.xlsx's "Known Data Quirks" sheet) -- this is
## substantial enough that it's flagged here rather than attempted
## inline. Tell me if you want this built out as its own extraction +
## analysis script.
message("\n---- 7. Broader NHANES 1999-2018 replication: NOT IMPLEMENTED ----")
message("NT-proBNP unavailable outside 1999-2004 -- see the comment above this message",
        " for what a partial (conventional-phenotype-only) replication would require.")

## ============================================================
## Table 5: Sensitivity Analyses (assembled, publication-ready layout)
## ============================================================
## Row structure follows the protocol's own Table 5 numbering (1-7).
## "Analysis" is the clean row label; "Specification" carries the
## technical detail (threshold value, cutoff, etc.) that used to be
## crammed into one long parenthetical string -- splitting these out is
## what actually makes the table readable at a glance.
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

## Footnotes -- what to caveat when writing each row into the manuscript.
## Kept as a separate lookup (printed below the table / as an Excel
## footnote block) rather than a table column, matching how most
## journals present a numbered Table 5 with lettered footnotes.
table5_footnotes <- c(
  "1" = "Threshold value (125 pg/mL) is an assumed clinical screening cutoff, not protocol-specified -- confirm before finalizing.",
  "2" = "Excludes advanced-CKD participants entirely; smaller N than primary -- check that any CI width difference is attributable to sample size, not just effect-size change.",
  "3" = "6-month window is an assumed cutoff, not protocol-specified.",
  "4a" = "Reference category collapsed ('No apparent risk' + 'Stage A' merged) for a clean side-by-side comparison with 4b -- see footnote 4b.",
  "4b" = "Reference category collapsed ('No apparent risk' + 'Stage A' merged) because obesity is itself a Stage A risk-factor component, making an unmerged 'obese + no apparent risk' cell near-empty (5 people, 0 deaths) and producing complete separation (unmerged HR ~148,000) before this fix.",
  "5" = "No multiple imputation -- drops anyone missing any Model 3 covariate; compare N against the primary analysis to gauge how much data this discards.",
  "6" = "Age is the underlying time scale (not a covariate); HR magnitude is not directly comparable to the primary analysis's per-year age adjustment, though direction/significance are what matters for this comparison."
)

message("\n============================================================")
message("Table 5. Sensitivity Analyses -- Stage B/pre-HF HR across specifications")
message("(Row 7, broader 1999-2018 replication, NOT included -- not implemented, see script comments)")
message("============================================================")
print(table5 %>% select(row_no, analysis, specification, n, HR_CI))
message("\nFootnotes:")
for (i in seq_along(table5_footnotes)) {
  message(names(table5_footnotes)[i], ". ", table5_footnotes[i])
}

## ---- Publication-ready Excel output (openxlsx -- supports styling; ---
## writexl does not) --------------------------------------------------
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

## Shade the primary-analysis (reference) row so it visually anchors the
## comparisons, matching how the other tables in this project mark
## reference rows.
addStyle(wb, "Table 5", createStyle(fgFill = "#D9EAD3", textDecoration = "bold"),
         rows = 4, cols = 1:ncol(display_tbl), gridExpand = TRUE, stack = TRUE)

setColWidths(wb, "Table 5", cols = 1:ncol(display_tbl), widths = c(6, 26, 34, 10, 20))

## Footnotes block below the table.
footnote_start_row <- 3 + nrow(display_tbl) + 2
writeData(wb, "Table 5", "Footnotes:", startRow = footnote_start_row, startCol = 1)
addStyle(wb, "Table 5", createStyle(textDecoration = "bold"), rows = footnote_start_row, cols = 1)
for (i in seq_along(table5_footnotes)) {
  writeData(wb, "Table 5",
            paste0(names(table5_footnotes)[i], ". ", table5_footnotes[i]),
            startRow = footnote_start_row + i, startCol = 1)
}

## Second sheet: full numeric detail (HR/lower/upper separately) + MI vs
## complete-case comparison across ALL terms, for anyone who wants the
## raw numbers behind the formatted HR (95% CI) strings above.
addWorksheet(wb, "Full Numeric Detail")
writeData(wb, "Full Numeric Detail", table5 %>% select(row_no, analysis, specification, n, HR, HR_lower, HR_upper))
addWorksheet(wb, "MI vs Complete-Case (all terms)")
writeData(wb, "MI vs Complete-Case (all terms)", comparison)

saveWorkbook(wb, "data/HF_continuum_Table5_sensitivity_analyses.xlsx", overwrite = TRUE)
message("\nSaved: data/HF_continuum_Table5_sensitivity_analyses.xlsx (styled, with footnotes)")