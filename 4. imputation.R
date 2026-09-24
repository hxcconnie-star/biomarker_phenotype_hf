## ============================================================
## MULTIPLE IMPUTATION
## ============================================================

## impute covariates only. mortstat, permth_exm, and
## hf_continuum_level are NEVER imputed 
##
## phenotype exposure:
##   Model 1: age, sex, race/ethnicity
##   Model 2: + education, poverty-income ratio, insurance
##   Model 3: + hypertension, diabetes, obesity, CKD, smoking, CHD/MI/stroke


mi_vars <- c("seqn", "sdmvpsu", "sdmvstra", "wtmec6yr",
             "mortstat", "permth_exm", "hf_continuum_level",
             ## Model 1
             "ridageyr", "riagendr", "ridreth1",
             ## Model 2 adds
             "dmdeduc2", "indfmpir", "hiq011",
             ## Model 3 adds
             "flag_hypertension", "flag_diabetes", "flag_obesity", "flag_ckd",
             "flag_smoking", "flag_chd_mi", "flag_stroke")
mi_data_raw <- cohort %>% select(all_of(mi_vars))
mi_data_raw$nelsonaalen <- mice::nelsonaalen(mi_data_raw, "permth_exm", "mortstat")

message("\n---- 8. Multiple imputation ----")
message("Variables being imputed (have missingness): ",
        paste(names(which(colSums(is.na(mi_data_raw)) > 0)), collapse = ", "))

continuous_mi_vars <- c("ridageyr", "indfmpir", "nelsonaalen")
categorical_mi_vars <- c("riagendr", "ridreth1", "dmdeduc2", "hiq011",
                         "flag_hypertension", "flag_diabetes", "flag_obesity",
                         "flag_ckd", "flag_smoking", "flag_chd_mi", "flag_stroke")


mi_data <- mi_data_raw %>% mutate(across(all_of(categorical_mi_vars), as.factor))
message("Forced factor conversion for: ", paste(categorical_mi_vars, collapse = ", "))

imp_setup <- mice(mi_data, maxit = 0, printFlag = FALSE)
meth <- imp_setup$method
never_impute <- c("seqn", "sdmvpsu", "sdmvstra", "wtmec6yr",
                  "mortstat", "permth_exm", "hf_continuum_level", "nelsonaalen")
meth[never_impute] <- ""

pred <- imp_setup$predictorMatrix

pred[, "seqn"] <- 0

categorical_mi_vars_with_missing <- categorical_mi_vars[
  sapply(mi_data[categorical_mi_vars], function(x) any(is.na(x)))
]
wrong_method <- categorical_mi_vars_with_missing[
  !meth[categorical_mi_vars_with_missing] %in% c("logreg", "polyreg", "polr")
]
if (length(wrong_method) > 0) {
  message("\n*** STOP: these categorical variables (which DO have missing values) did NOT",
          " get a categorical method even after factor conversion -- do not proceed: ***")
  print(meth[wrong_method])
} else {
  message("\nAll categorical variables with missingness correctly assigned a categorical method",
          " (logreg/polyreg/polr). Variables with 0% missingness (", 
          paste(setdiff(categorical_mi_vars, categorical_mi_vars_with_missing), collapse = ", "),
          ") correctly show method=\"\" since there's nothing to impute.")
}

verify_mi_setup <- tibble(
  variable = names(meth),
  method = meth,
  will_be_imputed = meth != "",
  pct_missing = round(100 * sapply(mi_data[names(meth)], function(x) mean(is.na(x))), 1)
)
message("\n---- MI setup verification ----")
print(verify_mi_setup, n = Inf)

n_imputations <- 20

imp <- mice(mi_data, m = n_imputations, method = meth, predictorMatrix = pred,
            seed = 12345, printFlag = FALSE)

message("Multiple imputation complete: ", n_imputations, " imputed datasets.")


imputed_var_names <- names(meth)[meth != ""]
residual_na_check <- sapply(imputed_var_names, function(v) {
  max(sapply(1:n_imputations, function(i) sum(is.na(complete(imp, i)[[v]]))))
})
if (any(residual_na_check > 0)) {
  message("\n*** WARNING: these variables still have residual NA in at least one",
          " completed dataset -- patching below before proceeding: ***")
  print(residual_na_check[residual_na_check > 0])
} else {
  message("Completeness check passed: no residual NA in any of the ", n_imputations,
          " completed datasets, for any imputed variable.")
}

## ---- 8b. Patch any residual NA left by mice ---
patch_residual_na <- function(df, vars) {
  df <- df %>% mutate(.age_grp = cut(ridageyr, breaks = c(19, 39, 59, 120),
                                     labels = c("20-39", "40-59", "60+")))
  for (v in vars) {
    if (sum(is.na(df[[v]])) == 0) next
    is_continuous <- v %in% continuous_mi_vars   # e.g. indfmpir -- median, not mode
    if (is_continuous) {
      fallback_by_group <- df %>%
        filter(!is.na(.data[[v]])) %>%
        group_by(.age_grp, riagendr) %>%
        summarise(.fallback = median(as.numeric(as.character(.data[[v]])), na.rm = TRUE), .groups = "drop")
    } else {
      fallback_by_group <- df %>%
        filter(!is.na(.data[[v]])) %>%
        group_by(.age_grp, riagendr) %>%
        summarise(.fallback = names(sort(table(as.character(.data[[v]])), decreasing = TRUE))[1],
                  .groups = "drop")
    }
    df <- df %>%
      left_join(fallback_by_group, by = c(".age_grp", "riagendr")) %>%
      mutate(!!v := ifelse(is.na(.data[[v]]), .fallback, as.character(.data[[v]]))) %>%
      select(-.fallback)

    if (is_continuous) {
      df <- df %>% mutate(!!v := as.numeric(.data[[v]]))
    }
  }
  df %>% select(-.age_grp)
}

vars_to_patch <- imputed_var_names[residual_na_check > 0]

## Merge each imputed dataset's covariates back onto the FULL cohort
imputed_vars <- c("dmdeduc2", "indfmpir", "hiq011", "flag_hypertension",
                  "flag_diabetes", "flag_obesity", "flag_ckd",
                  "flag_smoking", "flag_chd_mi", "flag_stroke")
imputed_cohorts <- purrr::map(1:n_imputations, function(i) {
  completed <- complete(imp, i)
  if (length(vars_to_patch) > 0) completed <- patch_residual_na(completed, vars_to_patch)
  cohort %>%
    select(-any_of(imputed_vars)) %>%   # drop originals, replace with imputed
    left_join(completed %>% select(seqn, all_of(imputed_vars)), by = "seqn")
})

## Final verification: confirm the patch actually achieved 100%
## completeness -- do NOT proceed to any Cox model below if this fails.
final_completeness_check <- sapply(imputed_vars, function(v) {
  max(sapply(imputed_cohorts, function(d) sum(is.na(d[[v]]))))
})
message("\n---- Final completeness check (after residual-NA patch) ----")
print(final_completeness_check)
if (any(final_completeness_check > 0)) {
  message("\n*** STOP: residual NA remains even after patching -- do NOT proceed",
          " to the Cox models below until this is resolved: ***")
  print(final_completeness_check[final_completeness_check > 0])
} else {
  message("All imputed variables are 100% complete across all ", n_imputations,
          " imputed cohorts. Safe to proceed.")
}

indfmpir_class <- sapply(imputed_cohorts, function(d) class(d$indfmpir)[1])
if (any(indfmpir_class != "numeric")) {
  message("\n*** STOP: indfmpir is not numeric in at least one imputed cohort",
          " (found: ", paste(unique(indfmpir_class), collapse = ", "), ") -- ",
          "this WILL cause a singular-matrix crash in the Cox models. Fix before proceeding. ***")
} else {
  message("indfmpir confirmed numeric in all imputed cohorts.")
}


## -------------- save imputation ----------------
saveRDS(
  list(
    imputed_cohorts = imputed_cohorts,
    cohort = cohort,
    # cox_formula = cox_formula,
    # cox_formula_m3 = cox_formula_m3,
    n_imputations = n_imputations
  ),
  file = "data/HF_continuum_imputed_cohorts.rds"
)
message("Saved: data/HF_continuum_imputed_cohorts.rds")