## ============================================================
##   >>> MULTIPLE IMPUTATION INSERTION POINT (Methods §8.2) <<<
##
##   Everything above this line (Results §1-3) used OBSERVED data only --
##   phenotype assignment, cohort construction, and baseline description
##   never touch imputed values. Everything below this line (Results §4,
##   §7) uses the imputed covariates for adjusted models, per the
##   agreed pipeline order:
##     Results §1-3 (observed data)
##       -> MI (covariates only, never phenotype/outcome)
##       -> Results §4 (trajectory analysis, adjusted models use MI)
##       -> Results §7 (sensitivity: complete-case vs MI comparison)
## ============================================================
## ============================================================

## ---- 8. Multiple Imputation (Methods §8.2) --------------------------
## ============================================================
## Per protocol §8.2: impute covariates only. mortstat, permth_exm, and
## hf_continuum_level are NEVER imputed -- they enter the imputation
## model as complete, fully-observed predictors (this is what "include
## the outcome indicator in the imputation model" means; it prevents the
## imputation from being biased with respect to survival). The
## Nelson-Aalen cumulative hazard estimate is included instead of raw
## follow-up time, per White & Royston (2009) -- standard practice for
## MI feeding into a Cox model, and exactly what your protocol names
## explicitly ("Nelson-Aalen cumulative hazard").
##
## Imputed covariates: expanded to cover ALL variables needed across
## protocol §8.4's three nested Cox models (Model 4 excluded for the main
## phenotype exposure -- see section 9.2 below for why):
##   Model 1: age, sex, race/ethnicity
##   Model 2: + education, poverty-income ratio, insurance
##   Model 3: + hypertension, diabetes, obesity, CKD, smoking, CHD/MI/stroke
## ridageyr/riagendr/ridreth1 are complete (0% missing) so mice passes
## them through unchanged; the rest have real missingness.

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

## CONFIRMED BUG, ROOT-CAUSE FIX: mice's OWN auto-detection of which
## method to use per variable is NOT reliable enough to build the factor-
## conversion list from -- that's exactly what went wrong last time.
## flag_obesity is just as binary (0/1) as flag_hypertension/flag_diabetes/
## etc, but mice's dry-run auto-assigned it "pmm" instead of "logreg" (an
## inconsistency in mice's own heuristic, not something we can predict),
## so a detection step that only looks for existing logreg/polyreg/polr
## assignments never catches it -- confirmed: rerunning that "fix"
## reproduced the EXACT same 899/23/205/62/... failure counts, proving
## flag_obesity was never actually converted.
##
## Fix: stop asking mice what it thinks each variable is. Classify
## explicitly, from what we already know about each variable's true
## nature, and force factor conversion for every categorical one --
## regardless of what mice's own dry-run method-detection would guess.
continuous_mi_vars <- c("ridageyr", "indfmpir", "nelsonaalen")
categorical_mi_vars <- c("riagendr", "ridreth1", "dmdeduc2", "hiq011",
                         "flag_hypertension", "flag_diabetes", "flag_obesity",
                         "flag_ckd", "flag_smoking", "flag_chd_mi", "flag_stroke")
## (seqn/sdmvpsu/sdmvstra/wtmec6yr/mortstat/permth_exm/hf_continuum_level
## are all in never_impute below regardless, so their type doesn't matter.)

mi_data <- mi_data_raw %>% mutate(across(all_of(categorical_mi_vars), as.factor))
message("Forced factor conversion for: ", paste(categorical_mi_vars, collapse = ", "))

## Dry run to get mice's default method/predictor matrix, then override:
## never impute mortstat/permth_exm/hf_continuum_level/design vars/
## nelsonaalen (they're complete anyway, but this is explicit and safe).
imp_setup <- mice(mi_data, maxit = 0, printFlag = FALSE)
meth <- imp_setup$method
never_impute <- c("seqn", "sdmvpsu", "sdmvstra", "wtmec6yr",
                  "mortstat", "permth_exm", "hf_continuum_level", "nelsonaalen")
meth[never_impute] <- ""

pred <- imp_setup$predictorMatrix
## Don't use seqn as a predictor for anything (it's just an ID).
pred[, "seqn"] <- 0

## Extra safety check: confirm every categorical_mi_vars entry that
## ACTUALLY HAS MISSING VALUES got a categorical method (logreg/polyreg/
## polr). A variable with 0% missingness correctly gets method="" from
## mice (nothing to impute) -- that's expected, not a bug, so only check
## variables where there's something to actually impute.
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

## Verification: which variables will actually be imputed (non-empty
## method) vs excluded (method==""), and how much missingness each has.
## "" with 0% missing = correctly complete, nothing to impute anyway.
## "" with >0% missing would be a real problem (something that SHOULD be
## imputed got excluded) -- none of the never_impute vars should show
## that, since mortstat/permth_exm/hf_continuum_level are guaranteed
## complete by cohort construction.
verify_mi_setup <- tibble(
  variable = names(meth),
  method = meth,
  will_be_imputed = meth != "",
  pct_missing = round(100 * sapply(mi_data[names(meth)], function(x) mean(is.na(x))), 1)
)
message("\n---- MI setup verification ----")
print(verify_mi_setup, n = Inf)

n_imputations <- 20   # rule-of-thumb: at least the % missing in the most
# incomplete imputed variable (dmdeduc2 ~50% here,
# so 20 is a reasonably conservative choice, not a
# strict requirement -- adjust if runtime is an issue)

imp <- mice(mi_data, m = n_imputations, method = meth, predictorMatrix = pred,
            seed = 12345, printFlag = FALSE)

message("Multiple imputation complete: ", n_imputations, " imputed datasets.")

## Automated completeness check: for every variable mice was supposed to
## impute (non-empty method), confirm EVERY completed dataset (not just
## #1) has zero residual NAs. This is exactly the check that would have
## caught the flag_obesity/PMM donor-matching failure immediately.
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

## ---- 8b. Patch any residual NA left by mice (Methods §8.2 addendum) ---
## CONFIRMED via diagnosis: flag_obesity has ~899 residual NA in EVERY
## completed dataset, REGARDLESS of imputation method (logreg AND cart
## independently failed on the identical rows -- ruled out an algorithm-
## specific numerical issue). dmdeduc2's 23 and hiq011's 62 residual-NA
## rows are BOTH fully contained within this same ~899-person set,
## confirming these are people with enough SIMULTANEOUS missingness
## across several Model 3 covariates that no chained-equations method can
## generate a stable prediction/donor for them. This is a genuine data-
## sparsity limitation (~6% of the cohort), not a configuration bug.
##
## Leaving these ~6% with real missingness breaks svycoxph()'s internal
## variance calculation downstream (the survey design object and the
## actual fitted N disagree -- the same "replacement length" crash found
## and fixed for complete_case_data in section 10). Patched here with a
## simple, transparent, reportable rule: the modal (most common) value
## within the same age-decade x sex stratum, computed separately within
## EACH completed dataset. This preserves mice's own imputation
## uncertainty for the ~94% it DID handle successfully across the m
## datasets; only the unresolved ~6% get this deterministic fallback
## (identical across all m for these specific people) -- disclose this as
## a limitation in the manuscript's methods/missing-data section.
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
    ## CONFIRMED BUG FIX: ifelse() above always returns character, even
    ## for continuous variables. Left as character, indfmpir (continuous,
    ## ~hundreds of distinct decimal values) would get silently treated as
    ## a FACTOR with hundreds of levels by coxph()'s formula machinery --
    ## this is exactly what caused the "system is exactly singular" /
    ## "variable ...500" crash: the design matrix exploded to hundreds of
    ## spurious dummy columns instead of one continuous term. Convert
    ## continuous variables back to numeric explicitly; categorical ones
    ## are fine left as character (downstream code re-applies as.factor()
    ## explicitly before every model fit anyway).
    if (is_continuous) {
      df <- df %>% mutate(!!v := as.numeric(.data[[v]]))
    }
  }
  df %>% select(-.age_grp)
}

vars_to_patch <- imputed_var_names[residual_na_check > 0]

## Merge each imputed dataset's covariates back onto the FULL cohort
## (which still has hf_continuum_level, phenotype labels, egfr, etc. --
## mi_data only carried the minimal imputation-model variables).
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

## Extra type check: indfmpir (the only continuous variable in
## imputed_vars) MUST be numeric, not character -- confirmed root cause
## of the earlier "system is exactly singular" crash was exactly this
## (character-typed indfmpir silently treated as a several-hundred-level
## factor by coxph's formula machinery). This would have caught it
## immediately instead of a cryptic Lapack error deep inside svycoxph.
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