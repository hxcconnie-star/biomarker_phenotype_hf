## ============================================================
## Data Cleaning: between Extraction and Analysis
## ============================================================
##
## Input:  data/HF_continuum_analytic_dataset.xlsx  (raw output of the
##         extraction script, nhanes_hf_continuum_extraction.R)
## Output: data/HF_continuum_cleaned_dataset.xlsx    (feed this into
##         hf_continuum_phenotype_analysis.R instead of the raw file)
##
##
## 2.1 variable profiling (data-driven, all columns) -> 2.2 harmonize
## categorical text (+ residual check for anything missed) -> 2.3 range/
## plausibility checks -> 2.4 logical consistency checks -> 2.5 duplicate
## SEQN check -> 2.6 missingness summary -> save.

## ---- 0. Packages -------------------------------------------------
pkgs <- c("dplyr", "purrr", "readxl", "writexl", "tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(purrr)
library(readxl)
library(writexl)
library(tibble)

## ---- 1. Load raw extracted data --------------------------------------
hf_continuum_data <- readxl::read_excel("data/HF_continuum_analytic_dataset.xlsx")
message("Loaded raw extracted data: ", nrow(hf_continuum_data), " rows x ",
        ncol(hf_continuum_data), " cols")

names(hf_continuum_data) <- tolower(names(hf_continuum_data))

valid_pct_before <- sapply(hf_continuum_data, function(x) round(100 * mean(!is.na(x)), 2))

## ---- 2.1 Variable profiling (first 50 rows) --------------------------
preview <- head(hf_continuum_data, 50)

profile_variable <- function(x, varname) {
  x_nomiss <- x[!is.na(x)]
  tibble(
    variable = varname,
    class = class(x)[1],
    n_distinct = dplyr::n_distinct(x_nomiss),
    pct_missing_in_preview = round(100 * mean(is.na(x)), 1),
    sample_values = paste(head(unique(x_nomiss), 8), collapse = " | ")
  )
}

profile_report <- purrr::imap_dfr(preview, profile_variable)
message("\n---- 2.1 Variable profile (class, distinct values, samples -- first 50 rows) ----")
print(profile_report, n = Inf, width = Inf)
message("\nReview the table above: any 'character' class column with sample values that",
        " look like Yes/No/category text (not plain numbers) needs a harmonizer in 2.2.",
        " The list below reflects what's been confirmed by this method so far -- if the",
        " profile shows something new/different, that's the signal to add/adjust one.")

## ---- 2.1b Full untruncated text for low-cardinality character columns --
low_card_char_cols <- profile_report %>%
  filter(class == "character", n_distinct > 0, n_distinct <= 10) %>%
  pull(variable)

message("\n---- 2.1b Full distinct values for low-cardinality character columns ----")
for (v in low_card_char_cols) {
  vals <- unique(na.omit(hf_continuum_data[[v]]))
  message("\n", v, " (", length(vals), " distinct value(s)):")
  for (val in vals) message("  - \"", val, "\"")
}

## ---- 2.2 Harmonize categorical text to clean numeric codes -----------
harmonize_yesno <- function(x, yes_pattern = "^1$|^yes$", no_pattern = "^2$|^no$") {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl(yes_pattern, x_chr)] <- 1
  out[grepl(no_pattern, x_chr)] <- 2
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_yesno(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_smq040 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|every day", x_chr)] <- 1
  out[grepl("^2$|some days", x_chr)] <- 2
  out[grepl("^3$|not at all", x_chr)] <- 3
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_smq040(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_riagendr <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|^male$", x_chr)] <- 1
  out[grepl("^2$|^female$", x_chr)] <- 2
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_riagendr(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_ridexprg <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|pregnant", x_chr)] <- 1
  out[grepl("^2$|not pregnant", x_chr)] <- 2
  out[grepl("^3$|cannot ascertain", x_chr)] <- 3
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_ridexprg(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_dmdeduc2 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|9th grade|less than 9", x_chr)] <- 1
  out[grepl("^2$|9-11th|11th grade", x_chr)] <- 2
  out[grepl("^3$|high school", x_chr)] <- 3
  out[grepl("^4$|some college|aa degree", x_chr)] <- 4
  out[grepl("^5$|college graduate", x_chr)] <- 5
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_dmdeduc2(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_huq030 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|^yes$", x_chr)] <- 1
  out[grepl("^2$|no place", x_chr)] <- 2
  out[grepl("^3$|more than one", x_chr)] <- 3
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_huq030(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_huq010 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|excellent", x_chr)] <- 1
  out[grepl("^3$|good", x_chr)] <- 3
  out[grepl("^2$|very good", x_chr)] <- 2
  out[grepl("^4$|fair", x_chr)] <- 4
  out[grepl("^5$|poor", x_chr)] <- 5
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_huq010(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_ridreth1 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|mexican american", x_chr)] <- 1
  out[grepl("^2$|other hispanic", x_chr)] <- 2
  out[grepl("^3$|non-hispanic white", x_chr)] <- 3
  out[grepl("^4$|non-hispanic black", x_chr)] <- 4
  out[grepl("^5$|other race", x_chr)] <- 5
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_ridreth1(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_paq180 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|not walk about", x_chr)] <- 1
  out[grepl("^2$|stand or walk", x_chr)] <- 2
  out[grepl("^3$|light load|climb stairs", x_chr)] <- 3
  out[grepl("^4$|heavy work|heavy loads", x_chr)] <- 4
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_paq180(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_ssbnpl <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^0$|within the detection", x_chr)] <- 0
  out[grepl("^1$|below lower", x_chr)] <- 1
  out[grepl("^2$|above upper", x_chr)] <- 2
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_ssbnpl(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

harmonize_sspris <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(x_chr))
  out[grepl("^1$|^pristine$", x_chr)] <- 1
  out[grepl("^0$|^non-pristine$", x_chr)] <- 0
  unclassified <- !is.na(x_chr) & is.na(out)
  if (any(unclassified)) {
    message("  harmonize_sspris(): ", sum(unclassified), " unclassified value(s), e.g.: ",
            paste(unique(x_chr[unclassified])[1:min(5, sum(unclassified))], collapse = ", "))
  }
  out
}

message("\n---- 2.2 Harmonizing categorical variables ----")
hf_continuum_data <- hf_continuum_data %>%
  mutate(
    hiq011  = harmonize_yesno(hiq011),
    kiq022  = harmonize_yesno(kiq022),
    across(c(mcq160b, mcq160c, mcq160d, mcq160e, mcq160f,
             diq010, bpq020, bpq050a, bpq080, smq020, pad200), harmonize_yesno),
    ## paq180's category order confirmed against CDC's own codebook pages
    ## (see harmonize_paq180 above) -- no longer a guess.
    paq180 = harmonize_paq180(paq180),
    smq040 = harmonize_smq040(smq040),
    riagendr = harmonize_riagendr(riagendr),
    ridexprg = harmonize_ridexprg(ridexprg),
    dmdeduc2 = harmonize_dmdeduc2(dmdeduc2),
    huq030   = harmonize_huq030(huq030),
    huq010   = harmonize_huq010(huq010),
    ridreth1 = harmonize_ridreth1(ridreth1),
    ssbnpl   = harmonize_ssbnpl(ssbnpl),
    sspris   = harmonize_sspris(sspris),
    fsdhh_num = suppressWarnings(as.numeric(fsdhh)),
    fsdhh = case_when(
      !is.na(fsdhh_num) ~ fsdhh_num,
      grepl("full", tolower(fsdhh)) ~ 1,
      grepl("marginal", tolower(fsdhh)) ~ 2,
      grepl("low", tolower(fsdhh)) & !grepl("very", tolower(fsdhh)) ~ 3,
      grepl("very low", tolower(fsdhh)) ~ 4,
      TRUE ~ NA_real_
    ),

    across(c(mcd180b, mcd180c, mcd180d, mcd180e, mcd180f,
             eligstat, mortstat), ~ suppressWarnings(as.numeric(.x))),

    across(c(ridageyr, bmxbmi, bmxwaist, indfmpir, ssbnp,
             bpxsy1, bpxsy2, bpxsy3, bpxsy4, bpxdi1, bpxdi2, bpxdi3, bpxdi4,
             lbxgh, lbxscr, lbxsal, urxuma, urxucr, lbxtr, lbdldl, lbxtc, lbdhdd,
             permth_int, permth_exm),
           ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  select(-fsdhh_num)

message("Post-harmonization check (should be clean 1/2, or 1-4 for fsdhh, 1-3 for smq040/ridexprg, 1-4 for paq180):")
for (v in c("hiq011", "kiq022", "fsdhh", "mcq160b", "riagendr", "ridexprg",
            "smq040", "dmdeduc2", "huq030", "huq010", "pad200", "paq180", "ssbnpl", "sspris", "ridreth1")) {
  if (v %in% names(hf_continuum_data)) print(table(hf_continuum_data[[v]], useNA = "always"))
}

## ---- 2.2b Residual check: any character column not yet handled? -------
still_character <- names(hf_continuum_data)[sapply(hf_continuum_data, is.character)]
if (length(still_character) > 0) {
  message("\n*** ", length(still_character), " column(s) still character after harmonizing -- ",
          "review and add a harmonizer for these: ***")
  print(profile_report %>% filter(variable %in% still_character))
  message("Note: not everything flagged here is necessarily wrong -- e.g. ucod_leading",
          " (cause-of-death description) is legitimately text, not a code to harmonize.",
          " Check the sample_values column above to judge which ones actually need fixing.")
} else {
  message("\nNo character columns remain -- everything categorical has been harmonized.")
}

## ---- 2.3 Range / plausibility checks for continuous variables --------
flag_out_of_range <- function(x, lower, upper, varname) {
  bad <- !is.na(x) & (x < lower | x > upper)
  if (any(bad)) {
    message("  ", varname, ": ", sum(bad), " value(s) outside plausible range [",
            lower, ", ", upper, "] -- set to NA (observed: ",
            paste(round(range(x[bad], na.rm = TRUE), 1), collapse = " to "), ")")
    x[bad] <- NA_real_
  }
  x
}


range_bounds <- list(
  ridageyr = c(0, 85),      # FIXED: was (20, 85) -- that's a cohort criterion, not a plausibility bound
  bmxbmi   = c(8, 95),     # widened slightly (was 10, 90) as a precaution against the same age-assumption issue
  bmxwaist = c(15, 250),   # FIXED: was (40, 200) -- same bug as ridageyr, assumed adult-only population. The raw data spans all ages; a child's genuine waist circumference can be well under 40cm.
  bpxsy1 = c(60, 300), bpxsy2 = c(60, 300), bpxsy3 = c(60, 300), bpxsy4 = c(60, 300),
  bpxdi1 = c(0, 180), bpxdi2 = c(0, 180), bpxdi3 = c(0, 180), bpxdi4 = c(0, 180),  # WIDENED: DBP=0 is a legitimate finding
  lbxgh    = c(2, 20),
  lbxscr   = c(0.1, 20),
  lbxsal   = c(1, 6),
  urxuma   = c(0, 20000),
  urxucr   = c(1, 2000),    # WIDENED: concentrated urine samples can genuinely exceed 500
  lbxtr    = c(10, 5000),
  lbdldl   = c(10, 1000),   # WIDENED: familial hypercholesterolemia can genuinely exceed 500
  lbxtc    = c(50, 800),
  lbdhdd   = c(5, 250),     # WIDENED: rare genetic conditions can genuinely exceed 150
  ssbnp    = c(0, 70000),
  indfmpir = c(0, 5),
  mcd180b = c(0, 85), mcd180c = c(0, 85), mcd180d = c(0, 85),
  ## NOTE: profiling confirmed mcd180c-f contain NHANES sentinel codes like
  ## 99999 ("Don't know") and 77777 ("Refused") mixed in with real ages --
  ## these get correctly caught and nulled out here since they're wildly
  ## outside the [0, 85] bound, same as any other implausible value. No
  ## special-case code needed; flagging this so it's not a surprise when
  ## the "value(s) outside plausible range" messages below mention them.
  mcd180e = c(0, 85), mcd180f = c(0, 85),
  permth_int = c(0, 300), permth_exm = c(0, 300)
)

message("\n---- 2.3 Range/plausibility checks ----")
for (v in names(range_bounds)) {
  if (v %in% names(hf_continuum_data)) {
    hf_continuum_data[[v]] <- flag_out_of_range(
      hf_continuum_data[[v]], range_bounds[[v]][1], range_bounds[[v]][2], v
    )
  }
}

## ---- 2.4 Logical consistency checks -----------------------------------
message("\n---- 2.4 Logical consistency checks ----")

## Age at diagnosis can't exceed current age.
check_age_consistency <- function(data, dx_var, age_var = "ridageyr") {
  bad <- !is.na(data[[dx_var]]) & !is.na(data[[age_var]]) & data[[dx_var]] > data[[age_var]]
  if (any(bad)) {
    message("  ", dx_var, ": ", sum(bad), " row(s) where age at diagnosis > current age -- set to NA")
    data[[dx_var]][bad] <- NA_real_
  }
  data
}
for (v in c("mcd180b", "mcd180c", "mcd180d", "mcd180e", "mcd180f")) {
  hf_continuum_data <- check_age_consistency(hf_continuum_data, v)
}

bad_preg <- !is.na(hf_continuum_data$ridexprg) & hf_continuum_data$ridexprg == 1 &
  !is.na(hf_continuum_data$riagendr) & hf_continuum_data$riagendr == 1
if (any(bad_preg)) {
  message("  WARNING: ", sum(bad_preg), " row(s) coded pregnant (ridexprg==1) AND male",
          " (riagendr==1) -- inconsistent, review manually. NOT auto-corrected.")
} else {
  message("  Pregnancy x sex consistency: OK (0 conflicts)")
}

## ---- 2.5 Duplicate SEQN check -----------------------------------------
dup_seqn <- sum(duplicated(hf_continuum_data$seqn))
message("\n---- 2.5 Duplicate SEQN check: ", dup_seqn,
        ifelse(dup_seqn > 0, " -- INVESTIGATE, should be 0", " (OK)"))

## ---- 2.6 Missingness summary ------------------------------------------
message("\n---- 2.6 Missingness summary (% NA per variable, sorted descending) ----")
missing_summary <- sort(sapply(hf_continuum_data, function(x) round(100 * mean(is.na(x)), 1)),
                        decreasing = TRUE)
print(missing_summary)

## ---- 2.7 Post-cleaning class check -------------------------------------
message("\n---- 2.7 Post-cleaning class check (everything harmonized should now be numeric) ----")
print(sapply(hf_continuum_data, class))

## ---- 4. Before/after data-loss audit -----------------------------------
## Compares valid (non-missing) % for every variable, raw vs. cleaned.
## Two kinds of decrease are EXPECTED and CORRECT (not bugs):
##  - Yes/No-style variables (mcq160*, diq010, bpq*, smq*, kiq022, huq030,
##    dmdeduc2, pad200, paq180, ridexprg): "Don't know"/"Refused"/
##    "Borderline" text was non-missing as raw text, but correctly becomes
##    NA after harmonizing (those aren't real Yes/No/category answers).
##  - mcd180b-f: NHANES sentinel codes 99999/"Don't know", 77777/"Refused"
##    were non-missing as raw text/numbers, correctly become NA.
## Anything else showing a meaningful drop is worth investigating --
## exactly the kind of check that would have caught the ridageyr bug
## immediately instead of requiring a manual "why did this change" question.
valid_pct_after <- sapply(hf_continuum_data, function(x) round(100 * mean(!is.na(x)), 2))

common_vars <- intersect(names(valid_pct_before), names(valid_pct_after))
data_loss_report <- tibble(
  variable = common_vars,
  valid_pct_before = valid_pct_before[common_vars],
  valid_pct_after = valid_pct_after[common_vars],
  pct_point_change = round(valid_pct_after[common_vars] - valid_pct_before[common_vars], 2)
) %>%
  arrange(pct_point_change)

expected_recode_vars <- c("mcq160b", "mcq160c", "mcq160d", "mcq160e", "mcq160f",
                          "diq010", "bpq020", "bpq050a", "bpq080", "smq020", "smq040",
                          "kiq022", "huq030", "dmdeduc2", "pad200", "paq180", "ridexprg",
                          "hiq011", "riagendr", "fsdhh", "ssbnpl", "sspris", "ridreth1",
                          "mcd180b", "mcd180c", "mcd180d", "mcd180e", "mcd180f")
data_loss_report <- data_loss_report %>%
  mutate(expected_drop = variable %in% expected_recode_vars)

message("\n---- 4. Before/after valid-% audit (sorted by biggest drop first) ----")
print(data_loss_report, n = Inf)

unexpected_drops <- data_loss_report %>% filter(!expected_drop & pct_point_change < -0.01)
if (nrow(unexpected_drops) > 0) {
  message("\n*** UNEXPECTED data loss in variable(s) not on the known-recode list -- investigate: ***")
  print(unexpected_drops, n = Inf)
} else {
  message("\nNo unexpected data loss -- every drop is accounted for by",
          " Don't-know/Refused recoding or sentinel-code cleanup (99999/77777).")
}

## ---- 5. Save cleaned dataset + data-loss report -------------------------
dir.create("data", showWarnings = FALSE)
writexl::write_xlsx(hf_continuum_data, "data/HF_continuum_cleaned_dataset.xlsx")
message("\nSaved: data/HF_continuum_cleaned_dataset.xlsx")
message("Use THIS file (not HF_continuum_analytic_dataset.xlsx) as the input to hf_continuum_phenotype_analysis.R")

writexl::write_xlsx(data_loss_report, "data/HF_continuum_cleaning_data_loss_report.xlsx")
message("Saved: data/HF_continuum_cleaning_data_loss_report.xlsx (before/after valid-% audit)")