## ============================================================
## Biomarker-Informed Phenotyping of the HF Continuum
## Descriptive Variable Summary: 3-stage comparison
## ============================================================
##
## Produces THREE parallel per-variable summary tables:
##   Stage 1 RAW      -- data/HF_continuum_analytic_dataset.xlsx, exactly
##                        as extracted (messy: mixed text/numeric,
##                        inconsistent punctuation/casing across cycles)
##   Stage 2 RECODED  -- same raw data, but standardized into clean,
##                        consistent labels ("Every day," and "Every day"
##                        both become "Every day") -- Don't know/Refused/
##                        etc. are KEPT VISIBLE as their own category,
##                        not converted to NA yet
##   Stage 3 CLEANED  -- data/HF_continuum_cleaned_dataset.xlsx, the final
##                        analysis-ready data (Don't know/Refused -> NA)
##
## Output: data/HF_continuum_variable_summary_3stage.xlsx, one sheet per
## stage. For numeric variables: n, missing, mean, sd, median, min, max.
## For categorical variables: one row per category with n and %.

## ---- 0. Packages -------------------------------------------------
pkgs <- c("dplyr", "purrr", "readxl", "writexl", "tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(purrr)
library(readxl)
library(writexl)
library(tibble)

## ---- 1. Load both source files -----------------------------------------
raw_data <- readxl::read_excel("data/HF_continuum_analytic_dataset.xlsx")
names(raw_data) <- tolower(names(raw_data))   # same fix as the cleaning script

cleaned_data <- readxl::read_excel("data/HF_continuum_cleaned_dataset.xlsx")

message("Raw data: ", nrow(raw_data), " rows x ", ncol(raw_data), " cols")
message("Cleaned data: ", nrow(cleaned_data), " rows x ", ncol(cleaned_data), " cols")

## ---- 2. Generic summary function (works for numeric or character) -----
categorical_vars <- c(
  "riagendr", "ridreth1", "dmdeduc2", "ridexprg", "huq030", "huq010",
  "mcq160b", "mcq160c", "mcq160d", "mcq160e", "mcq160f",
  "diq010", "bpq020", "bpq050a", "bpq080", "smq020", "smq040",
  "pad200", "paq180", "hiq011", "kiq022", "fsdhh",
  "ssbnpl", "sspris", "eligstat", "mortstat", "ucod_leading",
  "diabetes", "hyperten"
)

summarize_column <- function(x, varname, categorical_vars = character(0)) {
  is_categorical <- is.character(x) || varname %in% categorical_vars
  if (!is_categorical) {
    tibble(
      variable = varname, type = "numeric",
      category = NA_character_,
      n = sum(!is.na(x)), pct_of_total = round(100 * sum(!is.na(x)) / length(x), 1),
      mean = round(mean(x, na.rm = TRUE), 2), sd = round(sd(x, na.rm = TRUE), 2),
      median = round(median(x, na.rm = TRUE), 2),
      min = round(suppressWarnings(min(x, na.rm = TRUE)), 2),
      max = round(suppressWarnings(max(x, na.rm = TRUE)), 2)
    )
  } else {
    tab <- table(x, useNA = "ifany")
    tibble(
      variable = varname, type = "categorical",
      category = ifelse(is.na(names(tab)), "(missing)", names(tab)),
      n = as.integer(tab),
      pct_of_total = round(100 * as.integer(tab) / length(x), 1),
      mean = NA_real_, sd = NA_real_, median = NA_real_, min = NA_real_, max = NA_real_
    )
  }
}

summarize_dataset <- function(df, categorical_vars = character(0)) {
  purrr::imap_dfr(df, ~ summarize_column(.x, .y, categorical_vars))
}

## ---- 3. Stage 1: RAW summary --------------------------------------------
message("\nBuilding Stage 1 (raw) summary...")
stage1_summary <- summarize_dataset(raw_data, categorical_vars)

## ---- 4. Stage 2: RECODED (standardized labels, non-response visible) ---
recode_label_yesno <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|^yes$", x_chr)] <- "Yes"
  out[grepl("^2$|^no$", x_chr)] <- "No"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[grepl("refused", x_chr)] <- "Refused"
  out[grepl("borderline", x_chr)] <- "Borderline"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_smq040 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("every day", x_chr)] <- "Every day"
  out[grepl("some days", x_chr)] <- "Some days"
  out[grepl("not at all", x_chr)] <- "Not at all"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[grepl("refused", x_chr)] <- "Refused"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_riagendr <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|^male$", x_chr)] <- "Male"
  out[grepl("^2$|^female$", x_chr)] <- "Female"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_ridexprg <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|pregnant", x_chr)] <- "Pregnant"
  out[grepl("^2$|not pregnant", x_chr)] <- "Not pregnant"
  out[grepl("^3$|cannot ascertain", x_chr)] <- "Cannot ascertain"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_dmdeduc2 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|9th grade|less than 9", x_chr)] <- "Less than 9th grade"
  out[grepl("^2$|9-11th|11th grade", x_chr)] <- "9-11th grade"
  out[grepl("^3$|high school", x_chr)] <- "High school grad/GED"
  out[grepl("^4$|some college|aa degree", x_chr)] <- "Some college/AA degree"
  out[grepl("^5$|college graduate", x_chr)] <- "College graduate+"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[grepl("refused", x_chr)] <- "Refused"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_huq030 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|^yes$", x_chr)] <- "Yes"
  out[grepl("^2$|no place", x_chr)] <- "No place"
  out[grepl("^3$|more than one", x_chr)] <- "More than one place"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_huq010 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|excellent", x_chr)] <- "Excellent"
  out[grepl("^3$|good", x_chr)] <- "Good"
  out[grepl("^2$|very good", x_chr)] <- "Very good"
  out[grepl("^4$|fair", x_chr)] <- "Fair"
  out[grepl("^5$|poor", x_chr)] <- "Poor"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[grepl("refused", x_chr)] <- "Refused"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_pad200 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|^yes$", x_chr)] <- "Yes"
  out[grepl("^2$|^no$", x_chr)] <- "No"
  out[grepl("unable", x_chr)] <- "Unable to do activity"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[grepl("refused", x_chr)] <- "Refused"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_paq180 <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^1$|not walk about", x_chr)] <- "1: Sit, don't walk much"
  out[grepl("^2$|stand or walk", x_chr)] <- "2: Stand/walk a lot"
  out[grepl("^3$|light load|climb stairs", x_chr)] <- "3: Lift light loads/climb"
  out[grepl("^4$|heavy work|heavy loads", x_chr)] <- "4: Heavy work/carry heavy loads"
  out[grepl("don.t know", x_chr)] <- "Don't know"
  out[grepl("refused", x_chr)] <- "Refused"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_fsdhh <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  x_chr <- tolower(trimws(as.character(x)))
  case_when(
    !is.na(x_num) & x_num == 1 ~ "Full food security",
    !is.na(x_num) & x_num == 2 ~ "Marginal food security",
    !is.na(x_num) & x_num == 3 ~ "Low food security",
    !is.na(x_num) & x_num == 4 ~ "Very low food security",
    grepl("full", x_chr) ~ "Full food security",
    grepl("marginal", x_chr) ~ "Marginal food security",
    grepl("low", x_chr) & !grepl("very", x_chr) ~ "Low food security",
    grepl("very low", x_chr) ~ "Very low food security",
    !is.na(x_chr) ~ paste0("Other/unrecognized: ", x_chr),
    TRUE ~ NA_character_
  )
}

recode_label_ssbnpl <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("within the detection", x_chr)] <- "Within detection limits"
  out[grepl("below lower", x_chr)] <- "Below lower detection limit"
  out[grepl("above upper", x_chr)] <- "Above upper detection limit"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_sspris <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x_chr))
  out[grepl("^pristine$", x_chr)] <- "Pristine"
  out[grepl("non-pristine", x_chr)] <- "Non-pristine"
  out[is.na(out) & !is.na(x_chr)] <- paste0("Other/unrecognized: ", x_chr)
  out
}

recode_label_age_sentinel <- function(x) {
  x_chr <- trimws(as.character(x))
  x_num <- suppressWarnings(as.numeric(x_chr))
  case_when(
    !is.na(x_num) & x_num == 99999 ~ "Don't know",
    !is.na(x_num) & x_num == 77777 ~ "Refused",
    !is.na(x_num) ~ as.character(x_num),   # real age, kept as its own value
    !is.na(x_chr) & x_chr != "" ~ paste0("Other/unrecognized: ", x_chr),
    TRUE ~ NA_character_
  )
}

recoders <- list(
  hiq011 = recode_label_yesno, kiq022 = recode_label_yesno,
  mcq160b = recode_label_yesno, mcq160c = recode_label_yesno,
  mcq160d = recode_label_yesno, mcq160e = recode_label_yesno, mcq160f = recode_label_yesno,
  diq010 = recode_label_yesno, bpq020 = recode_label_yesno, bpq050a = recode_label_yesno,
  bpq080 = recode_label_yesno, smq020 = recode_label_yesno, pad200 = recode_label_pad200,
  smq040 = recode_label_smq040, riagendr = recode_label_riagendr,
  ridexprg = recode_label_ridexprg, dmdeduc2 = recode_label_dmdeduc2,
  huq030 = recode_label_huq030, huq010 = recode_label_huq010, paq180 = recode_label_paq180, fsdhh = recode_label_fsdhh,
  ssbnpl = recode_label_ssbnpl, sspris = recode_label_sspris,
  mcd180b = recode_label_age_sentinel, mcd180c = recode_label_age_sentinel,
  mcd180d = recode_label_age_sentinel, mcd180e = recode_label_age_sentinel,
  mcd180f = recode_label_age_sentinel
)

message("Building Stage 2 (recoded, non-response visible) summary...")
recoded_data <- raw_data
for (v in names(recoders)) {
  if (v %in% names(recoded_data)) recoded_data[[v]] <- recoders[[v]](recoded_data[[v]])
}

bin_age_labels <- function(x) {
  is_sentinel <- x %in% c("Don't know", "Refused") | grepl("^Other/unrecognized", x)
  x_num <- suppressWarnings(as.numeric(x))
  out <- x
  bucket <- cut(x_num, breaks = c(-1, 19, 39, 59, 79, 200),
                labels = c("0-19", "20-39", "40-59", "60-79", "80+"))
  out[!is_sentinel & !is.na(x_num)] <- as.character(bucket[!is_sentinel & !is.na(x_num)])
  out
}
for (v in c("mcd180b", "mcd180c", "mcd180d", "mcd180e", "mcd180f")) {
  if (v %in% names(recoded_data)) recoded_data[[v]] <- bin_age_labels(recoded_data[[v]])
}

stage2_summary <- summarize_dataset(recoded_data, categorical_vars)

## ---- 5. Stage 3: CLEANED summary ----------------------------------------
message("Building Stage 3 (cleaned) summary...")
stage3_summary <- summarize_dataset(cleaned_data, categorical_vars)

## ---- 6. Save all three as one workbook ----------------------------------
dir.create("data", showWarnings = FALSE)
writexl::write_xlsx(
  list(
    "Stage1_Raw" = stage1_summary,
    "Stage2_Recoded_with_NR" = stage2_summary,
    "Stage3_Cleaned" = stage3_summary
  ),
  "data/HF_continuum_variable_summary_3stage.xlsx"
)
message("\nSaved: data/HF_continuum_variable_summary_3stage.xlsx (3 sheets: Stage1_Raw, Stage2_Recoded_with_NR, Stage3_Cleaned)")