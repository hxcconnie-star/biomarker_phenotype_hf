## ============================================================
## Append-only snippet: labeled version of the Stage 3 cleaned dataset
## ============================================================
## NOTE: this is a SECONDARY/reporting output only -- for the next step
## (phenotype analysis), use the NUMERIC-CODED HF_continuum_cleaned_dataset.xlsx,
## NOT this labeled file. Text labels ("Yes"/"No") would break every ==
## comparison in the analysis script (e.g. mcq160b == 1), the exact bug
## class this whole project has been fixing -- this file exists purely
## for human-readable review/sharing, never as analysis input.
##
## Always freshly loads the file (no exists() check) to avoid using a
## stale/wrong object left over from a different R session state.
library(readxl)
library(writexl)
library(dplyr)

cleaned_data <- readxl::read_excel("data/HF_continuum_cleaned_dataset.xlsx")
message("Loaded fresh: ", nrow(cleaned_data), " rows x ", ncol(cleaned_data), " cols")

## Label dictionaries, matching the coding used throughout this project's
## harmonizers/derivation scripts. ridreth1 added now that the cleaning
## script numerically recodes it too (1=Mexican American, 2=Other
## Hispanic, 3=Non-Hispanic White, 4=Non-Hispanic Black, 5=Other Race).
label_maps <- list(
  riagendr = c(`1` = "Male", `2` = "Female"),
  ridreth1 = c(`1` = "Mexican American", `2` = "Other Hispanic",
               `3` = "Non-Hispanic White", `4` = "Non-Hispanic Black",
               `5` = "Other Race - Including Multi-Racial"),
  dmdeduc2 = c(`1` = "Less than 9th grade", `2` = "9-11th grade",
               `3` = "High school grad/GED", `4` = "Some college/AA degree",
               `5` = "College graduate+"),
  ridexprg = c(`1` = "Pregnant", `2` = "Not pregnant", `3` = "Cannot ascertain"),
  huq030   = c(`1` = "Yes", `2` = "No place", `3` = "More than one place"),
  mcq160b = c(`1` = "Yes", `2` = "No"), mcq160c = c(`1` = "Yes", `2` = "No"),
  mcq160d = c(`1` = "Yes", `2` = "No"), mcq160e = c(`1` = "Yes", `2` = "No"),
  mcq160f = c(`1` = "Yes", `2` = "No"),
  diq010   = c(`1` = "Yes", `2` = "No"),
  bpq020   = c(`1` = "Yes", `2` = "No"), bpq050a = c(`1` = "Yes", `2` = "No"),
  bpq080   = c(`1` = "Yes", `2` = "No"),
  smq020   = c(`1` = "Yes", `2` = "No"),
  smq040   = c(`1` = "Every day", `2` = "Some days", `3` = "Not at all"),
  pad200   = c(`1` = "Yes", `2` = "No"),
  paq180   = c(`1` = "Sit, don't walk much", `2` = "Stand/walk a lot",
               `3` = "Lift light loads/climb", `4` = "Heavy work/carry heavy loads"),
  hiq011   = c(`1` = "Yes", `2` = "No"),
  kiq022   = c(`1` = "Yes", `2` = "No"),
  fsdhh    = c(`1` = "Full food security", `2` = "Marginal food security",
               `3` = "Low food security", `4` = "Very low food security"),
  ssbnpl   = c(`0` = "Within detection limits", `1` = "Below lower detection limit",
               `2` = "Above upper detection limit"),
  sspris   = c(`0` = "Non-pristine", `1` = "Pristine"),
  mortstat = c(`0` = "Assumed alive", `1` = "Assumed deceased"),
  diabetes = c(`0` = "No", `1` = "Yes"),
  hyperten = c(`0` = "No", `1` = "Yes")
  ## eligstat and ucod_leading deliberately NOT mapped -- their exact
  ## coding wasn't independently verified against NCHS documentation in
  ## this project (unlike everything else here). Left as numeric codes
  ## rather than risk a wrong label.
)

labeled_data <- cleaned_data
for (v in names(label_maps)) {
  if (!v %in% names(labeled_data)) next
  n_before_nonmissing <- sum(!is.na(labeled_data[[v]]))
  labeled_data[[v]] <- unname(label_maps[[v]][as.character(labeled_data[[v]])])
  n_after_nonmissing <- sum(!is.na(labeled_data[[v]]))
  ## Verification: if labeling worked, non-missing count should be
  ## IDENTICAL before/after (every real code should have found a label).
  ## If it prints a mismatch here, something in label_maps doesn't match
  ## the actual codes in the data -- check the printed table below it.
  if (n_before_nonmissing != n_after_nonmissing) {
    message("  WARNING: ", v, " -- non-missing count changed from ",
            n_before_nonmissing, " to ", n_after_nonmissing,
            " after labeling. Some code(s) didn't match label_maps. Original values:")
    print(table(cleaned_data[[v]], useNA = "always"))
  }
}

message("\n---- Verification: sample labeled variable (mcq160b) ----")
print(table(labeled_data$mcq160b, useNA = "always"))

dir.create("data", showWarnings = FALSE)
writexl::write_xlsx(labeled_data, "data/HF_continuum_cleaned_dataset_labeled.xlsx")
message("Saved: data/HF_continuum_cleaned_dataset_labeled.xlsx",
        " (same as HF_continuum_cleaned_dataset.xlsx, but with readable labels",
        " instead of numeric codes for categorical variables)")
