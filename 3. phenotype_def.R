## ============================================================
## Biomarker-Informed Phenotyping of the HF Continuum
## Analysis Pipeline, Part 2: from hf_continuum_data through
## Results §4 (Trajectory analysis) and §7 (Sensitivity)
## ============================================================
##
## Covers:
##   Methods §3  Mapping the 2026 HF definition to NHANES variables
##   Methods §4  HF-continuum phenotype definition
##   Methods §5  Etiologic and social vulnerability characterization
##   Methods §6  Longitudinal trajectory analysis
##   Methods §8.2 Missing data (multiple imputation)
##   Results §1  Analytic cohort construction
##   Results §2  Distribution of HF-continuum phenotypes
##   Results §3  Baseline characteristics across phenotypes
##   Results §4  Mortality trajectories across phenotypes
##   Results §7  Sensitivity (complete-case vs MI) -- partial: only the
##               complete-case-vs-MI Cox comparison. Alternate NT-proBNP
##               thresholds and replication in the full 1999-2018 cohort
##               are NOT covered here.
##
## SCRIPT SECTION MAP (numbered sections in this file, in order):
##   0-2   Setup, load data, derived variables/phenotype construction
##   3-4   Results §1: analytic cohort construction + survey design
##   5     Results §2: distribution of HF-continuum phenotypes
##   6-7   Results §3: baseline characteristics (+ save)
##   ----- MULTIPLE IMPUTATION INSERTION POINT (Methods §8.2) -----
##   8     Multiple imputation (covariates only)
##   9.1-9.4  Results §4: mortality trajectories (KM, adjusted Cox, RMST)
##   10    Results §7: sensitivity (complete-case vs MI)
##   11    Save Results §4 and §7 outputs
##
## PIPELINE ORDER: nhanes_hf_continuum_extraction.R -> data cleaning
## (hf_continuum_data_cleaning.R) -> THIS SCRIPT. Data cleaning and
## validation (class audits, harmonizing text categories to numeric
## codes, range/plausibility checks, logical consistency checks) live in
## their own script now -- this one assumes hf_continuum_data is already
## clean and just builds derived variables/phenotypes/results on top.
##
## KEY ASSUMPTION (flag if this isn't what you intended): the four HF-
## continuum levels are treated as MUTUALLY EXCLUSIVE via a strict
## priority hierarchy: Stage C (self-reported HF) > Stage B (elevated
## NT-proBNP) > Stage A (>=1 risk factor) > Level 0 (no apparent risk).
## Someone with self-reported HF is Stage C regardless of biomarker or
## risk-factor status.
##
## Numeric thresholds below (BP, BMI, eGFR, age-risk cutoff, etc.) are
## reasonable defaults where your protocol specified a category but not
## an exact cutpoint. They're called out inline -- adjust as needed.
##
## SECOND KEY ASSUMPTION (flag if wrong): the Cox model adjustment set
## (age, sex, race/ethnicity, education, poverty-income ratio) is my
## choice, not specified by your protocol -- it's a standard baseline
## demographic/SES adjustment for a mortality model, but you may want
## more (e.g. individual comorbidities) or fewer covariates.
##
## WHERE MULTIPLE IMPUTATION FITS (per your question): MI is NOT used to
## assign the HF-continuum phenotype -- Section 3's inclusion criterion
## #5 ("available core variables needed to assign phenotype") means
## anyone missing those is excluded from the cohort entirely, not
## imputed. MI (Section 9 below) runs AFTER cohort construction and
## phenotyping (Results §1-3), and only imputes the ADJUSTMENT COVARIATES
## needed for the Cox models in Section 6/Results §4 -- never mortstat,
## permth_exm, or hf_continuum_level itself (all three enter the
## imputation model as complete predictors, per protocol §8.2, but are
## never themselves imputed).

## ---- 0. Packages -------------------------------------------------
pkgs <- c("dplyr", "purrr", "readxl", "writexl", "survey", "tableone", "stringr",
          "mice", "mitools", "survival", "ggplot2", "scales")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(purrr)
library(readxl)
library(writexl)
library(survey)
library(tableone)
library(stringr)
library(mice)
library(mitools)
library(survival)
library(ggplot2)
library(scales)
library(patchwork)

## ---- 1. Load CLEANED data ---------------------------------------------
## Reads the output of hf_continuum_data_cleaning.R, NOT the raw
## extraction output -- run that script first if this file doesn't exist.
hf_continuum_data <- readxl::read_excel("./data/HF_continuum_cleaned_dataset.xlsx")
message("Loaded cleaned hf_continuum_data: ", nrow(hf_continuum_data), " rows x ",
        ncol(hf_continuum_data), " cols")

## ---- 2. Derived variables (Methods §3 mapping / §4 phenotype / §5 etiologic domains) ----

## 3.1 eGFR: CKD-EPI 2021 race-free creatinine equation (Inker et al 2021)
compute_egfr <- function(scr, age, sex_male) {
  kappa <- ifelse(sex_male, 0.9, 0.7)
  alpha <- ifelse(sex_male, -0.302, -0.241)
  sex_mult <- ifelse(sex_male, 1, 1.012)
  142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^-1.200 * 0.9938^age * sex_mult
}

## 3.2 UACR + 3.3 CKD category (simplified KDIGO risk tiers from eGFR x UACR)
kdigo_risk_category <- function(egfr, uacr) {
  g_cat <- case_when(
    egfr >= 90 ~ "G1", egfr >= 60 ~ "G2", egfr >= 45 ~ "G3a",
    egfr >= 30 ~ "G3b", egfr >= 15 ~ "G4", TRUE ~ "G5"
  )
  a_cat <- case_when(uacr < 30 ~ "A1", uacr < 300 ~ "A2", TRUE ~ "A3")
  ## Simplified KDIGO heat-map -> 4-tier risk category
  case_when(
    g_cat %in% c("G1", "G2") & a_cat == "A1" ~ "Low risk",
    g_cat %in% c("G1", "G2") & a_cat == "A2" ~ "Moderately increased risk",
    g_cat %in% c("G1", "G2") & a_cat == "A3" ~ "High risk",
    g_cat == "G3a" & a_cat == "A1" ~ "Moderately increased risk",
    g_cat == "G3a" & a_cat %in% c("A2", "A3") ~ "High risk",
    g_cat == "G3b" & a_cat == "A1" ~ "High risk",
    g_cat %in% c("G3b", "G4", "G5") ~ "Very high risk",
    TRUE ~ NA_character_
  )
}

hf_continuum_data <- hf_continuum_data %>%
  mutate(
    male = riagendr == 1,
    egfr = compute_egfr(lbxscr, ridageyr, male),
    uacr = (urxuma / urxucr) * 100,
    ckd_category = kdigo_risk_category(egfr, uacr)
  )

## 3.4 NT-proBNP transforms: log, age/sex-specific percentile among
## participants WITHOUT self-reported HF, and the elevated indicator.
## Age-sex strata: 10-year age bands x sex (adjust bin width if your
## sample size in a stratum gets too small for a stable 90th percentile).
hf_continuum_data <- hf_continuum_data %>%
  mutate(
    nt_probnp_log = log(ssbnp),
    age_band = cut(ridageyr, breaks = c(19, 29, 39, 49, 59, 69, 79, Inf),
                   labels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))
  ) %>%
  group_by(age_band, riagendr) %>%
  mutate(
    ## Reference distribution = same age-sex stratum, no self-reported HF.
    ## IMPORTANT: mcq160b == 2 evaluates to NA (not FALSE) when mcq160b is
    ## itself NA, and R keeps NA positions when subsetting with a logical
    ## index containing NA -- so without !is.na(mcq160b) first, missing-HF-
    ## status rows leak NA values into ref, and if a stratum is small/odd
    ## enough that ref ends up ALL NA, ecdf() errors outright even though
    ## length(ref) looked >= 10. Also guard the NA age_band group itself
    ## (people with missing/out-of-range age) -- not a meaningful stratum.
    nt_probnp_pctile = {
      if (is.na(age_band[1])) {
        rep(NA_real_, n())
      } else {
        ref <- ssbnp[!is.na(mcq160b) & mcq160b == 2 & !is.na(ssbnp)]
        if (length(ref) < 10) {
          rep(NA_real_, n())   # stratum too small for a stable percentile
        } else {
          ecdf(ref)(ssbnp) * 100
        }
      }
    }
  ) %>%
  ungroup() %>%
  mutate(nt_probnp_elevated = as.numeric(nt_probnp_pctile > 90))

## 3.5 Disease duration (HF)
hf_continuum_data <- hf_continuum_data %>%
  mutate(disease_duration_hf = ridageyr - mcd180b)

## 3.5b Cardiovascular/heart-disease mortality (Protocol §7, Outcome #4)
## Was PLANNED in the codebook's "Derived Variables" sheet early on but
## never actually implemented -- caught during a full variable-list
## review against Section 7. UCOD_LEADING coding (NCHS's standard 10-
## category leading-cause-of-death recode, confirmed via two independent
## sources, consistent with the 1/2/3/6/10 values seen in this project's
## own data profiling): 1 = Diseases of heart, 5 = Cerebrovascular
## diseases. mortality_cvd is only meaningful for mortstat==1 (deceased);
## it's NA for people who are alive or not eligible for mortality linkage.
hf_continuum_data <- hf_continuum_data %>%
  mutate(mortality_cvd = case_when(
    is.na(mortstat) | mortstat == 0 ~ NA_real_,
    ucod_leading %in% c(1, 5) ~ 1,
    !is.na(ucod_leading) ~ 0,
    TRUE ~ NA_real_
  ))

message("\n---- mortality_cvd check (should only be non-NA for decedents) ----")
print(table(hf_continuum_data$mortality_cvd, hf_continuum_data$mortstat, useNA = "always"))

## 3.6 Dyslipidemia composite
hf_continuum_data <- hf_continuum_data %>%
  mutate(
    dyslipidemia = as.numeric(
      bpq080 == 1 |
        lbxtc >= 200 |
        lbdldl >= 160 |
        (male & lbdhdd < 40) | (!male & lbdhdd < 50) |
        lbxtr >= 150
    )
  )

## 3.7 Stage A risk-factor components + composite flag
## Per your explicit instruction: use direct self-reported/simple-threshold
## variables rather than derived clinical measures. mean_sbp/mean_dbp are
## still computed below (needed for Table 1 / trend tests -- protocol
## Section 7 lists Systolic/Diastolic BP as their own reported variables)
## but are NO LONGER used to trigger the hypertension flag.
##
## TWO ASSUMPTIONS I made that need your confirmation:
##  1. flag_diabetes now uses DIQ010 alone (dropped LBXGH/HbA1c) -- you
##     didn't explicitly say to drop HbA1c, but it matches the same
##     self-report-only pattern as hypertension/CKD. If you want HbA1c to
##     still count (either combined with DIQ010, or as its own separate
##     Stage A trigger item), tell me and I'll add it back.
##  2. flag_albuminuria uses URXUMA (urine albumin, ug/mL) >= 30 as a raw
##     concentration threshold. Note URXUMA (ug/mL) and URXUMASI (mg/L)
##     are numerically identical (1 ug/mL = 1 mg/L), so only one is
##     needed -- urxuma already covers both. The 30 threshold mirrors the
##     conventional microalbuminuria value of 30, but that value is
##     normally applied to an ALBUMIN-CREATININE RATIO (mg/g), not a raw
##     concentration (mg/L) -- using it on raw URXUMA is a simplification
##     and not a standard clinical cutoff. Confirm this is what you want,
##     or tell me the threshold you'd prefer.
##
## flag_ckd and flag_albuminuria are kept as SEPARATE flags (matching your
## list, which has them as separate Stage A items #4 and #5) -- unlike the
## etiologic "Kidney-related" domain in 3.9, which per protocol Section 5
## combines "reduced eGFR or albuminuria" into one domain and still uses
## the eGFR-based definition there.
hf_continuum_data <- hf_continuum_data %>%
  rowwise() %>%
  mutate(
    mean_sbp = mean(c(bpxsy1, bpxsy2, bpxsy3, bpxsy4), na.rm = TRUE),
    mean_dbp = mean(c(bpxdi1, bpxdi2, bpxdi3, bpxdi4), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    mean_sbp = ifelse(is.nan(mean_sbp), NA_real_, mean_sbp),
    mean_dbp = ifelse(is.nan(mean_dbp), NA_real_, mean_dbp),
    
    flag_hypertension = as.numeric(bpq020 == 1),
    flag_diabetes      = as.numeric(diq010 == 1),
    flag_obesity       = as.numeric(bmxbmi >= 30),
    flag_ckd           = as.numeric(kiq022 == 1),
    flag_albuminuria   = as.numeric(urxuma >= 30),
    flag_chd_mi        = as.numeric(mcq160c == 1 | mcq160e == 1),
    flag_stroke        = as.numeric(mcq160f == 1),
    flag_smoking       = as.numeric(smq040 %in% c(1, 2)),
    flag_older_age     = as.numeric(ridageyr >= 65),
    flag_social_vuln   = as.numeric(
      hiq011 == 2 | dmdeduc2 <= 2 | indfmpir < 1 | huq030 == 2
    ),
    
    stage_a_risk_factor = as.numeric(
      flag_hypertension == 1 | flag_diabetes == 1 | flag_obesity == 1 |
        flag_ckd == 1 | flag_albuminuria == 1 | flag_chd_mi == 1 | flag_stroke == 1 |
        flag_smoking == 1 | flag_older_age == 1 | flag_social_vuln == 1
    )
  )

## 3.8 HF-continuum phenotype (0-3) -- see priority-order note at top
## CONFIRMED BUG #1, FIXED: the original logic only treated a person as
## "indeterminate" (NA) when BOTH nt_probnp_elevated AND
## stage_a_risk_factor were missing. A person with NO self-reported HF,
## a CONFIRMED absence of Stage A risk factors, and MISSING NT-proBNP
## fell through to `TRUE ~ 0L` -- silently assigned "No apparent HF
## risk" even though their true Stage B status was never known. This is
## the failure mode flagged in review: missing NT-proBNP treated as if
## it were a confirmed non-elevated value.
##
## CONFIRMED BUG #2, FIXED (caught on further review): the first fix
## checked `stage_a_risk_factor == 1` BEFORE resolving NT-proBNP
## missingness. This let someone with missing NT-proBNP but a confirmed
## OTHER risk factor (e.g. diagnosed diabetes) be classified as Stage A
## -- but the classification hierarchy is explicitly priority-ordered
## C > B > A > 0, meaning a person can only be assigned to a LOWER-
## priority level once every HIGHER-priority level has been positively
## RULED OUT, not merely left unevaluated. Checking Stage A before
## Stage B's missingness was resolved effectively let missing data
## silently reverse the stated priority order for that person -- someone
## whose true phenotype might have been Stage B (had NT-proBNP been
## measured and elevated) could be assigned Stage A instead, simply
## because Stage A happened to be evaluable and Stage B was not.
##
## Fix: NT-proBNP missingness is now resolved immediately after the
## Stage B check and BEFORE Stage A is ever evaluated. A missing
## NT-proBNP value (in anyone not already resolved as Stage C or
## confirmed-elevated Stage B) is always NA, regardless of Stage A risk
## factor status. Stage A is only evaluated once NT-proBNP is CONFIRMED
## non-elevated, and Level 0 is only reached once NT-proBNP is confirmed
## non-elevated AND stage_a_risk_factor is confirmed 0.

hf_continuum_data <- hf_continuum_data %>%
  mutate(
    hf_continuum_level = case_when(
      mcq160b == 1               ~ 3L,             # Stage C (always directly evaluable via self-report)
      is.na(mcq160b)             ~ NA_integer_,
      nt_probnp_elevated == 1    ~ 2L,             # Stage B / pre-HF (confirmed elevated)
      is.na(nt_probnp_elevated)  ~ NA_integer_,    # Stage B cannot be ruled out -- indeterminate regardless of Stage A status
      stage_a_risk_factor == 1   ~ 1L,             # Stage A (only reached once Stage B is confirmed ruled out)
      is.na(stage_a_risk_factor) ~ NA_integer_,    # Stage A itself indeterminate
      TRUE                       ~ 0L               # No apparent HF risk (Stage B and Stage A both confirmed ruled out)
    ),
    hf_continuum_label = factor(
      hf_continuum_level,
      levels = 0:3,
      labels = c("No apparent HF risk", "Stage A (at risk)",
                 "Stage B (pre-HF, biomarker)", "Stage C (self-reported HF)")
    )
  )

## 3.9 Etiologic domain flags (Methods §5) -- secondary characterization,
## kept subordinate to the main 4-level phenotype above.
## NOTE: these are BROADER than the Stage A risk-factor flags in 3.7,
## matching your plan's Section 5 domain table exactly rather than
## reusing Section 4's (now self-report-only) Stage A risk-factor flags:
##   - etio_hypertensive uses self-report OR measured BP OR medication use
##     (Section 5: "Self-reported hypertension, measured elevated BP, or
##     antihypertensive use"), unlike flag_hypertension (3.7), which is
##     now self-report-only (BPQ020) per your instruction.
##   - etio_ischemic adds angina (mcq160d), which Section 5 explicitly
##     includes ("CHD, angina, or MI") but Section 4's Stage A list does
##     not (only "prior CHD or MI").
##   - etio_kidney uses reduced eGFR OR albuminuria (Section 5: "Reduced
##     eGFR or albuminuria"), unlike flag_ckd (3.7), which is now self-
##     report-only (KIQ022) per your instruction.
##   - etio_social_vuln adds food insecurity (fsdhh = Low or Very Low,
##     i.e. 3/4), which Section 5 explicitly lists but Section 4's Stage A
##     list does not.
## Deliberately NOT changing the Stage A flags themselves (3.7) or
## stage_a_risk_factor -- that would silently widen the Stage A/HF-
## continuum-level definition itself, which Section 4 does not call for.
hf_continuum_data <- hf_continuum_data %>%
  mutate(
    etio_hypertensive = as.numeric(
      bpq020 == 1 | bpq050a == 1 | mean_sbp >= 130 | mean_dbp >= 80
    ),
    etio_ischemic      = as.numeric(mcq160c == 1 | mcq160d == 1 | mcq160e == 1),
    etio_metabolic     = as.numeric(flag_diabetes == 1 | flag_obesity == 1 | dyslipidemia == 1),
    etio_kidney        = as.numeric(egfr < 60 | flag_albuminuria == 1),
    etio_vascular      = flag_stroke,
    etio_social_vuln   = as.numeric(flag_social_vuln == 1 | fsdhh %in% c(3, 4))
  )

## ============================================================
## RESULTS §1: Analytic cohort construction
## ============================================================
## ---- 3. Analytic cohort construction (Results §1) -------------------
## Applies protocol Section 3 inclusion/exclusion criteria in order,
## printing the running N at each step -- use this to build your CONSORT
## flow diagram.
n_start <- nrow(hf_continuum_data)
message("\n---- Analytic cohort construction ----")
message("Starting N (all loaded records): ", n_start)

cohort <- hf_continuum_data

cohort <- cohort %>% filter(ridageyr >= 20)
message("After age >= 20: ", nrow(cohort))

cohort <- cohort %>% filter(!is.na(mcq160b))
message("After completed medical conditions questionnaire (non-missing mcq160b): ", nrow(cohort))

## CONFIRMED BUG, FIXED: the 6-year pooled weight was previously built as
## WTMEC2YR/3 for every participant -- i.e. from the GENERAL MEC exam
## weight, uniformly divided by 3. This is the standard, documented way
## to pool three ordinary 2-year NHANES cycles, but it is NOT correct
## for this analysis: NT-proBNP was measured only in the NHANES
## surplus-sera biomarker sub-study, which has its OWN subsample-specific
## weights (WTSSCB2Y, WTSSCB4Y) that CDC explicitly documents must be
## used instead of WTMEC2YR for any analysis of this biomarker.
##
## Correct construction (per NHANES documentation for SSBNP_A/SSCARD_A/
## SSTROP_A, confirmed in review): WTSSCB4Y already pools the 1999-2000
## and 2001-2002 cycles into one 4-year subsample weight; WTSSCB2Y is the
## standalone 2003-2004 2-year subsample weight. Each participant has
## exactly one of the two populated, never both. To combine these into a
## single 6-year pooled weight, each is scaled by the fraction of the
## pooled 6-year window its own cycle represents: the 4-year group
## contributes 4/6 = 2/3, the 2-year group contributes 2/6 = 1/3 -- NOT a
## uniform division by 3 of a single weight variable.
cohort <- cohort %>%
  mutate(
    wtmec6yr = case_when(
      !is.na(wtsscb4y) ~ wtsscb4y * (2 / 3),
      !is.na(wtsscb2y) ~ wtsscb2y * (1 / 3),
      TRUE ~ NA_real_
    )
  )
message("wtmec6yr constructed from biomarker subsample weights: ",
        sum(!is.na(cohort$wtmec6yr)), " of ", nrow(cohort), " have a valid weight",
        " (", sum(!is.na(cohort$wtsscb4y)), " from the 1999-2002 4-year group, ",
        sum(!is.na(cohort$wtsscb2y)), " from the 2003-2004 2-year group)")

## CONFIRMED via review: because WTSSCB2Y/4Y are only defined for
## participants actually selected into the surplus-sera biomarker
## sub-study, switching to this weight means the entire analytic cohort
## -- not just participants who need NT-proBNP for phenotype assignment
## -- must be restricted to biomarker sub-study members. This resolves
## an issue flagged in review: a participant with self-reported HF
## (Stage C) does not need NT-proBNP to be classified, but if they were
## never selected into the biomarker sub-study, they have no valid
## WTSSCB weight and cannot be validly included in any survey-weighted
## analysis using this weight. The filter below excludes such
## participants explicitly and reports how many are lost this way,
## rather than letting them fail silently later with a missing weight.
n_before_biomarker_subsample_filter <- nrow(cohort)
n_stageC_without_biomarker_weight <- cohort %>%
  filter(is.na(wtmec6yr) & mcq160b == 1) %>%
  nrow()
message("Of participants not in the biomarker sub-study (missing wtmec6yr): ",
        n_stageC_without_biomarker_weight,
        " have self-reported HF (Stage C) and would otherwise not have needed NT-proBNP",
        " for phenotype assignment, but are excluded here because they lack a valid",
        " survey weight for this analysis.")

cohort <- cohort %>% filter(!is.na(sdmvpsu) & !is.na(sdmvstra) & !is.na(wtmec6yr))
message("After available exam data + survey design variables + biomarker sub-study",
        " membership (non-missing wtmec6yr): ", nrow(cohort),
        " (excluded ", n_before_biomarker_subsample_filter - nrow(cohort), ")")

cohort <- cohort %>% filter(!is.na(eligstat))
message("After non-missing mortality-linkage eligibility: ", nrow(cohort))
cohort <- cohort %>% filter(eligstat == 1)
message("After eligible for public-use mortality linkage (eligstat==1): ", nrow(cohort))

## Diagnostic breakdown BEFORE filtering -- separately counts how many
## of the pending exclusions are attributable specifically to missing
## NT-proBNP (biomarker availability) versus other reasons, addressing
## the request to report biomarker availability as its own step rather
## than folding it into a single generic "core variables" exclusion.
## Given the priority-order fix above (NT-proBNP missingness is now
## resolved before Stage A is ever evaluated), n_missing_ntprobnp_only is
## expected to capture the large majority of indeterminate cases here --
## this is the intended, priority-order-consistent consequence discussed
## in the comment above, not a sign that something is newly broken.
n_before_core_filter <- nrow(cohort)
n_missing_ntprobnp_only <- cohort %>%
  filter(is.na(hf_continuum_level) & is.na(ssbnp) & !is.na(mcq160b)) %>%
  nrow()
n_missing_other_core <- cohort %>%
  filter(is.na(hf_continuum_level) & !is.na(ssbnp)) %>%
  nrow()
message("Of participants with indeterminate HF-continuum phenotype: ",
        n_missing_ntprobnp_only, " attributable to missing NT-proBNP (not selected for, or missing",
        " within, the surplus-sera biomarker sub-study), and ", n_missing_other_core,
        " attributable to other missing core variables.")

cohort <- cohort %>% filter(!is.na(hf_continuum_level))
message("After available core variables needed to assign HF-continuum phenotype",
        " (including NT-proBNP availability where it was the deciding factor): ", nrow(cohort),
        " (excluded ", n_before_core_filter - nrow(cohort), ")")

## CONFIRMED via review: a distinct, independent missingness issue from
## the NT-proBNP step above -- PERMTH_EXM (follow-up time anchored to the
## MEC examination date) can be missing even for participants who are
## eligible for mortality linkage (ELIGSTAT==1) and who have a
## successfully assigned HF-continuum phenotype. This happens because
## NHANES's public-use mortality file tracks two independent follow-up
## clocks: PERMTH_INT (anchored to the household interview date) and
## PERMTH_EXM (anchored to the MEC exam date) -- ELIGSTAT==1 only
## guarantees interview-level linkage eligibility, not that an
## exam-based follow-up time was successfully computed. All primary
## survival/Cox analyses in this project use PERMTH_EXM specifically
## (not PERMTH_INT) because exposure status -- NT-proBNP, phenotype
## assignment -- is only known as of the exam date; anchoring follow-up
## time to the (earlier) interview date would let the observation window
## begin before exposure was actually measured. Confirmed: coxph() was
## silently dropping these participants in every model fit via its
## default na.action, without this ever appearing as an explicit,
## reported exclusion step -- this filter makes that exclusion explicit
## and keeps the reported cohort N consistent with the N actually used in
## every downstream Cox model.
n_before_permth_filter <- nrow(cohort)
cohort <- cohort %>% filter(!is.na(permth_exm))
message("After excluding missing exam-based follow-up time (PERMTH_EXM): ", nrow(cohort),
        " (excluded ", n_before_permth_filter - nrow(cohort), ")")

## Exclusions
cohort <- cohort %>% filter(is.na(ridexprg) | ridexprg != 1)
message("After excluding pregnant at examination: ", nrow(cohort))

message("Final analytic cohort N: ", nrow(cohort), " (excluded ", n_start - nrow(cohort), " of ", n_start, ")")

## NOTE: "missing key covariates after multiple imputation rules" (protocol
## §3, exclusion #4) is intentionally NOT applied here -- that's a
## post-imputation criterion for the main models (Methods §8.2), not part
## of building this base analytic cohort. Complete-case vs. MI sensitivity
## comparison happens in Results §7.

## ---- 4. Survey design object (Results §1 setup) -----------------------
## wtmec6yr was already constructed above (from WTSSCB2Y/WTSSCB4Y, the
## biomarker subsample weights) before the cohort was filtered, so it can
## be used directly here -- no further construction needed.
nhanes_design <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  nest = TRUE,
  data = cohort
)
message("\nSurvey design object created. Design df: ", degf(nhanes_design))

## ============================================================
## RESULTS §2: Distribution of HF-continuum phenotypes
## ============================================================
## ---- 5. Distribution of HF-continuum phenotypes (Results §2) --------
phenotype_dist_unweighted <- cohort %>%
  count(hf_continuum_label) %>%
  mutate(pct_unweighted = round(100 * n / sum(n), 1))

phenotype_dist_weighted <- svymean(~hf_continuum_label, nhanes_design, na.rm = TRUE)
phenotype_dist_weighted_ci <- confint(phenotype_dist_weighted)

round(cbind(estimate = coef(phenotype_dist_weighted) * 100,
            phenotype_dist_weighted_ci * 100), 1)

message("\n---- Distribution of HF-continuum phenotypes ----")
message("Unweighted:")
print(phenotype_dist_unweighted)
message("\nSurvey-weighted (proportion, 95% CI):")
print(round(cbind(estimate = coef(phenotype_dist_weighted) * 100,
                  phenotype_dist_weighted_ci * 100), 1))

## ---- Figure 3: weighted phenotype prevalence stratified by 6 variables ----
## Doesn't need MI (uses observed data only, matching everything in
## Results §1-3). Stratifiers per your Figure 3 spec: age group, sex,
## race/ethnicity, poverty-income ratio, CKD status, diabetes status.
cohort_fig3 <- cohort %>%
  mutate(
    pir_cat = cut(indfmpir, breaks = c(-Inf, 1, 2, 4, Inf),
                  labels = c("<1", "1-2", "2-4", ">4")),
    ckd_status = ifelse(flag_ckd == 1, "CKD present", "No CKD"),
    diabetes_status = ifelse(flag_diabetes == 1, "Diabetes present", "No diabetes"),
    sex_label = ifelse(riagendr == 1, "Male", "Female"),
    race_label = case_when(
      ridreth1 == 1 ~ "Mexican\nAmerican", ridreth1 == 2 ~ "Other\nHispanic",
      ridreth1 == 3 ~ "Non-Hispanic\nWhite", ridreth1 == 4 ~ "Non-Hispanic\nBlack",
      ridreth1 == 5 ~ "Other\nRace", TRUE ~ NA_character_
    )
  )

fig3_design <- svydesign(ids = ~sdmvpsu, strata = ~sdmvstra, weights = ~wtmec6yr,
                         nest = TRUE, data = cohort_fig3)

## Manual (not svyby-based) approach: subset the design to each stratum
## level in turn and compute svymean() directly -- more explicit/robust
## than relying on svyby()'s auto-reshaped output format, which isn't
## independently verified in this project.
compute_strat_prevalence <- function(design, stratvar, stratlabel) {
  dat <- design$variables
  levels_present <- sort(unique(na.omit(dat[[stratvar]])))
  purrr::map_dfr(levels_present, function(lvl) {
    sub_design <- subset(design, dat[[stratvar]] == lvl)
    m <- svymean(~hf_continuum_label, sub_design, na.rm = TRUE)
    tibble(
      stratifier = stratlabel,
      stratum_level = as.character(lvl),
      phenotype = sub("^hf_continuum_label", "", names(m)),
      pct = as.numeric(coef(m)) * 100
    )
  })
}

strat_specs <- list(
  age_band = "Age Group", sex_label = "Sex", race_label = "Race/Ethnicity",
  pir_cat = "Poverty-Income Ratio", ckd_status = "CKD Status",
  diabetes_status = "Diabetes Status"
)

fig3_data <- purrr::imap_dfr(strat_specs, function(label, var) {
  compute_strat_prevalence(fig3_design, var, label)
})

## Explicit facet (subplot) order -- matches your original Figure 3 spec
## order (Age group, Sex, Race/ethnicity, PIR, CKD, Diabetes). Without
## this, facet_wrap() defaults to alphabetical, which doesn't match that
## order (e.g. "CKD Status" would come before "Age Group").
stratifier_order <- c("Age Group", "Sex", "Race/Ethnicity", "Poverty-Income Ratio",
                      "CKD Status", "Diabetes Status")
fig3_data <- fig3_data %>%
  mutate(stratifier = factor(stratifier, levels = stratifier_order))

## Explicit within-facet (x-axis, left-to-right) category order for each
## stratifier. stratum_level came out of compute_strat_prevalence() as
## plain character (as.character(lvl)), which discards any intended order
## (e.g. pir_cat's cut()-defined "<1,1-2,2-4,>4" sequence) and falls back
## to alphabetical -- this is what made PIR's "<"/">" symbols sort
## strangely. All 6 stratifiers' levels are combined into ONE factor
## here; facet_wrap(scales="free_x") automatically shows only the levels
## present in each panel, in this specified order, and drops the rest.
stratum_level_order <- c(
  ## Age Group: chronological (already sorts this way alphabetically by
  ## coincidence, restated explicitly to lock it in regardless)
  "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+",
  ## Sex
  "Male", "Female",
  ## Race/Ethnicity: matches NHANES RIDRETH1's own numeric code order
  ## (1=Mexican American...5=Other Race), not alphabetical
  "Mexican\nAmerican", "Other\nHispanic", "Non-Hispanic\nWhite", "Non-Hispanic\nBlack", "Other\nRace",
  ## Poverty-Income Ratio: low to high, matching the cut() breaks
  "<1", "1-2", "2-4", ">4",
  ## CKD Status: reference (no condition) first
  "No CKD", "CKD present",
  ## Diabetes Status: reference (no condition) first
  "No diabetes", "Diabetes present"
)
fig3_data <- fig3_data %>%
  mutate(stratum_level = factor(stratum_level, levels = stratum_level_order))

## Shortened, print-friendly phenotype labels for the legend (the full
## hf_continuum_label text is long and repeats "HF" awkwardly in a
## legend). Order explicitly as a factor so both the stacking order and
## the legend order follow disease severity (No risk -> Stage C), and the
## fill scale below lines up with this same order.
phenotype_short_labels <- c(
  "No apparent HF risk" = "No apparent risk",
  "Stage A (at risk)" = "Stage A (at risk)",
  "Stage B (pre-HF, biomarker)" = "Stage B (pre-HF)",
  "Stage C (self-reported HF)" = "Stage C (HF)"
)
fig3_data <- fig3_data %>%
  mutate(phenotype = factor(phenotype_short_labels[phenotype],
                            levels = phenotype_short_labels))

## Colors per your own adjustment.
phenotype_colors <- c(
  "No apparent risk"   = "#F8766D",
  "Stage A (at risk)"  = "#7CAE00",
  "Stage B (pre-HF)"   = "#00BFC4",
  "Stage C (HF)"       = "#C77CFF"
)

# fig3_plot <- ggplot2::ggplot(fig3_data, ggplot2::aes(x = stratum_level, y = pct, fill = phenotype)) +
#   ## Fixed proportional width (of the category band) instead of ggplot's
#   ## default -- keeps bars from looking inconsistently thin/thick across
#   ## facets that have different numbers of stratum levels.
#   ggplot2::geom_col(position = "stack", width = 0.65) +
#   ## Percentage label centered within each stacked segment. Segments
#   ## under 3% are left blank (empty string, not just tiny text) --
#   ## labeling every sliver clutters the chart and the numbers become
#   ## unreadable at that size anyway.
#   ggplot2::geom_text(
#     ggplot2::aes(label = ifelse(pct >= 3, paste0(round(pct, 1), "%"), "")),
#     position = ggplot2::position_stack(vjust = 0.5),
#     size = 3.2, color = "#1A1A1A"
#   ) +
#   ## Both x and y scales free per facet now (was free_x only) -- each
#   ## stratifier's panel scales independently on both axes.
#   ggplot2::facet_wrap(~stratifier, scales = "free", ncol = 2) +
#   ggplot2::scale_fill_manual(values = phenotype_colors, drop = FALSE) +
#   ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0, 0.02, 0)) +
#   ggplot2::labs(x = NULL, y = "Weighted prevalence", fill = "HF-continuum phenotype") +
#   ggplot2::theme_minimal(base_size = 12) +
#   ggplot2::theme(
#     legend.position = "bottom",
#     axis.text.x = ggplot2::element_text(size = 11),
#     axis.title.y = ggplot2::element_text(size = 13),
#     panel.grid.major.x = ggplot2::element_blank(),
#     strip.text = ggplot2::element_text(face = "bold", size = 12),
#     plot.title = ggplot2::element_text(face = "bold", size = 12),
#     legend.text = ggplot2::element_text(size = 13), legend.title = ggplot2::element_text(size = 13)
#   )

fig3_data <- fig3_data %>%
  group_by(stratifier) %>%
  mutate(n_categories = n_distinct(stratum_level)) %>%
  ungroup() %>%
  mutate(bar_width = 0.65 * n_categories / max(n_categories))  # scale bar width

fig3_plot <- ggplot2::ggplot(fig3_data, ggplot2::aes(x = stratum_level, y = pct, fill = phenotype)) +
  ggplot2::geom_col(ggplot2::aes(width = bar_width), position = "stack") + 
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(pct >= 3, paste0(round(pct, 1), "%"), "")),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.2, color = "#1A1A1A"
  ) +
  ggplot2::facet_wrap(~stratifier, scales = "free", ncol = 2) +
  ggplot2::scale_fill_manual(values = phenotype_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1), expand = c(0, 0, 0.02, 0)) +
  ggplot2::labs(x = NULL, y = "Weighted prevalence", fill = "HF-continuum phenotype") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(size = 11),
    axis.title.y = ggplot2::element_text(size = 13),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    legend.text = ggplot2::element_text(size = 13), legend.title = ggplot2::element_text(size = 13)
  )

ggplot2::ggsave("data/HF_continuum_Figure3_stratified_prevalence.png", fig3_plot,
                width = 12, height = 9, dpi = 300)
message("Saved: data/HF_continuum_Figure3_stratified_prevalence.png")


## ============================================================
## RESULTS §3: Baseline characteristics across phenotypes
## ============================================================
## ---- 6. Baseline characteristics across phenotypes (Results §3) -----
## Survey-weighted Table 1 by HF-continuum phenotype: weighted means/props
## with design-based tests, PLUS an explicit linear trend-across-stages
## p-value (protocol §8.3 point 4), which svyCreateTableOne's overall test
## does not itself provide.

table1_vars <- c(
  "ridageyr", "riagendr", "ridreth1", "dmdeduc2", "indfmpir", "hiq011", "huq030", "fsdhh",
  "flag_hypertension", "mean_sbp", "mean_dbp", "flag_diabetes", "lbxgh",
  "bmxbmi", "bmxwaist", "dyslipidemia", "smq020", "smq040",
  "flag_ckd", "flag_albuminuria", "egfr", "uacr", "ckd_category",
  "flag_chd_mi", "flag_stroke",
  "ssbnp", "nt_probnp_log", "huq010"
)
factor_vars <- c("riagendr", "ridreth1", "dmdeduc2", "hiq011", "huq030", "fsdhh",
                 "flag_hypertension", "flag_diabetes", "dyslipidemia", "smq020", "smq040",
                 "flag_ckd", "flag_albuminuria", "ckd_category", "flag_chd_mi", "flag_stroke",
                 "huq010")

## huq010 is about to become a factor (below, for Table 1) -- keep a
## numeric copy first, since the trend test treats it as ordinal/
## continuous (1=Excellent...5=Poor) and svyglm() can't do that on a
## factor response.
cohort <- cohort %>% mutate(huq010_numeric = as.numeric(huq010))

cohort <- cohort %>% mutate(across(all_of(factor_vars), as.factor))

## Rebuild the design from the factor-converted cohort -- update() on a
## survey design object is for adding/replacing individual variables via
## name=expression, not for swapping in a whole new data frame, so it's
## safer to just re-run svydesign() with the same ids/strata/weights.
nhanes_design <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  nest = TRUE,
  data = cohort
)

table1 <- svyCreateTableOne(
  vars = table1_vars,
  strata = "hf_continuum_label",
  data = nhanes_design,
  factorVars = factor_vars,
  test = TRUE,
  smd = TRUE
)
table1_printed <- print(table1, smd = TRUE, printToggle = FALSE, showAllLevels = FALSE)

message("\n---- Baseline characteristics across HF-continuum phenotypes (WEIGHTED) ----")
print(table1_printed)

## Unweighted version: same variables/strata, but on the plain cohort data
## frame (no survey design) -- this is what gives ACTUAL SAMPLE COUNTS
## (hundreds/thousands) rather than population-weighted estimates
## (millions). Standard practice is to report unweighted N alongside
## weighted %/means -- this table is what supplies that N.
table1_unweighted <- CreateTableOne(
  vars = table1_vars,
  strata = "hf_continuum_label",
  data = cohort,
  factorVars = factor_vars,
  test = TRUE,
  smd = TRUE
)
table1_unweighted_printed <- print(table1_unweighted, smd = TRUE, printToggle = FALSE, showAllLevels = FALSE)

message("\n---- Baseline characteristics across HF-continuum phenotypes (UNWEIGHTED) ----")
print(table1_unweighted_printed)

## ---- 6b. Explicit linear trend test across ordered stages (Results §3) ----
## Continuous vars: svyglm(var ~ level), trend p = Wald test on the slope.
## Binary/categorical vars: same approach on a 0/1-recoded indicator
## (logistic trend). This directly answers protocol §8.3's "trend tests
## across ordered stages," which svyCreateTableOne's test=TRUE does not.
trend_test_continuous <- function(varname, design) {
  f <- as.formula(paste0(varname, " ~ hf_continuum_level"))
  m <- tryCatch(svyglm(f, design = design), error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  coef(summary(m))["hf_continuum_level", "Pr(>|t|)"]
}

trend_test_binary <- function(varname, design, event_level = "1") {
  design_tmp <- update(design, .tmp_bin = as.numeric(get(varname) == event_level))
  m <- tryCatch(svyglm(.tmp_bin ~ hf_continuum_level, design = design_tmp, family = quasibinomial()),
                error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  coef(summary(m))["hf_continuum_level", "Pr(>|t|)"]
}

## huq010_numeric (self-rated health, 1=Excellent...5=Poor) treated as
## ordinal/continuous for the linear trend test -- standard practice for
## a 5-level ordered scale. Uses the preserved numeric copy since huq010
## itself is now a factor (for Table 1).
continuous_trend_vars <- c("ridageyr", "indfmpir", "mean_sbp", "mean_dbp", "lbxgh",
                           "bmxbmi", "bmxwaist", "egfr", "uacr", "ssbnp", "nt_probnp_log",
                           "huq010_numeric")
binary_trend_vars <- c("flag_hypertension", "flag_diabetes", "dyslipidemia",
                       "flag_ckd", "flag_albuminuria", "flag_chd_mi", "flag_stroke")

trend_results <- bind_rows(
  tibble(variable = continuous_trend_vars,
         trend_p = map_dbl(continuous_trend_vars, trend_test_continuous, design = nhanes_design)),
  tibble(variable = binary_trend_vars,
         trend_p = map_dbl(binary_trend_vars, trend_test_binary, design = nhanes_design))
)

message("\n---- Linear trend p-values across ordered HF-continuum stages ----")
print(trend_results)

## ---- 7. Save Results §1-3 outputs --------------------------------------
dir.create("data", showWarnings = FALSE)
writexl::write_xlsx(cohort, "data/HF_continuum_analytic_cohort.xlsx")

## Bypass as.data.frame()'s automatic row-name sanitization (make.names()
## turns "   High risk" into "X...High.risk" -- spaces become periods, and
## an "X" gets prepended since the result no longer starts with a letter).
## Grabbing rownames(table1_printed) directly keeps the original readable
## text; trimws() removes tableone's leading indentation spaces.
table1_df <- as.data.frame(table1_printed, stringsAsFactors = FALSE)
table1_df <- tibble::tibble(Variable = trimws(rownames(table1_printed)), table1_df)

table1_unweighted_df <- as.data.frame(table1_unweighted_printed, stringsAsFactors = FALSE)
table1_unweighted_df <- tibble::tibble(
  Variable = trimws(rownames(table1_unweighted_printed)), table1_unweighted_df
)

## ---- Add a readable "label" column ------------------------------------
## tableone's raw row text mixes several formats: a continuous variable is
## one row ("ridageyr (mean (SD))"); a 2-level factor collapses to one row
## ("hiq011 = 2 (%)"); a >2-level factor is a header row ("dmdeduc2 (%)")
## followed by one bare-code row per category ("1", "2", ...) with NO
## variable name repeated -- so labeling those requires remembering which
## header row we most recently passed. add_readable_labels() below does
## that by walking down the rows in order.

var_labels <- c(
  ridageyr = "Age, years", riagendr = "Sex", ridreth1 = "Race/Ethnicity",
  dmdeduc2 = "Education", indfmpir = "Poverty-Income Ratio",
  hiq011 = "Health Insurance Coverage", huq030 = "Usual Source of Care",
  fsdhh = "Household Food Security", flag_hypertension = "Hypertension",
  mean_sbp = "Systolic BP, mmHg", mean_dbp = "Diastolic BP, mmHg",
  flag_diabetes = "Diabetes", lbxgh = "HbA1c, %", bmxbmi = "BMI, kg/m2",
  bmxwaist = "Waist Circumference, cm", dyslipidemia = "Dyslipidemia",
  smq020 = "Ever Smoked >=100 Cigarettes", smq040 = "Current Smoking Frequency",
  flag_ckd = "Chronic Kidney Disease (self-report)",
  flag_albuminuria = "Albuminuria", egfr = "eGFR, mL/min/1.73m2",
  uacr = "Urine Albumin-Creatinine Ratio, mg/g", ckd_category = "CKD Risk Category",
  flag_chd_mi = "Coronary Heart Disease / MI", flag_stroke = "Stroke",
  ssbnp = "NT-proBNP, pg/mL", nt_probnp_log = "Log NT-proBNP", n = "Sample Size",
  huq010 = "Self-Rated Health"
)

## (variable, category code) -> readable category label, for multi-level
## factors' sub-rows and for the single value shown by 2-level factors.
category_labels <- list(
  riagendr = c(`1` = "Male", `2` = "Female"),
  ridreth1 = c(`1` = "Mexican American", `2` = "Other Hispanic",
               `3` = "Non-Hispanic White", `4` = "Non-Hispanic Black",
               `5` = "Other Race - Including Multi-Racial"),
  dmdeduc2 = c(`1` = "Less than 9th grade", `2` = "9-11th grade",
               `3` = "High school grad/GED", `4` = "Some college/AA degree",
               `5` = "College graduate+"),
  hiq011 = c(`1` = "Insured", `2` = "Uninsured"),
  huq030 = c(`1` = "Has a usual place", `2` = "No usual place",
             `3` = "More than one place"),
  fsdhh = c(`1` = "Full food security", `2` = "Marginal food security",
            `3` = "Low food security", `4` = "Very low food security"),
  huq010 = c(`1` = "Excellent", `2` = "Very good", `3` = "Good",
             `4` = "Fair", `5` = "Poor"),
  flag_hypertension = c(`0` = "No", `1` = "Yes"),
  flag_diabetes = c(`0` = "No", `1` = "Yes"),
  dyslipidemia = c(`0` = "No", `1` = "Yes"),
  smq020 = c(`1` = "Yes", `2` = "No"),
  smq040 = c(`1` = "Every day", `2` = "Some days", `3` = "Not at all"),
  flag_ckd = c(`0` = "No", `1` = "Yes"),
  flag_albuminuria = c(`0` = "No", `1` = "Yes"),
  flag_chd_mi = c(`0` = "No", `1` = "Yes"),
  flag_stroke = c(`0` = "No", `1` = "Yes")
  ## ckd_category is NOT listed here -- it's already stored as readable
  ## text ("Low risk", "High risk", etc.), not a numeric code, so its
  ## sub-rows are already human-readable as-is.
)

add_readable_labels <- function(df) {
  var_names_sorted <- names(var_labels)[order(-nchar(names(var_labels)))]  # longest first, avoids partial-name collisions
  current_var <- NA_character_
  labels_out <- character(nrow(df))
  
  for (i in seq_len(nrow(df))) {
    row_text <- df$Variable[i]
    
    ## Try to match this row to a known variable name as its header/own row.
    matched_var <- var_names_sorted[startsWith(row_text, var_names_sorted)][1]
    
    if (!is.na(matched_var)) {
      current_var <- matched_var
      suffix <- trimws(sub(paste0("^", matched_var), "", row_text))
      if (grepl("^= ", suffix)) {
        ## 2-level factor collapsed to one row, e.g. "hiq011 = 2 (%)"
        code <- trimws(gsub("=|\\(%\\)", "", suffix))
        cat_lab <- category_labels[[matched_var]][code]
        labels_out[i] <- if (!is.na(cat_lab)) paste0(var_labels[matched_var], ": ", cat_lab)
        else var_labels[matched_var]
      } else {
        ## Continuous var row, or multi-level factor HEADER row
        labels_out[i] <- unname(var_labels[matched_var])
      }
    } else if (!is.na(current_var) && row_text %in% names(category_labels[[current_var]])) {
      ## Bare category-code sub-row belonging to the last-seen multi-level header
      labels_out[i] <- category_labels[[current_var]][row_text]
    } else if (!is.na(current_var) && current_var == "ckd_category") {
      ## ckd_category's sub-rows are already readable text -- pass through
      labels_out[i] <- row_text
    } else {
      labels_out[i] <- NA_character_   # unmatched -- see verification print below
    }
  }
  df$label <- labels_out
  df %>% select(Variable, label, everything())
}

table1_df <- add_readable_labels(table1_df)
table1_unweighted_df <- add_readable_labels(table1_unweighted_df)

unmatched <- table1_df %>% filter(is.na(label) & !Variable %in% c("", "p", "test", "SMD"))
if (nrow(unmatched) > 0) {
  message("\n*** Rows that could NOT be auto-labeled -- check these manually: ***")
  print(unmatched %>% select(Variable))
} else {
  message("\nAll Table 1 rows successfully labeled.")
}

writexl::write_xlsx(
  list(
    "Baseline (Weighted)" = table1_df,
    "Baseline (Unweighted)" = table1_unweighted_df,
    "Trend Tests" = trend_results,
    "Phenotype Dist (Unwtd)" = phenotype_dist_unweighted
  ),
  "data/HF_continuum_baseline_characteristics.xlsx"
)

message("\nSaved: data/HF_continuum_analytic_cohort.xlsx")
message("Saved: data/HF_continuum_baseline_characteristics.xlsx (weighted + unweighted sheets)")


## ============================================================
## NT-proBNP: 按phenotype分组，同时输出 arithmetic mean 和 geometric mean
## （几何均值 = exp(log值的加权算术均值)，CI在log尺度上算好再exp回来，
## 这是右偏分布变量算几何均值置信区间的标准做法，不能直接在原始尺度上
## 套算术均值的CI公式）
## ============================================================
## 前提：nhanes_design（或您自己命名的survey design对象）、cohort 已经
## 在当前session中存在，cohort里有 ssbnp（原始NT-proBNP）和
## nt_probnp_log（log转换后的值）。

strata_levels <- levels(cohort$hf_continuum_label)

compute_ntprobnp_summary <- function(design, lvl) {
  sub_design <- subset(design, design$variables$hf_continuum_label == lvl)
  
  ## ---- Arithmetic mean（原始尺度，标签准确反映实际算法）----
  m_arith <- svymean(~ssbnp, sub_design, na.rm = TRUE)
  ci_arith <- confint(m_arith)
  
  ## ---- Geometric mean：在log尺度上算均值+CI，最后再exp回原始单位 ----
  ## exp(加权算术均值(log(x))) = 几何均值；CI同理先在log尺度算好再exp，
  ## 不能直接对原始尺度的CI做log/exp转换（那样是错的）。
  m_log <- svymean(~nt_probnp_log, sub_design, na.rm = TRUE)
  ci_log <- confint(m_log)
  geo_mean <- exp(coef(m_log)[1])
  geo_ci_lower <- exp(ci_log[1, 1])
  geo_ci_upper <- exp(ci_log[1, 2])
  
  tibble(
    hf_continuum_label = lvl,
    arithmetic_mean = round(coef(m_arith)[1], 2),
    arithmetic_ci_lower = round(ci_arith[1, 1], 2),
    arithmetic_ci_upper = round(ci_arith[1, 2], 2),
    geometric_mean = round(geo_mean, 2),
    geometric_ci_lower = round(geo_ci_lower, 2),
    geometric_ci_upper = round(geo_ci_upper, 2)
  )
}

ntprobnp_summary <- purrr::map_dfr(strata_levels, function(lvl) {
  compute_ntprobnp_summary(nhanes_design, lvl)
})

message("\n---- NT-proBNP: Arithmetic mean vs Geometric mean, by phenotype ----")
print(ntprobnp_summary)

## ---- 格式化成表格里可以直接用的文字，两个版本都给 ----
ntprobnp_formatted <- ntprobnp_summary %>%
  mutate(
    label_arithmetic = sprintf("%.1f (%.1f\u2013%.1f)", arithmetic_mean, arithmetic_ci_lower, arithmetic_ci_upper),
    label_geometric = sprintf("%.1f (%.1f\u2013%.1f)", geometric_mean, geometric_ci_lower, geometric_ci_upper)
  ) %>%
  select(hf_continuum_label, label_arithmetic, label_geometric)

print(ntprobnp_formatted)

## ============================================================
## Table 1 重新格式化：unweighted N + weighted % (95% CI) / weighted mean (95% CI)
## ============================================================
## 前提：nhanes_design（或您自己命名的survey design对象）、cohort、
## table1_vars、factor_vars 都已经在当前session里存在。

strata_levels <- levels(cohort$hf_continuum_label)

## ---- 分类变量：每个类别 -> "未加权N (加权% [95% CI])" ----
format_categorical_var <- function(var, design, cohort) {
  purrr::map_dfr(strata_levels, function(lvl) {
    sub_design <- subset(design, design$variables$hf_continuum_label == lvl)
    sub_cohort <- cohort %>% filter(hf_continuum_label == lvl)
    
    ## 加权百分比+CI：对这个变量的每个类别分别算
    form <- as.formula(paste0("~", var))
    m <- svymean(form, sub_design, na.rm = TRUE)
    ci <- confint(m)
    
    var_levels <- sub(paste0("^", var), "", names(m))
    weighted_pct <- as.numeric(coef(m)) * 100
    ci_lower <- ci[, 1] * 100
    ci_upper <- ci[, 2] * 100
    
    ## 未加权N：这个类别在这个phenotype组里的实际人数
    unweighted_n <- sapply(var_levels, function(lv) {
      sum(as.character(sub_cohort[[var]]) == lv, na.rm = TRUE)
    })
    
    tibble(
      variable = var,
      hf_continuum_label = lvl,
      level = var_levels,
      formatted = sprintf("%d (%.1f%% [%.1f\u2013%.1f])",
                          unweighted_n, weighted_pct, ci_lower, ci_upper)
    )
  })
}

## ---- 连续变量：直接 "加权均值 (95% CI)" ----
format_continuous_var <- function(var, design, cohort) {
  purrr::map_dfr(strata_levels, function(lvl) {
    sub_design <- subset(design, design$variables$hf_continuum_label == lvl)
    
    form <- as.formula(paste0("~", var))
    m <- svymean(form, sub_design, na.rm = TRUE)
    ci <- confint(m)
    
    tibble(
      variable = var,
      hf_continuum_label = lvl,
      level = NA_character_,
      formatted = sprintf("%.1f (%.1f\u2013%.1f)", coef(m)[1], ci[1, 1], ci[1, 2])
    )
  })
}

## ---- 对table1_vars里所有变量跑一遍，自动区分分类/连续 ----
table1_formatted_long <- purrr::map_dfr(table1_vars, function(var) {
  if (var %in% factor_vars) {
    format_categorical_var(var, nhanes_design, cohort)
  } else {
    format_continuous_var(var, nhanes_design, cohort)
  }
})

## ---- 转成宽格式：行=变量(类别)，列=4个phenotype ----
table1_formatted_wide <- table1_formatted_long %>%
  mutate(row_label = ifelse(is.na(level), variable, paste0(variable, ": ", level))) %>%
  select(row_label, hf_continuum_label, formatted) %>%
  tidyr::pivot_wider(names_from = hf_continuum_label, values_from = formatted)

print(table1_formatted_wide, n = Inf)

writexl::write_xlsx(table1_formatted_wide, "data/HF_continuum_Table1_weighted_formatted.xlsx")
message("已保存: data/HF_continuum_Table1_weighted_formatted.xlsx")