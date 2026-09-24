## ============================================================
## Phenotype definition
## ============================================================

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
hf_continuum_data <- readxl::read_excel("./data/HF_continuum_cleaned_dataset.xlsx")
message("Loaded cleaned hf_continuum_data: ", nrow(hf_continuum_data), " rows x ",
        ncol(hf_continuum_data), " cols")

## ---- 2. Derived variables --------------------------------------------

## 3.1 eGFR
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

## 3.5b Cardiovascular/heart-disease mortality
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

## 3.8 HF-continuum phenotype (0-3) 
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


## ============================================================
## Analytic cohort construction
## ============================================================
## ---- 3. Analytic cohort construction (Results §1) -------------------
## Applies protocol Section 3 inclusion/exclusion criteria in order
n_start <- nrow(hf_continuum_data)
message("\n---- Analytic cohort construction ----")
message("Starting N (all loaded records): ", n_start)

cohort <- hf_continuum_data

cohort <- cohort %>% filter(ridageyr >= 20)
message("After age >= 20: ", nrow(cohort))

cohort <- cohort %>% filter(!is.na(mcq160b))
message("After completed medical conditions questionnaire (non-missing mcq160b): ", nrow(cohort))


cohort <- cohort %>%
  mutate(
    wtmec6yr = case_when(
      !is.na(wtsscb4y) ~ wtsscb4y * (2 / 3),
      !is.na(wtsscb2y) ~ wtsscb2y * (1 / 3),
      TRUE ~ NA_real_
    )
  )

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

n_before_permth_filter <- nrow(cohort)
cohort <- cohort %>% filter(!is.na(permth_exm))
message("After excluding missing exam-based follow-up time (PERMTH_EXM): ", nrow(cohort),
        " (excluded ", n_before_permth_filter - nrow(cohort), ")")

## Exclusions
cohort <- cohort %>% filter(is.na(ridexprg) | ridexprg != 1)
message("After excluding pregnant at examination: ", nrow(cohort))

message("Final analytic cohort N: ", nrow(cohort), " (excluded ", n_start - nrow(cohort), " of ", n_start, ")")


## ---- 4. Survey design object -----------------------
nhanes_design <- svydesign(
  ids = ~sdmvpsu,
  strata = ~sdmvstra,
  weights = ~wtmec6yr,
  nest = TRUE,
  data = cohort
)
message("\nSurvey design object created. Design df: ", degf(nhanes_design))

## ============================================================
## Distribution of HF-continuum phenotypes
## ============================================================
## ---- 5. Distribution of HF-continuum phenotypes--------
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

stratifier_order <- c("Age Group", "Sex", "Race/Ethnicity", "Poverty-Income Ratio",
                      "CKD Status", "Diabetes Status")
fig3_data <- fig3_data %>%
  mutate(stratifier = factor(stratifier, levels = stratifier_order))


stratum_level_order <- c(
  ## Age Group
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

cohort <- cohort %>% mutate(huq010_numeric = as.numeric(huq010))

cohort <- cohort %>% mutate(across(all_of(factor_vars), as.factor))

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

## ---- 6b. Explicit linear trend test across ordered stages ----
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


table1_df <- as.data.frame(table1_printed, stringsAsFactors = FALSE)
table1_df <- tibble::tibble(Variable = trimws(rownames(table1_printed)), table1_df)

table1_unweighted_df <- as.data.frame(table1_unweighted_printed, stringsAsFactors = FALSE)
table1_unweighted_df <- tibble::tibble(
  Variable = trimws(rownames(table1_unweighted_printed)), table1_unweighted_df
)

## ---- Add a readable "label" column ------------------------------------
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
## NT-proBNP
## ============================================================
strata_levels <- levels(cohort$hf_continuum_label)

compute_ntprobnp_summary <- function(design, lvl) {
  sub_design <- subset(design, design$variables$hf_continuum_label == lvl)
  
  ## ---- Arithmetic mean----
  m_arith <- svymean(~ssbnp, sub_design, na.rm = TRUE)
  ci_arith <- confint(m_arith)
  
  ## ---- Geometric mean----
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

ntprobnp_formatted <- ntprobnp_summary %>%
  mutate(
    label_arithmetic = sprintf("%.1f (%.1f\u2013%.1f)", arithmetic_mean, arithmetic_ci_lower, arithmetic_ci_upper),
    label_geometric = sprintf("%.1f (%.1f\u2013%.1f)", geometric_mean, geometric_ci_lower, geometric_ci_upper)
  ) %>%
  select(hf_continuum_label, label_arithmetic, label_geometric)

print(ntprobnp_formatted)

## ============================================================
## Table 1 unweighted N + weighted % (95% CI) / weighted mean (95% CI)
## ============================================================
strata_levels <- levels(cohort$hf_continuum_label)

format_categorical_var <- function(var, design, cohort) {
  purrr::map_dfr(strata_levels, function(lvl) {
    sub_design <- subset(design, design$variables$hf_continuum_label == lvl)
    sub_cohort <- cohort %>% filter(hf_continuum_label == lvl)

    form <- as.formula(paste0("~", var))
    m <- svymean(form, sub_design, na.rm = TRUE)
    ci <- confint(m)
    
    var_levels <- sub(paste0("^", var), "", names(m))
    weighted_pct <- as.numeric(coef(m)) * 100
    ci_lower <- ci[, 1] * 100
    ci_upper <- ci[, 2] * 100

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

table1_formatted_long <- purrr::map_dfr(table1_vars, function(var) {
  if (var %in% factor_vars) {
    format_categorical_var(var, nhanes_design, cohort)
  } else {
    format_continuous_var(var, nhanes_design, cohort)
  }
})

table1_formatted_wide <- table1_formatted_long %>%
  mutate(row_label = ifelse(is.na(level), variable, paste0(variable, ": ", level))) %>%
  select(row_label, hf_continuum_label, formatted) %>%
  tidyr::pivot_wider(names_from = hf_continuum_label, values_from = formatted)

print(table1_formatted_wide, n = Inf)

writexl::write_xlsx(table1_formatted_wide, "data/HF_continuum_Table1_weighted_formatted.xlsx")
message("Saved: data/HF_continuum_Table1_weighted_formatted.xlsx")