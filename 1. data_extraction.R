library(nhanesdata)
library(nhanesA)
library(dplyr)
library(purrr)
library(stringr)
library(writexl)
library(readxl)

## Analytic cycles: NHANES 1999-2000, 2001-2002, 2003-2004
## In nhanesdata's `year` column these are 1999, 2001, 2003
analytic_years <- c(1999, 2001, 2003)


## ---- 1. Variable map ------------------------------------------------
## Domain                    Dataset      Variables                              Notes
## ------------------------  -----------  -------------------------------------  --------------------------------
## Identifiers/design         demo         seqn, sdmvpsu, sdmvstra,
##                                         wtint2yr, wtmec2yr                     survey design + weights
## Demographics                demo         ridageyr, riagendr, ridreth1,
##                                         dmdeduc2, indfmpir, ridexprg           age/sex/race/educ/PIR/pregnancy
## Insurance                   (see section 2C)  hiq011                        loaded as HID010 for 1999-2004, renamed to hiq011
## Usual source of care        huq          huq030                                routine place for healthcare
## Food security                (see section 2C)  fsdhh                         2003-2004=FSDHH; 1999-2002=HHFDSEC
## HF / CVD history (Q)        mcq          mcq160b (CHF), mcq160c (CHD), mcq160d (angina),
##                                         mcq160e (MI), mcq160f (stroke)        yes/no items only -- see below for age-at-dx
## Age at HF/CVD diagnosis      (see section 2C)  mcd180b-f                     loaded via nhanesA; ALL sourced as
##                                                                                MCQ180B/C/D/E/F for 1999-2004 (MCD180*
##                                                                                prefix only exists in later cycles)
## Diabetes (Q)                 diq          diq010
## Hypertension/Dyslipidemia(Q) bpq          bpq020, bpq050a, bpq080 (self-report high cholesterol)
## Kidney disease (Q)           (see section 2C)  kiq022                        2001-2004=KIQ022; 1999-2000=KIQ020
## Smoking (Q)                  smq          smq020, smq040
## Physical activity (Q)        paq          pad200, paq180                       consistent across 1999-2004 specifically
## Blood pressure (exam)        bpx          bpxsy1-bpxsy4, bpxdi1-bpxdi4
## Body measures (exam)         bmx          bmxbmi, bmxwaist
## NT-proBNP (special/surplus)  ssbnp        ssbnp, ssbnpl, sspris                loaded separately in section 3
##                                                                                 (single pooled file, nhanesA::nhanes("SSBNP_A"))
##
## NOT in var_map -- nhanesdata has 1999-2004 coverage gaps
## and/or these variables were literally renamed for later cycles, so both
## groups are loaded robustly via nhanesA::nhanesSearchVarName() instead,
## which looks up the exact correct table name/year per variable directly
## from CDC (section 2B for labs, section 2C for insurance/age-at-dx/food
## security/kidney disease):
##   ghb (lbxgh) | biopro (lbxscr, lbxsal) | alb_cr (urxuma, urxucr)
##   trigly (lbxtr, lbdldl) | tchol (lbxtc, lbdhdd)
##   hiq011 (sourced as HID010) | mcd180b-f (all sourced as MCQ180B/C/D/E/F)
##   fsdhh (sourced as FSDHH or HHFDSEC) | kiq022 (sourced as KIQ022 or KIQ020)

var_map <- list(
  demo   = c("seqn", "sdmvpsu", "sdmvstra", "wtint2yr", "wtmec2yr",
             "ridageyr", "riagendr", "ridreth1", "dmdeduc2", "indfmpir", "ridexprg"),
  huq    = c("seqn", "huq030", "huq010"),
  mcq    = c("seqn", "mcq160b", "mcq160c", "mcq160d", "mcq160e", "mcq160f"),
  diq    = c("seqn", "diq010"),
  bpq    = c("seqn", "bpq020", "bpq050a", "bpq080"),
  smq    = c("seqn", "smq020", "smq040"),
  paq    = c("seqn", "pad200", "paq180"),
  bpx    = c("seqn", "bpxsy1", "bpxsy2", "bpxsy3", "bpxsy4",
             "bpxdi1", "bpxdi2", "bpxdi3", "bpxdi4"),
  bmx    = c("seqn", "bmxbmi", "bmxwaist")
)

## ---- 2. Load each dataset, filter to 1999-2004, keep needed vars ----
safe_read <- function(dataset_name, keep_vars) {
  message("Loading: ", dataset_name)
  df <- tryCatch(nhanesdata::read_nhanes(dataset_name), error = function(e) {
    warning("Could not read '", dataset_name, "' from nhanesdata: ", conditionMessage(e),
            "\n  -> try nhanesdata::term_search() to find the correct dataset name.")
    NULL
  })
  if (is.null(df)) return(NULL)
  
  df <- df %>% filter(year %in% analytic_years)
  
  keep_vars <- intersect(keep_vars, names(df))
  missing <- setdiff(keep_vars, names(df))
  if (length(missing) > 0) {
    message("  Note: not found in ", dataset_name, ": ", paste(missing, collapse = ", "))
  }
  df %>% select(seqn, year, all_of(setdiff(keep_vars, c("seqn", "year"))))
}

person_level_list <- purrr::imap(var_map, ~ safe_read(.y, .x))
person_level_list <- purrr::compact(person_level_list)   # drop any that failed to load

## ---- 2B. Lab variables via nhanesA (robust to 1999-2004 legacy naming) ----
manual_table_overrides <- list(
  biopro = list(`2001` = "L40_B")   # confirmed to exist on CDC; scrape missed it
)

## Which variables to pull for each lab domain (this got accidentally
## dropped in an earlier edit -- restored here).
lab_groups <- list(
  ghb    = c("lbxgh"),
  biopro = c("lbxscr", "lbxsal"),
  alb_cr = c("urxuma", "urxucr"),
  trigly = c("lbxtr", "lbdldl"),
  tchol  = c("lbxtc", "lbdhdd")
)

var_aliases <- list(
  lbxscr = c("lbxscr", "lbdscr"),          # L40_B (2001-2002) uses LBDSCR instead of LBXSCR
  lbdhdd = c("lbdhdd", "lbdhdl", "lbxhdd") # 1999-2000 & 2001-2002 = LBDHDL (two-method-
  # combined, bias-corrected); 2003-2004 = LBXHDD
  # (direct-method-only, uncorrected -- bias was
  # already acceptable that cycle). All three
  # variable names refer to the same underlying
  # concept (this respondent's HDL-cholesterol,
  # mg/dL) -- confirmed against CDC's variable
  # index (LBXHDD / l13_c / 2003-2004).
)

resolve_var <- function(target, df_names) {
  ## Returns the actual column name in df_names matching target or one of
  ## its known aliases, or NA if none found.
  candidates <- c(target, var_aliases[[target]])
  hit <- candidates[candidates %in% df_names]
  if (length(hit) == 0) NA_character_ else hit[1]
}

extract_vars <- function(df, varnames, tbl) {
  ## Resolves each requested variable (checking aliases), renames to the
  ## canonical (target) name, and reports if a table loaded but matched
  ## nothing -- printing its columns so a new alias is easy to spot.
  resolved <- purrr::map_chr(varnames, ~ resolve_var(.x, names(df)))
  found <- !is.na(resolved)
  if (!any(found)) {
    message("    NOTE: none of [", paste(varnames, collapse = ", "), "] (or known aliases) found in ",
            tbl, ". Its columns are: ", paste(names(df), collapse = ", "))
    return(NULL)
  }
  out <- df[, resolved[found], drop = FALSE]
  names(out) <- varnames[found]
  out[] <- lapply(out, function(x) suppressWarnings(as.numeric(as.character(x))))
  out
}

load_lab_group <- function(varnames, label, ystart = 1999, ystop = 2004) {
  primary <- toupper(varnames[1])
  info <- tryCatch(
    nhanesA::nhanesSearchVarName(primary, ystart = ystart, ystop = ystop, namesonly = FALSE),
    error = function(e) NULL
  )
  if (is.null(info) || nrow(info) == 0) {
    warning("No tables found containing '", primary, "' for ", ystart, "-", ystop)
    info <- data.frame(Data.File.Name = character(0), Begin.Year = numeric(0))
  } else {
    ## Exclude restricted-access (RDC only) and "second exam session"
    ## replicate files, which can duplicate SEQN within the same year.
    if ("Use.Constraints" %in% names(info)) {
      info <- info[is.na(info$Use.Constraints) | info$Use.Constraints == "None", ]
    }
    if ("Data.File.Description" %in% names(info)) {
      info <- info[!grepl("second", info$Data.File.Description, ignore.case = TRUE), ]
    }
  }
  
  out <- purrr::map_dfr(seq_len(nrow(info)), function(i) {
    tbl <- info$Data.File.Name[i]
    yr  <- suppressWarnings(as.numeric(info$Begin.Year[i]))
    message("  Loading ", tbl, " (", yr, ") for: ", paste(varnames, collapse = ", "))
    df <- tryCatch(nhanesA::nhanes(tbl), error = function(e) {
      warning("  Could not load ", tbl, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(df)) return(NULL)
    names(df) <- tolower(names(df))
    keep <- extract_vars(df, tolower(varnames), tbl)
    if (is.null(keep)) return(NULL)
    dplyr::bind_cols(seqn = as.numeric(df$seqn), year = yr, keep)
  })
  
  ## Apply confirmed manual overrides for years the scrape is known to miss.
  overrides <- manual_table_overrides[[label]]
  years_present <- if (nrow(out) > 0) unique(out$year) else numeric(0)
  for (yr_chr in names(overrides)) {
    yr <- as.numeric(yr_chr)
    if (yr %in% years_present) next
    tbl <- overrides[[yr_chr]]
    message("  Patching ", label, "/", yr, " with confirmed table: ", tbl)
    df <- tryCatch(nhanesA::nhanes(tbl), error = function(e) NULL)
    if (is.null(df)) next
    names(df) <- tolower(names(df))
    keep <- extract_vars(df, tolower(varnames), tbl)
    if (is.null(keep)) next
    patch <- dplyr::bind_cols(seqn = as.numeric(df$seqn), year = yr, keep)
    out <- dplyr::bind_rows(out, patch)
  }
  
  if (nrow(out) == 0) return(NULL)
  out <- out %>% distinct(seqn, year, .keep_all = TRUE)   # safety net against residual duplicates
  
  missing_years <- setdiff(c(ystart, ystart + 2, ystart + 4), unique(out$year))
  if (length(missing_years) > 0) {
    warning("Lab group '", label, "' is missing year(s): ", paste(missing_years, collapse = ", "),
            " -- verify manually against CDC's data file pages; this may need a manual_table_overrides entry.")
  }
  
  out
}

lab_list <- purrr::imap(lab_groups, function(vars, label) {
  message("Loading lab group: ", label, " (", paste(vars, collapse = ", "), ")")
  load_lab_group(vars, label)
})
lab_list <- purrr::compact(lab_list)

person_level_list <- c(person_level_list, lab_list)

## ---- 2C. Insurance + MCQ age-at-diagnosis variables via nhanesA -------
legacy_vars <- list(
  hiq011  = "HID010",
  mcd180b = "MCQ180B",
  mcd180c = "MCQ180C",
  mcd180d = "MCQ180D",
  mcd180e = "MCQ180E",
  mcd180f = "MCQ180F",
  fsdhh   = c("FSDHH", "HHFDSEC"),   # 2003-2004 = FSDHH; 1999-2000 & 2001-2002 = HHFDSEC
  kiq022  = c("KIQ022", "KIQ020")    # 2001-2004 = KIQ022; 1999-2000 = KIQ020
)

manual_table_overrides_legacy <- list(
  fsdhh = list(`1999` = "FSQ", `2001` = "FSQ_B")   # confirmed: both contain HHFDSEC
)

load_renamed_var <- function(source_names, output_name, ystart = 1999, ystop = 2004) {
  ## Try every candidate name, since some variables use a different name in
  ## different cycles rather than one uniform alternate name.
  info <- dplyr::bind_rows(purrr::map(source_names, function(sn) {
    tryCatch(
      nhanesA::nhanesSearchVarName(sn, ystart = ystart, ystop = ystop, namesonly = FALSE),
      error = function(e) NULL
    )
  }))
  if (nrow(info) == 0) {
    warning("No tables found containing any of [", paste(source_names, collapse = ", "),
            "] for ", ystart, "-", ystop)
    info <- data.frame(Data.File.Name = character(0), Begin.Year = numeric(0))
  } else {
    if ("Use.Constraints" %in% names(info)) {
      info <- info[is.na(info$Use.Constraints) | info$Use.Constraints == "None", ]
    }
    if ("Data.File.Description" %in% names(info)) {
      info <- info[!grepl("second", info$Data.File.Description, ignore.case = TRUE), ]
    }
    info <- info[!duplicated(info$Data.File.Name), ]   # same table can match >1 candidate name
  }
  
  out <- purrr::map_dfr(seq_len(nrow(info)), function(i) {
    tbl <- info$Data.File.Name[i]
    yr  <- suppressWarnings(as.numeric(info$Begin.Year[i]))
    message("  Loading ", tbl, " (", yr, ") for: ", paste(source_names, collapse = "/"), " -> ", output_name)
    df <- tryCatch(nhanesA::nhanes(tbl), error = function(e) {
      warning("  Could not load ", tbl, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(df)) return(NULL)
    names(df) <- tolower(names(df))
    sn_hit <- tolower(source_names)[tolower(source_names) %in% names(df)]
    if (length(sn_hit) == 0) {
      message("    NOTE: none of [", paste(source_names, collapse = ", "), "] found in ", tbl,
              ". Its columns are: ", paste(names(df), collapse = ", "))
      return(NULL)
    }
    sn <- sn_hit[1]
    df %>%
      transmute(seqn = as.numeric(seqn), year = yr, !!output_name := as.character(.data[[sn]]))
  })
  
  ## Apply confirmed manual overrides for years the scrape is known to miss.
  overrides <- manual_table_overrides_legacy[[output_name]]
  years_present <- if (nrow(out) > 0) unique(out$year) else numeric(0)
  for (yr_chr in names(overrides)) {
    yr <- as.numeric(yr_chr)
    if (yr %in% years_present) next
    tbl <- overrides[[yr_chr]]
    message("  Patching ", output_name, "/", yr, " with confirmed table: ", tbl)
    df <- tryCatch(nhanesA::nhanes(tbl), error = function(e) NULL)
    if (is.null(df)) next
    names(df) <- tolower(names(df))
    sn_hit <- tolower(source_names)[tolower(source_names) %in% names(df)]
    if (length(sn_hit) == 0) next
    sn <- sn_hit[1]
    patch <- df %>% transmute(seqn = as.numeric(seqn), year = yr, !!output_name := as.character(.data[[sn]]))
    out <- dplyr::bind_rows(out, patch)
  }
  
  if (nrow(out) == 0) return(NULL)
  out <- out %>% distinct(seqn, year, .keep_all = TRUE)
  
  missing_years <- setdiff(c(ystart, ystart + 2, ystart + 4), unique(out$year))
  if (length(missing_years) > 0) {
    warning("'", output_name, "' is missing year(s): ", paste(missing_years, collapse = ", "),
            " -- verify manually against CDC's data file pages.")
  }
  
  out
}

legacy_list <- purrr::imap(legacy_vars, function(source_names, output_name) {
  message("Loading legacy-named variable: ", paste(source_names, collapse = "/"), " -> ", output_name)
  load_renamed_var(source_names, output_name)
})
legacy_list <- purrr::compact(legacy_list)

person_level_list <- c(person_level_list, legacy_list)

## ---- 3. NT-proBNP: load directly via nhanesA -----------------------
message("Loading NT-proBNP (pooled 1999-2004) from CDC: SSBNP_A")
ssbnp <- tryCatch(nhanesA::nhanes("SSBNP_A"), error = function(e) {
  warning("Could not load SSBNP_A via nhanesA: ", conditionMessage(e))
  NULL
})

if (!is.null(ssbnp)) {
  names(ssbnp) <- tolower(names(ssbnp))
  ssbnp <- ssbnp %>%
    select(seqn, any_of(c("ssbnp", "ssbnpl", "sspris", "wtsscb2y", "wtsscb4y"))) %>%
    mutate(seqn = as.numeric(seqn))
  message("NT-proBNP loaded: ", nrow(ssbnp), " rows (pooled across 1999-2004)")

  if (!all(c("wtsscb2y", "wtsscb4y") %in% names(ssbnp))) {
    warning("wtsscb2y/wtsscb4y not found in SSBNP_A as returned by nhanesA -- ",
            "the biomarker subsample weight construction downstream will fail.",
            " Check nhanesA::nhanesTableVars('SSBNP_A') for the current column names.")
  } else {
    n_4y <- sum(!is.na(ssbnp$wtsscb4y))
    n_2y <- sum(!is.na(ssbnp$wtsscb2y))
    n_both <- sum(!is.na(ssbnp$wtsscb4y) & !is.na(ssbnp$wtsscb2y))
    message("  wtsscb4y non-missing: ", n_4y, " | wtsscb2y non-missing: ", n_2y,
            " | both non-missing (should be 0): ", n_both)
  }
} else {
  warning("NT-proBNP could not be loaded -- check manually before proceeding.",
          " This is the central biomarker for the project's Stage B/pre-HF definition.")
}

## ---- 4. Merge all year+seqn-keyed person-level files -----------------
nhanes_analytic <- purrr::reduce(
  person_level_list,
  ~ dplyr::full_join(.x, .y, by = c("seqn", "year"), relationship = "one-to-one")
)

if (!is.null(ssbnp)) {
  nhanes_analytic <- nhanes_analytic %>%
    mutate(seqn = as.numeric(seqn)) %>%
    dplyr::left_join(ssbnp, by = "seqn", relationship = "one-to-one")
}

message("Analytic NHANES 1999-2004 dataset: ", nrow(nhanes_analytic), " rows x ",
        ncol(nhanes_analytic), " cols")

## ---- 4.5. Diagnostic: per-year data coverage audit --------------------
check_df <- nhanes_analytic %>% select(-any_of(c("sdmvpsu", "sdmvstra")), -seqn, -year)
years <- sort(unique(nhanes_analytic$year))

coverage <- sapply(years, function(yr) {
  round(100 * colMeans(!is.na(check_df[nhanes_analytic$year == yr, , drop = FALSE])), 1)
})
colnames(coverage) <- years

message("\n---- Per-year data coverage (% non-missing) ----")
print(coverage)

flagged <- coverage[apply(coverage == 0, 1, any), , drop = FALSE]
if (nrow(flagged) > 0) {
  message("\n*** FLAGGED: 0% present in at least one analytic year -- likely another naming/coverage quirk: ***")
  print(flagged)
} else {
  message("\nNo variable is fully missing (0%) within any single analytic year.")
}
message("Note: a variable can legitimately have PARTIAL missingness in a cycle (e.g. TRIGLY's ~50%",
        " fasting-subsample-only coverage) without being a bug -- only 0% in a cycle where the other",
        " two cycles have real data is a red flag worth checking against CDC's documentation.")

## ---- 5. Merge with mortality ----------------------------------------
## Option A: nhanesdata's built-in harmonized mortality linkage
mortality_pkg <- tryCatch(nhanesdata::read_nhanes("mortality"), error = function(e) NULL)

## Option B: your own Mortality_merged.xlsx built earlier from the raw
## NCHS .dat files -- use this as a cross-check, or as the primary
## source if you'd rather not depend on nhanesdata's mortality build.
mortality_own <- tryCatch(
  readxl::read_excel("data/Mortality_merged.xlsx"),
  error = function(e) NULL
)

mortality_source <- if (!is.null(mortality_pkg)) {
  message("Using nhanesdata::read_nhanes('mortality') for mortality linkage.")
  mortality_pkg %>% rename_with(tolower)
} else if (!is.null(mortality_own)) {
  message("nhanesdata mortality unavailable -- using your own Mortality_merged.xlsx instead.")
  mortality_own %>% rename_with(tolower)
} else {
  stop("No mortality source available. Check nhanesdata::read_nhanes('mortality') or the path to Mortality_merged.xlsx.")
}

nhanes_analytic <- nhanes_analytic %>% mutate(seqn = as.numeric(seqn))
mortality_source <- mortality_source %>% mutate(seqn = as.numeric(seqn))

hf_continuum_data <- dplyr::left_join(nhanes_analytic, mortality_source, by = "seqn")

message("Final HF-continuum analytic dataset: ", nrow(hf_continuum_data), " rows x ",
        ncol(hf_continuum_data), " cols")

## ---- 6. Save --------------------------------------------------------
dir.create("data", showWarnings = FALSE)
writexl::write_xlsx(hf_continuum_data, "data/HF_continuum_analytic_dataset.xlsx")
message("Saved: data/HF_continuum_analytic_dataset.xlsx")