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
## dataset = nhanesdata dataset name (VERIFY with term_search()/get_url()
## if a read_nhanes() call below fails -- catalog names aren't all
## independently confirmed here)
## vars    = lowercase variable names to keep from that dataset
##
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
## The following are NOT in var_map -- nhanesdata has 1999-2004 coverage gaps
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
## nhanesdata has gaps for these older-format lab files (see note above), so
## instead of guessing legacy CDC file names (LAB10, L40_2, etc.), we look
## them up dynamically per variable using nhanesA::nhanesSearchVarName(),
## which returns the exact correct table name AND begin year for every table
## containing that variable -- no guessing, self-correcting regardless of
## naming era.
##
## TWO SEPARATE QUIRKS confirmed in L40_B (2001-2002) alone:
##  1. nhanesSearchVarName()'s live scrape of CDC's variable index silently
##     dropped this table for lbxscr/lbxsal even though it genuinely exists
##     (L40_B.htm on CDC's site) -- likely CDC site instability during the
##     2026 government shutdown rather than a real gap. Fixed via
##     manual_table_overrides below.
##  2. Even once loaded, L40_B names creatinine "LBDSCR", not "LBXSCR" --
##     breaking the usual LBX=measured/LBD..SI=derived pattern used in
##     every other cycle. Fixed via var_aliases below.
## Given a single file had two different quirks, there may be others we
## haven't hit yet. Safeguards: (1) var_aliases for the one confirmed so
## far, (2) if a table loads but NONE of a group's variables (or aliases)
## are found in it, we print that table's full column list so a new alias
## is immediately visible instead of silently missing, and (3) a loud
## warning if any analytic year is still missing at the end.
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
  ## Force numeric: these are continuous lab values, but nhanesA sometimes
  ## returns a column as a labelled factor in one cycle's table and as a
  ## plain double in another's -- bind_rows() can't combine those types when
  ## stacking cycles together later, so standardize here (safe for
  ## continuous values; factor levels for these are numeric-as-text anyway).
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
  
  ## Loud warning if any analytic year is still missing after overrides --
  ## surfaces new gaps instead of letting them pass silently.
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

## Fold lab_list into the same year+seqn-keyed collection used in section 4
person_level_list <- c(person_level_list, lab_list)

## ---- 2C. Insurance + MCQ age-at-diagnosis variables via nhanesA -------
## Also came back empty, for two different reasons:
##  - Health insurance: the questionnaire was completely redesigned in
##    2005-2006. For 1999-2004 "covered by health insurance" is named
##    HID010, not HIQ011 (CDC: "HIQ011 is comparable to HID010 in
##    2003-2004") -- it never existed as HIQ011 in our analytic window.
##    We pull HID010 and rename it to hiq011 so downstream code/codebook
##    stay consistent.
##  - MCQ age-at-diagnosis follow-ups: same nhanesdata 1999-2004 coverage
##    gap as the lab files (section 2B), PLUS ALL FIVE age-at-diagnosis
##    variables are named MCQ180B/C/D/E/F (not MCD180*) for 1999-2004 --
##    confirmed by testing: MCD180C/D returned "no tables found" for this
##    window. The MCD180* prefix was only introduced in later cycles.
##  - fsdhh: confirmed (via coverage audit) 0% in 1999 and 2001, ~95% in
##    2003. HHFDSEC is the right column name (confirmed directly:
##    names(nhanesA::nhanes("FSQ")) and names(nhanesA::nhanes("FSQ_B"))
##    both contain it) -- but nhanesSearchVarName()'s live scrape missed
##    both tables anyway, the same failure mode as biopro/2001 (L40_B).
##    Fixed with a direct manual_table_overrides_legacy entry below,
##    bypassing the unreliable search entirely for this variable.
##  - kiq022: confirmed 0% in 1999, ~49% in 2001/2003 (partial coverage in
##    those years is expected -- likely a skip pattern, not a bug).
##    1999-2000 names this KIQ020, not KIQ022.
## fsdhh/kiq022 need a DIFFERENT source name in different cycles (unlike
## hiq011/mcd180*, which use one alternate name across all three), so each
## entry below is a vector of candidate names tried across the analytic
## window; whichever candidate exists in a given cycle's table is used.
##
## NOTE: self-rated health (protocol §7) is sourced as huq010, added
## directly to the "huq" entry in var_map above (not routed through this
## legacy_vars/nhanesA mechanism) -- it's in the same HUQ file as huq030,
## which already loads cleanly via nhanesdata with no known issues, so no
## special handling is needed. (An earlier version of this script tried
## the MEC-based HSD010/HSQ file instead; switched to HUQ010 per
## confirmation that it's present and consistently named across all three
## 1999-2004 cycles, and it avoids a whole separate file/loading path.)
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

## Confirmed table names for cases where nhanesSearchVarName()'s live scrape
## misses a table that genuinely exists (verified directly via
## names(nhanesA::nhanes(tbl))). Bypasses the search entirely for these.
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
    ## Coerce to character: nhanesA::nhanes() sometimes returns a column as a
    ## labelled factor in one cycle's table and as a plain numeric/character
    ## in another (confirmed: FSQ_C's FSDHH came back as a factor, FSQ's
    ## HHFDSEC as a double) -- bind_rows() can't combine those types later,
    ## so standardize here. Recode/convert as needed during analysis.
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
## CONFIRMED (via user testing): nhanes("SSBNP_A") returns the FULL
## pooled 1999-2004 NT-proBNP dataset in one call -- CDC released this
## one-time retrospective (archived-serum) study as a single combined
## file rather than per-cycle files, unlike standard NHANES components.
## There is no "SSBNP_B"/"SSBNP_C" -- SSBNP_A already has everyone.
##
## CONFIRMED per NHANES documentation (SSBNP_A/SSCARD_A/SSTROP_A share
## the same subsample): this file also carries the subsample-specific
## survey weights WTSSCB2Y and WTSSCB4Y, which MUST be used for this
## biomarker subsample instead of the general WTMEC2YR exam weight --
## using WTMEC2YR here was flagged in review as methodologically
## incorrect for this dataset. Each person has exactly ONE of the two
## weight columns populated, never both: WTSSCB4Y for participants
## examined during 1999-2002 (this subsample's 4-year weight already
## pools the 1999-2000 and 2001-2002 cycles together), WTSSCB2Y for
## participants examined during 2003-2004 alone. The correct 6-year
## pooled weight is constructed downstream (in the main analysis script)
## as WTSSCB4Y*(2/3) for the first group and WTSSCB2Y*(1/3) for the
## second, not a uniform /3 of a single weight variable.
##
## This file has NO year column. That's fine: NHANES SEQN values never
## repeat across cycles (each 2-year cycle uses a distinct, non-
## overlapping SEQN range), so merging by seqn alone is safe here --
## no risk of matching the wrong person from a different cycle.
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
  ## Defensive check: confirm the weight columns actually came through --
  ## if nhanesA ever changes what it returns for this file, this stops
  ## the pipeline here with a clear message instead of silently
  ## producing all-NA weights downstream.
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

## NT-proBNP joins by seqn only (pooled file, no year column -- see
## section 3 for why this is safe).
if (!is.null(ssbnp)) {
  nhanes_analytic <- nhanes_analytic %>%
    mutate(seqn = as.numeric(seqn)) %>%
    dplyr::left_join(ssbnp, by = "seqn", relationship = "one-to-one")
}

message("Analytic NHANES 1999-2004 dataset: ", nrow(nhanes_analytic), " rows x ",
        ncol(nhanes_analytic), " cols")

## ---- 4.5. Diagnostic: per-year data coverage audit --------------------
## The loading-time checks above (sections 2/2B/2C) don't catch everything:
##  - A multi-variable group only warns if it matches NONE of its target
##    variables in a table, not if it's missing just one (exactly what
##    LBDSCR was, before we added that alias).
##  - Variables still sourced via nhanesdata (section 2, e.g. mcq160b-f,
##    bpq*, diq010, smq*, paq*, bmx*, bpx*, huq030, demographics) never
##    got the per-year presence check built for the nhanesA-routed
##    variables -- they could still have an undiscovered 1999-2004 gap.
## This audits the ACTUAL final merged data instead of anticipating
## failure modes in advance, so it catches anything the checks above
## missed. Read this before trusting the dataset.
##
## Simple base-R version: split by year, compute % non-missing per column,
## bind into one table (rows = variable, columns = year).
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