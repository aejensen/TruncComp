load_local_trunccomp2 <- function(repo_root) {
  local_pkg <- file.path(repo_root, "packages", "TruncComp2")
  use_local <- !identical(tolower(Sys.getenv("TRUNCCOMP_USE_LOCAL_PACKAGE", "true")), "false")
  if (use_local && dir.exists(local_pkg) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(
      local_pkg,
      quiet = TRUE,
      export_all = FALSE,
      helpers = FALSE,
      attach_testthat = FALSE
    )
    return(invisible(TRUE))
  }

  if (!requireNamespace("TruncComp2", quietly = TRUE)) {
    stop(
      paste(
        "The TruncComp2 package is required to build the manuscript application.",
        "Install it or make packages/TruncComp2 available."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.application_data_dir <- function(manuscript_dir) {
  ensure_dir(file.path(manuscript_dir, "application-data"))
}

ist3_data_urls <- function() {
  list(
    record = "https://datashare.ed.ac.uk/handle/10283/1931",
    doi = "https://doi.org/10.7488/ds/1350",
    data = "https://datashare.ed.ac.uk/bitstream/handle/10283/1931/ist3.dat?sequence=15&isAllowed=y",
    sas_syntax = "https://datashare.ed.ac.uk/bitstream/handle/10283/1931/ist3.sas?sequence=16&isAllowed=y",
    sas_formats = "https://datashare.ed.ac.uk/bitstream/handle/10283/1931/sas_fmts.sas?sequence=17&isAllowed=y",
    dua = "https://datashare.ed.ac.uk/bitstream/handle/10283/1931/IST-3_Data_Use_Agreement%20v2.1%20.pdf?sequence=21&isAllowed=y"
  )
}

.application_ist3_standardized_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "ist3-standardized.csv")
}

.application_ist3_raw_candidates <- function(manuscript_dir) {
  c(
    file.path(.application_data_dir(manuscript_dir), "ist3.dat"),
    "/tmp/trunccomp-dataset-hunt/ist3/ist3.dat"
  )
}

.application_ist3_raw_path <- function(manuscript_dir) {
  candidates <- .application_ist3_raw_candidates(manuscript_dir)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) {
    return(existing[[1]])
  }

  NULL
}

download_ist3_application_data <- function(manuscript_dir) {
  data_dir <- .application_data_dir(manuscript_dir)
  dest <- file.path(data_dir, "ist3.dat")

  if (file.exists(dest)) {
    return(dest)
  }

  urls <- ist3_data_urls()
  message("Downloading IST-3 fixed-width data from Edinburgh DataShare.")
  utils::download.file(urls$data, destfile = dest, mode = "wb", quiet = TRUE)
  if (!file.exists(dest) || file.size(dest) == 0) {
    stop("The IST-3 download did not create a non-empty ist3.dat file.", call. = FALSE)
  }

  dest
}

.parse_ist3_fixed_width <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("IST-3 fixed-width file not found: %s", path), call. = FALSE)
  }

  lines <- readLines(path, warn = FALSE)
  read_num <- function(first, last) {
    suppressWarnings(as.numeric(trimws(substr(lines, first, last))))
  }
  read_chr <- function(first, last) {
    trimws(substr(lines, first, last))
  }

  data.frame(
    row_id = seq_along(lines),
    itt_treat = read_num(27, 27),
    dead6mo = read_num(309, 309),
    euroqol6 = read_num(341, 343),
    treatment = read_chr(503, 509),
    stringsAsFactors = FALSE
  )
}

.validate_ist3_coding <- function(raw) {
  coding <- stats::xtabs(~ itt_treat + treatment, data = raw)
  if (!("0" %in% rownames(coding)) || !("1" %in% rownames(coding))) {
    stop("The IST-3 allocated-treatment variable did not contain both documented codes 0 and 1.", call. = FALSE)
  }
  if (!("rt-PA" %in% colnames(coding)) || !("Placebo" %in% colnames(coding))) {
    stop("The IST-3 treatment text field did not contain both documented labels rt-PA and Placebo.", call. = FALSE)
  }
  if (coding["0", "rt-PA"] == 0 || coding["1", "Placebo"] == 0) {
    stop("The IST-3 treatment coding check failed: expected itt_treat 0 = rt-PA and 1 = Placebo.", call. = FALSE)
  }
  if (any(coding["0", setdiff(colnames(coding), "rt-PA")] > 0) ||
      any(coding["1", setdiff(colnames(coding), "Placebo")] > 0)) {
    stop("The IST-3 treatment coding is not one-to-one with the treatment text field.", call. = FALSE)
  }

  invisible(coding)
}

.ist3_analysis_frame <- function(raw) {
  .validate_ist3_coding(raw)

  raw$R <- ifelse(raw$itt_treat == 0, 1L, ifelse(raw$itt_treat == 1, 0L, NA_integer_))
  raw$group <- ifelse(raw$R == 1L, "rt-PA", ifelse(raw$R == 0L, "Placebo/control", NA_character_))
  raw$is_dead6mo <- raw$dead6mo == 1
  raw$is_alive6mo <- raw$dead6mo == 0
  raw$is_alive_missing_euroqol6 <- raw$is_alive6mo & is.na(raw$euroqol6)
  raw$included <- !is.na(raw$R) & (raw$is_dead6mo | (raw$is_alive6mo & is.finite(raw$euroqol6)))

  analysis <- raw[raw$included, , drop = FALSE]
  analysis$Y <- ifelse(analysis$is_dead6mo, -1, analysis$euroqol6)
  analysis$A <- as.integer(analysis$Y != -1)
  analysis$atom <- -1
  analysis <- analysis[, c(
    "row_id", "Y", "R", "A", "atom", "group", "itt_treat", "treatment",
    "dead6mo", "euroqol6", "is_dead6mo"
  )]

  if (!identical(sort(unique(analysis$R)), c(0L, 1L))) {
    stop("The IST-3 analysis data must contain treatment groups coded 0/1.", call. = FALSE)
  }
  if (any(tapply(analysis$A == 1L, analysis$R, sum) < 20)) {
    stop("The IST-3 analysis data do not contain enough non-atom observations in both groups.", call. = FALSE)
  }
  if (any(!is.finite(analysis$Y))) {
    stop("The IST-3 analysis data contain missing combined outcomes.", call. = FALSE)
  }

  list(raw = raw, data = analysis)
}

.application_metadata_ist3 <- function(raw, data, raw_path, standardized_path) {
  urls <- ist3_data_urls()
  randomized_counts <- as.integer(table(factor(raw$R, levels = 0:1)))
  included_counts <- as.integer(table(factor(data$R, levels = 0:1)))
  excluded_alive_missing <- as.integer(tapply(raw$is_alive_missing_euroqol6, factor(raw$R, levels = 0:1), sum, na.rm = TRUE))

  list(
    dataset_key = "ist3",
    dataset_name = "IST-3 stroke trial",
    dataset_description = "the third International Stroke Trial (IST-3), a randomized trial of intravenous alteplase (rt-PA) versus control in acute ischemic stroke",
    group_labels = c("Placebo/control", "rt-PA"),
    outcome_label = "death by 6 months plus 6-month EuroQol visual analogue score among survivors",
    outcome_short = "6-month EQ-VAS with death coded as -1",
    atom = -1,
    atom_label = "death by 6 months",
    atom_meaning = "death before the 6-month EuroQol assessment",
    observed_event_label = "being alive with a recorded EuroQol visual analogue score at 6 months",
    continuous_label = "EuroQol visual analogue score among 6-month survivors",
    display_x_label = "Six-month endpoint: death = -1; survivor EQ-VAS = 0 to 100",
    surface_resolution = 55,
    raw_path = raw_path,
    standardized_path = standardized_path,
    source_record_url = urls$record,
    source_data_url = urls$data,
    source_sas_url = urls$sas_syntax,
    source_dua_url = urls$dua,
    doi = "10.7488/ds/1350",
    citation_text = "Sandercock, Wardlaw, Lindley, Cohen, and Whiteley (2016)",
    randomized_n = nrow(raw),
    randomized_group_n = randomized_counts,
    analysis_group_n = included_counts,
    excluded_alive_missing_euroqol6 = sum(raw$is_alive_missing_euroqol6, na.rm = TRUE),
    excluded_alive_missing_euroqol6_by_group = excluded_alive_missing,
    dua_note = "The dataset is publicly downloadable from Edinburgh DataShare after the embargo period, but reuse is conditional on the IST-3 Data Use Agreement."
  )
}

load_ist3_application_data <- function(manuscript_dir) {
  standardized_path <- .application_ist3_standardized_path(manuscript_dir)
  raw_path <- .application_ist3_raw_path(manuscript_dir)
  if (is.null(raw_path)) {
    raw_path <- download_ist3_application_data(manuscript_dir)
  }
  parsed <- .ist3_analysis_frame(.parse_ist3_fixed_width(raw_path))

  if (file.exists(standardized_path)) {
    standardized <- utils::read.csv(standardized_path, stringsAsFactors = FALSE)
    if (!all(c("Y", "R", "A", "atom") %in% names(standardized))) {
      stop("The standardized IST-3 CSV must contain Y, R, A, and atom.", call. = FALSE)
    }
    check_cols <- c("Y", "R", "A", "atom")
    same_standardized <- nrow(standardized) == nrow(parsed$data) &&
      isTRUE(all.equal(standardized[, check_cols], parsed$data[, check_cols], check.attributes = FALSE))
    if (!same_standardized) {
      message("Refreshing stale IST-3 standardized CSV from the raw fixed-width file.")
      utils::write.csv(parsed$data, standardized_path, row.names = FALSE)
      standardized <- parsed$data
    }

    return(list(
      data = standardized,
      raw = parsed$raw,
      metadata = .application_metadata_ist3(parsed$raw, standardized, raw_path, standardized_path)
    ))
  }

  utils::write.csv(parsed$data, standardized_path, row.names = FALSE)

  list(
    data = parsed$data,
    raw = parsed$raw,
    metadata = .application_metadata_ist3(parsed$raw, parsed$data, raw_path, standardized_path)
  )
}

load_application_data <- function(manuscript_dir) {
  load_ist3_application_data(manuscript_dir)
}

.application_liver_standardized_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "liver-appendix-standardized.csv")
}

.load_joiner_liver_raw <- function() {
  if (!requireNamespace("joineR", quietly = TRUE)) {
    stop(
      paste(
        "The joineR package is required to build the liver appendix application.",
        "Install joineR from CRAN to load joineR::liver."
      ),
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())
  utils::data("liver", package = "joineR", envir = env)
  raw <- get("liver", envir = env)
  raw <- as.data.frame(raw)
  raw
}

.validate_joiner_liver <- function(raw) {
  required <- c("id", "prothrombin", "time", "treatment", "survival", "cens")
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop(
      sprintf("The joineR liver data are missing required columns: %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!identical(sort(unique(raw$treatment)), c(0L, 1L))) {
    stop("The joineR liver treatment variable must contain documented codes 0 = placebo and 1 = prednisone.", call. = FALSE)
  }
  if (!all(unique(raw$cens) %in% c(0L, 1L))) {
    stop("The joineR liver censoring indicator must use documented codes 0 = censored and 1 = died.", call. = FALSE)
  }

  subject_checks <- lapply(split(raw, raw$id), function(x) {
    c(
      treatment = length(unique(x$treatment)),
      survival = length(unique(x$survival)),
      cens = length(unique(x$cens))
    )
  })
  subject_checks <- do.call(rbind, subject_checks)
  if (any(subject_checks != 1L)) {
    stop("The joineR liver data have subject-level variables that vary within patient id.", call. = FALSE)
  }

  if (any(!is.finite(raw$prothrombin)) || any(!is.finite(raw$time))) {
    stop("The joineR liver prothrombin and time columns must be finite for all longitudinal rows.", call. = FALSE)
  }

  invisible(TRUE)
}

.liver_appendix_subject_frame <- function(raw, landmark = 2) {
  .validate_joiner_liver(raw)

  rows <- lapply(split(raw, raw$id), function(x) {
    x <- x[order(x$time), , drop = FALSE]
    id <- x$id[[1]]
    tr <- x$treatment[[1]]
    surv <- x$survival[[1]]
    cens <- x$cens[[1]]
    died_by_landmark <- cens == 1L && surv <= landmark
    censored_before_landmark <- cens == 0L && surv < landmark
    prior <- x[x$time <= landmark, , drop = FALSE]

    if (died_by_landmark) {
      return(data.frame(
        id = id,
        R = tr,
        group = ifelse(tr == 1L, "Prednisone", "Placebo"),
        Y = 0,
        A = 0L,
        atom = 0,
        prothrombin = NA_real_,
        measurement_time = NA_real_,
        measurement_lag = NA_real_,
        survival = surv,
        cens = cens,
        included = TRUE,
        reason = "death before two-year landmark",
        stringsAsFactors = FALSE
      ))
    }

    if (censored_before_landmark) {
      return(data.frame(
        id = id,
        R = tr,
        group = ifelse(tr == 1L, "Prednisone", "Placebo"),
        Y = NA_real_,
        A = NA_integer_,
        atom = 0,
        prothrombin = NA_real_,
        measurement_time = NA_real_,
        measurement_lag = NA_real_,
        survival = surv,
        cens = cens,
        included = FALSE,
        reason = "censored before two-year landmark",
        stringsAsFactors = FALSE
      ))
    }

    if (!nrow(prior)) {
      return(data.frame(
        id = id,
        R = tr,
        group = ifelse(tr == 1L, "Prednisone", "Placebo"),
        Y = NA_real_,
        A = NA_integer_,
        atom = 0,
        prothrombin = NA_real_,
        measurement_time = NA_real_,
        measurement_lag = NA_real_,
        survival = surv,
        cens = cens,
        included = FALSE,
        reason = "no prothrombin before two-year landmark",
        stringsAsFactors = FALSE
      ))
    }

    idx <- which.max(prior$time)
    y <- prior$prothrombin[[idx]]
    measurement_time <- prior$time[[idx]]
    data.frame(
      id = id,
      R = tr,
      group = ifelse(tr == 1L, "Prednisone", "Placebo"),
      Y = y,
      A = 1L,
      atom = 0,
      prothrombin = y,
      measurement_time = measurement_time,
      measurement_lag = landmark - measurement_time,
      survival = surv,
      cens = cens,
      included = TRUE,
      reason = "known alive at two-year landmark",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.liver_appendix_analysis_frame <- function(raw, landmark = 2) {
  subjects <- .liver_appendix_subject_frame(raw, landmark = landmark)
  analysis <- subjects[subjects$included, , drop = FALSE]
  analysis <- analysis[, c(
    "id", "Y", "R", "A", "atom", "group", "prothrombin",
    "measurement_time", "measurement_lag", "survival", "cens", "reason"
  )]

  if (!identical(sort(unique(analysis$R)), c(0L, 1L))) {
    stop("The joineR liver analysis data must contain treatment groups coded 0/1.", call. = FALSE)
  }
  if (any(tapply(analysis$A == 1L, analysis$R, sum) < 20)) {
    stop("The joineR liver analysis data do not contain enough non-atom observations in both groups.", call. = FALSE)
  }
  if (any(!is.finite(analysis$Y))) {
    stop("The joineR liver analysis data contain missing combined outcomes.", call. = FALSE)
  }
  if (any(analysis$A == 1L & analysis$Y <= 0)) {
    stop("The joineR liver prothrombin component must be strictly positive away from the atom.", call. = FALSE)
  }

  list(raw = raw, subjects = subjects, data = analysis)
}

.liver_reason_counts <- function(subjects, reason) {
  as.integer(tapply(subjects$reason == reason, factor(subjects$R, levels = 0:1), sum, na.rm = TRUE))
}

.application_metadata_liver <- function(raw, subjects, data, standardized_path, landmark = 2) {
  randomized_counts <- as.integer(table(factor(subjects$R, levels = 0:1)))
  included_counts <- as.integer(table(factor(data$R, levels = 0:1)))
  excluded_censored <- .liver_reason_counts(subjects, "censored before two-year landmark")
  excluded_no_prior <- .liver_reason_counts(subjects, "no prothrombin before two-year landmark")

  list(
    dataset_key = "joiner-liver",
    dataset_name = "joineR liver cirrhosis trial",
    dataset_description = "the liver cirrhosis prednisone trial data distributed with the joineR R package",
    group_labels = c("Placebo", "Prednisone"),
    outcome_label = "death by 2 years plus latest prothrombin index at or before 2 years among known 2-year survivors",
    outcome_short = "2-year death/prothrombin endpoint",
    atom = 0,
    atom_label = "death by 2 years",
    atom_meaning = "death before or at the two-year landmark",
    observed_event_label = "being known alive at two years with a recorded prothrombin index at or before the landmark",
    continuous_label = "latest prothrombin index (%) at or before 2 years among known 2-year survivors",
    display_x_label = "Two-year endpoint: death = 0; survivor prothrombin index (%)",
    surface_resolution = 55,
    standardized_path = standardized_path,
    source_package = "joineR",
    source_data_object = "liver",
    source_package_url = "https://github.com/graemeleehickey/joineR/",
    package_version = as.character(utils::packageVersion("joineR")),
    citation_text = "Philipson et al. (2018)",
    landmark = landmark,
    randomized_n = length(unique(raw$id)),
    randomized_group_n = randomized_counts,
    analysis_group_n = included_counts,
    excluded_censored_before_landmark = sum(excluded_censored),
    excluded_censored_before_landmark_by_group = excluded_censored,
    excluded_no_prior_measurement = sum(excluded_no_prior),
    excluded_no_prior_measurement_by_group = excluded_no_prior,
    availability_note = "The data are distributed with the CRAN joineR package and can be loaded with data(liver, package = \"joineR\")."
  )
}

load_liver_appendix_data <- function(manuscript_dir, landmark = 2) {
  standardized_path <- .application_liver_standardized_path(manuscript_dir)
  parsed <- .liver_appendix_analysis_frame(.load_joiner_liver_raw(), landmark = landmark)

  if (file.exists(standardized_path)) {
    standardized <- utils::read.csv(standardized_path, stringsAsFactors = FALSE)
    if (!all(c("Y", "R", "A", "atom") %in% names(standardized))) {
      stop("The standardized joineR liver CSV must contain Y, R, A, and atom.", call. = FALSE)
    }
    check_cols <- c("Y", "R", "A", "atom", "prothrombin", "measurement_time", "measurement_lag")
    same_standardized <- nrow(standardized) == nrow(parsed$data) &&
      isTRUE(all.equal(standardized[, check_cols], parsed$data[, check_cols], check.attributes = FALSE))
    if (!same_standardized) {
      message("Refreshing stale joineR liver standardized CSV from the package data.")
      utils::write.csv(parsed$data, standardized_path, row.names = FALSE)
      standardized <- parsed$data
    }

    return(list(
      data = standardized,
      raw = parsed$raw,
      subjects = parsed$subjects,
      metadata = .application_metadata_liver(parsed$raw, parsed$subjects, standardized, standardized_path, landmark = landmark)
    ))
  }

  utils::write.csv(parsed$data, standardized_path, row.names = FALSE)

  list(
    data = parsed$data,
    raw = parsed$raw,
    subjects = parsed$subjects,
    metadata = .application_metadata_liver(parsed$raw, parsed$subjects, parsed$data, standardized_path, landmark = landmark)
  )
}

.application_licorice_standardized_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "licorice-appendix-standardized.csv")
}

.load_medicaldata_licorice_raw <- function() {
  if (!requireNamespace("medicaldata", quietly = TRUE)) {
    stop(
      paste(
        "The medicaldata package is required to build the licorice appendix application.",
        "Install medicaldata from CRAN to load medicaldata::licorice_gargle."
      ),
      call. = FALSE
    )
  }

  env <- new.env(parent = emptyenv())
  utils::data("licorice_gargle", package = "medicaldata", envir = env)
  raw <- get("licorice_gargle", envir = env)
  raw <- as.data.frame(raw)
  raw$row_id <- seq_len(nrow(raw))
  raw
}

.validate_medicaldata_licorice <- function(raw, outcome_var = "pacu30min_swallowPain") {
  required <- c("row_id", "treat", outcome_var)
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop(
      sprintf("The medicaldata licorice_gargle data are missing required columns: %s.", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!setequal(sort(unique(raw$treat)), c(0, 1))) {
    stop("The licorice gargle treatment variable must contain documented codes 0 = sugar-water and 1 = licorice.", call. = FALSE)
  }
  observed <- raw[[outcome_var]][!is.na(raw[[outcome_var]])]
  if (!all(is.finite(observed)) || any(observed < 0) || any(observed > 10)) {
    stop("The selected licorice gargle pain endpoint must be finite on the documented 0 to 10 scale.", call. = FALSE)
  }
  if (length(setdiff(unique(observed), 0:10))) {
    stop("The selected licorice gargle pain endpoint contains undocumented score values outside the 0 to 10 integer scale.", call. = FALSE)
  }

  invisible(TRUE)
}

.licorice_appendix_analysis_frame <- function(raw, outcome_var = "pacu30min_swallowPain") {
  .validate_medicaldata_licorice(raw, outcome_var = outcome_var)

  raw$R <- as.integer(raw$treat)
  raw$group <- ifelse(raw$R == 1L, "Licorice", "Sugar-water")
  raw$outcome_value <- raw[[outcome_var]]
  raw$included <- !is.na(raw$R) & is.finite(raw$outcome_value)
  raw$reason <- ifelse(raw$included, "observed 30-minute swallowing-pain score", "missing 30-minute swallowing-pain score")

  analysis <- raw[raw$included, , drop = FALSE]
  analysis$Y <- analysis$outcome_value
  analysis$A <- as.integer(analysis$Y != 0)
  analysis$atom <- 0
  analysis <- analysis[, c(
    "row_id", "Y", "R", "A", "atom", "group", "treat",
    "outcome_value", "pacu30min_swallowPain", "pacu30min_throatPain",
    "pacu90min_throatPain", "postOp4hour_throatPain", "pod1am_throatPain",
    "reason"
  )]

  if (!identical(sort(unique(analysis$R)), c(0L, 1L))) {
    stop("The licorice gargle analysis data must contain treatment groups coded 0/1.", call. = FALSE)
  }
  if (any(tapply(analysis$A == 1L, analysis$R, sum) < 20)) {
    stop("The licorice gargle analysis data do not contain enough positive-score observations in both groups.", call. = FALSE)
  }
  if (any(!is.finite(analysis$Y))) {
    stop("The licorice gargle analysis data contain missing combined outcomes.", call. = FALSE)
  }
  if (any(analysis$A == 1L & analysis$Y <= 0)) {
    stop("The licorice gargle positive pain component must be strictly positive away from the atom.", call. = FALSE)
  }

  list(raw = raw, data = analysis)
}

.licorice_missing_outcome_counts <- function(raw, outcome_var = "pacu30min_swallowPain") {
  as.integer(tapply(is.na(raw[[outcome_var]]), factor(raw$treat, levels = 0:1), sum, na.rm = TRUE))
}

.application_metadata_licorice <- function(raw, data, standardized_path, outcome_var = "pacu30min_swallowPain") {
  randomized_counts <- as.integer(table(factor(raw$treat, levels = 0:1)))
  included_counts <- as.integer(table(factor(data$R, levels = 0:1)))
  missing_counts <- .licorice_missing_outcome_counts(raw, outcome_var = outcome_var)

  list(
    dataset_key = "licorice-gargle",
    dataset_name = "licorice gargle randomized trial",
    dataset_description = "the licorice versus sugar-water gargle randomized clinical trial data distributed with the medicaldata R package",
    group_labels = c("Sugar-water", "Licorice"),
    outcome_label = "30-minute PACU swallowing-pain score with no pain as the atom",
    outcome_short = "30-minute swallowing-pain endpoint",
    atom = 0,
    atom_label = "no swallowing pain at 30 minutes",
    atom_meaning = "no sore-throat pain during swallowing at 30 minutes after PACU arrival",
    observed_event_label = "having any sore-throat pain during swallowing at 30 minutes after PACU arrival",
    continuous_label = "positive 30-minute swallowing-pain score on an 11-point Likert scale",
    display_x_label = "Thirty-minute swallowing-pain endpoint: no pain = 0; pain score = 1 to 10",
    surface_resolution = 55,
    continuous_support = "positive_real",
    standardized_path = standardized_path,
    source_package = "medicaldata",
    source_data_object = "licorice_gargle",
    source_package_url = "https://higgi13425.github.io/medicaldata/",
    source_repository_url = "https://github.com/higgi13425/medicaldata/",
    source_portal_url = "https://www.causeweb.org/tshs/licorice-gargle/",
    source_paper_url = "https://doi.org/10.1213/ANE.0b013e318299a650",
    package_version = as.character(utils::packageVersion("medicaldata")),
    package_license = utils::packageDescription("medicaldata")$License,
    citation_text = "Nowacki (2017), Ruetzler et al. (2013), and Higgins (2021)",
    outcome_var = outcome_var,
    randomized_n = nrow(raw),
    randomized_group_n = randomized_counts,
    analysis_group_n = included_counts,
    excluded_missing_outcome = sum(missing_counts),
    excluded_missing_outcome_by_group = missing_counts,
    availability_note = paste(
      "The analysis data are distributed with the CRAN medicaldata package and can be loaded with data(licorice_gargle, package = \"medicaldata\").",
      "The original TSHS portal is educational and states that publication reuse permission is not granted by the portal itself."
    )
  )
}

load_licorice_appendix_data <- function(manuscript_dir, outcome_var = "pacu30min_swallowPain") {
  standardized_path <- .application_licorice_standardized_path(manuscript_dir)
  parsed <- .licorice_appendix_analysis_frame(.load_medicaldata_licorice_raw(), outcome_var = outcome_var)

  if (file.exists(standardized_path)) {
    standardized <- utils::read.csv(standardized_path, stringsAsFactors = FALSE)
    if (!all(c("Y", "R", "A", "atom") %in% names(standardized))) {
      stop("The standardized licorice gargle CSV must contain Y, R, A, and atom.", call. = FALSE)
    }
    check_cols <- c("Y", "R", "A", "atom", "outcome_value")
    same_standardized <- nrow(standardized) == nrow(parsed$data) &&
      isTRUE(all.equal(standardized[, check_cols], parsed$data[, check_cols], check.attributes = FALSE))
    if (!same_standardized) {
      message("Refreshing stale licorice gargle standardized CSV from the package data.")
      utils::write.csv(parsed$data, standardized_path, row.names = FALSE)
      standardized <- parsed$data
    }

    return(list(
      data = standardized,
      raw = parsed$raw,
      metadata = .application_metadata_licorice(parsed$raw, standardized, standardized_path, outcome_var = outcome_var)
    ))
  }

  utils::write.csv(parsed$data, standardized_path, row.names = FALSE)

  list(
    data = parsed$data,
    raw = parsed$raw,
    metadata = .application_metadata_licorice(parsed$raw, parsed$data, standardized_path, outcome_var = outcome_var)
  )
}
