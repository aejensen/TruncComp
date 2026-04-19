load_local_trunccomp2 <- function(repo_root) {
  local_pkg <- file.path(repo_root, "packages", "TruncComp2")
  use_local <- identical(tolower(Sys.getenv("TRUNCCOMP_USE_LOCAL_PACKAGE", "false")), "true")
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
