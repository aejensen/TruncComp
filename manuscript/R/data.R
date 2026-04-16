load_local_trunccomp <- function(repo_root) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required to load the local TruncComp package.", call. = FALSE)
  }

  pkgload::load_all(
    file.path(repo_root, "packages", "TruncComp"),
    quiet = TRUE,
    export_all = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE
  )

  invisible(TRUE)
}

load_example_data <- function(repo_root) {
  data_env <- new.env(parent = emptyenv())
  load(file.path(repo_root, "packages", "TruncComp", "data", "TruncCompExample.RData"), envir = data_env)
  data_env$TruncCompExample
}

.application_data_dir <- function(manuscript_dir) {
  ensure_dir(file.path(manuscript_dir, "application-data"))
}

.application_dryad_standardized_path <- function(manuscript_dir) {
  file.path(.application_data_dir(manuscript_dir), "dryad-ards-standardized.csv")
}

.application_dryad_raw_path <- function(manuscript_dir) {
  data_dir <- .application_data_dir(manuscript_dir)
  candidates <- c("dryad-ards-data.xlsx", "data.xlsx")
  matches <- file.path(data_dir, candidates)
  existing <- matches[file.exists(matches)]

  if (!length(existing)) {
    return(NULL)
  }

  existing[[1]]
}

.application_dryad_zip_path <- function(manuscript_dir) {
  data_dir <- .application_data_dir(manuscript_dir)
  candidates <- c(
    "doi_10_5061_dryad_7d8k0__v20170929.zip",
    "dryad-ards-data.zip",
    "data.zip"
  )
  matches <- file.path(data_dir, candidates)
  existing <- matches[file.exists(matches)]

  if (length(existing)) {
    return(existing[[1]])
  }

  zips <- list.files(data_dir, pattern = "\\.zip$", full.names = TRUE)
  if (!length(zips)) {
    return(NULL)
  }

  zips[[1]]
}

.normalize_application_names <- function(x) {
  tolower(gsub("[^[:alnum:]]+", " ", x))
}

.find_application_column <- function(data, patterns, label) {
  normalized <- .normalize_application_names(names(data))

  for (pattern in patterns) {
    idx <- which(grepl(pattern, normalized))
    if (length(idx)) {
      return(names(data)[idx[[1]]])
    }
  }

  stop(
    sprintf("Unable to locate a %s column. Checked patterns: %s", label, paste(patterns, collapse = ", ")),
    call. = FALSE
  )
}

.normalize_binary_group <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }

  if (is.numeric(x)) {
    vals <- sort(unique(stats::na.omit(x)))
    if (identical(vals, c(0, 1))) {
      return(as.integer(x))
    }
    if (length(vals) == 2) {
      return(as.integer(x == max(vals)))
    }
  }

  x_chr <- tolower(trimws(as.character(x)))
  x_chr[x_chr %in% c("", "na", "nan")] <- NA_character_
  non_drug <- grepl("non", x_chr) & grepl("drug", x_chr)
  drug <- grepl("drug", x_chr) & !non_drug
  out <- ifelse(non_drug, 0L, ifelse(drug, 1L, NA_integer_))

  if (length(unique(stats::na.omit(out))) == 2) {
    return(out)
  }

  stop("Could not coerce the application group variable to a binary 0/1 indicator.", call. = FALSE)
}

.validate_application_data <- function(data, dataset_label) {
  if (!all(c("Y", "R") %in% names(data))) {
    stop(sprintf("%s must contain columns named Y and R.", dataset_label), call. = FALSE)
  }

  keep <- is.finite(data$Y) & !is.na(data$R)
  data <- data[keep, , drop = FALSE]

  if (!nrow(data)) {
    stop(sprintf("%s did not contain any complete observations.", dataset_label), call. = FALSE)
  }

  data$R <- as.integer(data$R)
  if (!identical(sort(unique(data$R)), c(0L, 1L))) {
    stop(sprintf("%s must contain exactly two groups coded 0/1.", dataset_label), call. = FALSE)
  }

  positive_counts <- tapply(data$Y > 0, data$R, sum)
  if (any(positive_counts < 2)) {
    stop(sprintf("%s must contain at least two positive outcomes in each group.", dataset_label), call. = FALSE)
  }

  data
}

.application_metadata_dryad <- function(data) {
  list(
    dataset_key = "dryad-ards",
    dataset_name = "the public Dryad ARDS dataset",
    dataset_description = "the public Dryad ARDS comparison between non-drug-associated and drug-associated acute respiratory distress syndrome",
    group_labels = c("Non-drug-associated ARDS", "Drug-associated ARDS"),
    outcome_label = "Ventilator-free days",
    outcome_short = "ventilator-free days",
    atom_meaning = "no ventilator-free days",
    observed_event_label = "having positive ventilator-free days",
    display_x_label = "Ventilator-free days",
    histogram_breaks = seq(-0.5, max(28, ceiling(max(data$Y, na.rm = TRUE))) + 0.5, by = 1),
    histogram_xlim = c(0, max(28, ceiling(max(data$Y, na.rm = TRUE)))),
    histogram_cap = NULL,
    surface_resolution = 60
  )
}

.application_metadata_randhealth <- function() {
  list(
    dataset_key = "randhealth-inpdol-0-vs-50",
    dataset_name = "the RAND Health Insurance Experiment",
    dataset_description = "the public RAND Health Insurance Experiment year-1 comparison between the 0% and 50% coinsurance plans",
    group_labels = c("0% coinsurance", "50% coinsurance"),
    outcome_label = "Year 1 inpatient expenditure (USD)",
    outcome_short = "year 1 inpatient expenditure",
    atom_meaning = "no inpatient expenditure",
    observed_event_label = "having positive inpatient expenditure",
    display_x_label = "Year 1 inpatient expenditure (USD, values above 3,000 capped for display)",
    histogram_breaks = seq(0, 3000, by = 250),
    histogram_xlim = c(0, 3000),
    histogram_cap = 3000,
    surface_resolution = 60
  )
}

.load_dryad_application_standardized <- function(path) {
  raw <- utils::read.csv(path, stringsAsFactors = FALSE)
  data <- .validate_application_data(raw[, c("Y", "R")], "The standardized Dryad ARDS CSV")

  list(
    data = data,
    metadata = .application_metadata_dryad(data)
  )
}

.load_dryad_application_raw <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop(
      "The readxl package is required to use the local Dryad workbook. Install readxl or provide dryad-ards-standardized.csv.",
      call. = FALSE
    )
  }

  sheets <- readxl::excel_sheets(path)
  last_error <- NULL

  for (sheet in sheets) {
    standardized <- tryCatch({
      raw <- as.data.frame(readxl::read_excel(path, sheet = sheet))
      outcome_col <- .find_application_column(
        raw,
        patterns = c("(^| )vfd($| )", "ventilator.*free", "free.*ventilator"),
        label = "ventilator-free days"
      )
      group_col <- .find_application_column(
        raw,
        patterns = c("(^| )dards($| )", "drug.*associated", "cause", "etiolog", "(^| )group($| )"),
        label = "drug-associated group"
      )
      .validate_application_data(
        data.frame(
          Y = as.numeric(raw[[outcome_col]]),
          R = .normalize_binary_group(raw[[group_col]])
        ),
        "The local Dryad ARDS workbook"
      )
    }, error = function(e) {
      last_error <<- e
      NULL
    })

    if (!is.null(standardized)) {
      return(list(data = standardized, metadata = .application_metadata_dryad(standardized)))
    }
  }

  stop(
    paste(
      "Unable to parse the local Dryad workbook.",
      "Rename the relevant columns to Y and R and save it as dryad-ards-standardized.csv if auto-detection fails.",
      if (!is.null(last_error)) conditionMessage(last_error) else NULL
    ),
    call. = FALSE
  )
}

load_local_dryad_application_data <- function(manuscript_dir) {
  standardized_path <- .application_dryad_standardized_path(manuscript_dir)
  if (file.exists(standardized_path)) {
    message("Using local standardized Dryad ARDS application data.")
    return(.load_dryad_application_standardized(standardized_path))
  }

  raw_path <- .application_dryad_raw_path(manuscript_dir)
  if (!is.null(raw_path) && file.exists(raw_path)) {
    message("Using local Dryad ARDS workbook for the manuscript application.")
    return(.load_dryad_application_raw(raw_path))
  }

  zip_path <- .application_dryad_zip_path(manuscript_dir)
  if (!is.null(zip_path) && file.exists(zip_path)) {
    message("Using local Dryad ARDS zip archive for the manuscript application.")
    extract_dir <- tempfile("dryad-ards-")
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    extracted <- utils::unzip(zip_path, files = "data.xlsx", exdir = extract_dir)
    if (!length(extracted) || !file.exists(extracted[[1]])) {
      stop("Unable to extract data.xlsx from the local Dryad zip archive.", call. = FALSE)
    }
    return(.load_dryad_application_raw(extracted[[1]]))
  }

  NULL
}

.load_randhealth_application_data <- function() {
  url <- "https://vincentarelbundock.github.io/Rdatasets/csv/camerondata/randhealth.csv"
  raw <- utils::read.csv(url)
  raw <- subset(raw, year == 1 & coins %in% c(0, 50), select = c("inpdol", "coins"))
  data <- .validate_application_data(
    data.frame(
      Y = raw$inpdol,
      R = as.integer(raw$coins == 50)
    ),
    "The RAND HIE fallback application data"
  )

  list(
    data = data,
    metadata = .application_metadata_randhealth()
  )
}

load_application_data <- function(manuscript_dir) {
  dryad <- load_local_dryad_application_data(manuscript_dir)
  if (!is.null(dryad)) {
    return(dryad)
  }

  message("Local Dryad application data not found; using the RAND HIE fallback dataset.")
  .load_randhealth_application_data()
}
