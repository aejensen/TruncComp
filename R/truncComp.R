truncComp <- function(formula, atom, data, method, conf.level = 0.95, init = NULL) {
  if(!inherits(formula, "formula")) {
    stop("The formula must be a formula.")
  }

  if(!(method == "LRT" | method == "SPLRT")) {
    stop("Only LRT or SPLRT supported as methods.")
  }

  if(length(attr(terms(formula), "term.labels")) != 1) {
    stop("The current implementation must have one covariate in the formula.")
  }

  #Get the variables from the formula
  variables <- get_all_vars(formula, data = data)
  outcome <- variables[,1]
  treatment <- variables[,2]

  if(!all(sort(unique(treatment)) == c(0,1))) {
    stop("The covariate must be binary (0/1) indicating the two treatments.")
  }

  if(!(length(atom) == 1 & is.numeric(atom))) {
    stop("atom must be a single numeric value.")
  }

  alive <- as.numeric(outcome != atom)

  if(all(alive[treatment == 0] == 0) | all(alive[treatment == 1] == 0)) {
    stop("Nothing has been observed in one of the groups. Cannot do estimation.")
  }

  if(all(alive == 1)) {
    stop("Everything seems to have been observed. You should use a different method.")
  }


  truncComp.default(outcome, alive, treatment, method, conf.level, init)
}


truncComp.default <- function(y, a, r, method, conf.level = 0.95, init = NULL) {
  d <- data.frame(Y = y, A = a, R = r)

  if(!isDataOkay(d)) {
    error <- "Estimation failed due to data error."
    warning(error)
    out <- returnErrorData(error, method, conf.level)
  }

  if(method == "LRT") {
    if(is.null(init)) {
      init <- getLRstartingValues(d)
    }

    out <- LRT(d, init, conf.level)

    if(is.null(out)) {
      error <- "Estimation failed due to optimization error."
      warning(error)
      out <- returnErrorData(error)
    }
  } else if(method == "SPLRT") {
    #SPLRT cannot fail? Yes it can.
    out <- SPLRT(d, conf.level)
  }

  out$data <- d
  out
}

