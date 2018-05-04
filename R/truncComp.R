truncComp.formula <- function(formula, atom, data, method, conf.level = 0.95, init = NULL) {
  formulaVars <- all.vars(formula)
  d <- data.frame(Y = data[, formulaVars[1]], A = 1, R = data[, formulaVars[2]])
  d$A[d$Y == atom] <- 0

  truncComp.default(d$Y, d$A, d$R, method, conf.level, init)
}


truncComp.default <- function(y, a, r, method, conf.level = 0.95, init = NULL) {
  returnErrorData <- function(error) {
    out <- list(muDelta = NULL,
                muDeltaCI = NULL,
                alphaDelta = NULL,
                alphaDeltaCI = NULL,
                W = NULL,
                p = NULL,
                method = method,
                conf.level = conf.level,
                success = FALSE,
                error = error,
                init = NULL,
                data = NULL)
    class(out) <- append(class(out), "CompositeOutcomeEstimation")
    out
  }

  d <- data.frame(Y = y, A = a, R = r)

  if(!isDataOkay(d)) {
    error <- "Estimation failed due to data error."
    warning(error)
    out <- returnErrorData(error)
  }

  if(method != "LRT" & method != "SPLRT") {
    stop("method should be either LRT or SPLRT")
  }

  if(method == "LRT") {
    if(is.null(init)) {
      #message("Calculating starting values for likelihood optimization")
      init <- getLRstartingValues(d)
    }

    out <- LRT(d, init, conf.level)

    if(is.null(out)) {
      error <- "Estimation failed due to optimization error."
      warning(error)
      out <- returnErrorData(error)
    }
  } else if(method == "SPLRT") {
    #SPLRT cannot fail?
    out <- SPLRT(d, conf.level)
  }

  out$data <- d
  out
}

