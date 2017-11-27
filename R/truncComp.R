truncComp <- function(y, a, z, method, conf.level = 0.95, init = NULL) {
  returnErrorData <- function(error) {
    out <- list(muDelta = NULL,
                muDeltaCI = NULL,
                alphaDelta = NULL,
                alphaDeltaCI = NULL,
                W = NULL,
                p = NULL,
                method = method,
                conf.level = conf.level,
                success = TRUE,
                error = error,
                init = NULL)
    class(out) <- append(class(out), "CompositeOutcomeEstimation")
    out
  }

  d <- data.frame(Y = y, A = a, Z = z)

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
      message("Calculating starting values for likelihood optimization")
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

  out
}

