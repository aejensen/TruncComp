isDataOkay <- function(d) {
  dAlive <- subset(d, d$A == 1)

  yAlive1 <- dAlive$Y[dAlive$R == 0]
  yAlive2 <- dAlive$Y[dAlive$R == 1]

  if(length(yAlive1) < 2 | length(yAlive2) < 2) {
    return(FALSE)
  }

  TRUE
}

isValid <- function(truncCompObj) {
  truncCompObj$success
}

truncCompMethod <- function(method) {
  if(identical(method, "LRT") || identical(method, "Parametric Likelihood Ratio Test")) {
    return("Parametric Likelihood Ratio Test")
  }

  if(identical(method, "SPLRT") || identical(method, "Semi-empirical Likelihood Ratio Test")) {
    return("Semi-empirical Likelihood Ratio Test")
  }

  method
}

newTruncComp <- function(muDelta = NULL, muDeltaCI = NULL,
                         alphaDelta = NULL, alphaDeltaCI = NULL,
                         Delta = NULL, DeltaCI = NULL,
                         W = NULL, p = NULL,
                         method, conf.level, success,
                         error = "", init = NULL, data = NULL) {
  out <- list(muDelta = muDelta,
              muDeltaCI = muDeltaCI,
              alphaDelta = alphaDelta,
              alphaDeltaCI = alphaDeltaCI,
              Delta = Delta,
              DeltaCI = DeltaCI,
              W = W,
              p = p,
              method = truncCompMethod(method),
              conf.level = conf.level,
              success = success,
              error = error,
              init = init,
              data = data)
  class(out) <- c("TruncComp", "list")
  out
}

returnErrorData <- function(error, method, conf.level, init = NULL, data = NULL) {
  newTruncComp(method = method,
               conf.level = conf.level,
               success = FALSE,
               error = error,
               init = init,
               data = data)
}
