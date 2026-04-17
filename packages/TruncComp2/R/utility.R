isDataOkay <- function(d) {
  dAlive <- subset(d, d$A == 1)

  yAlive1 <- dAlive$Y[dAlive$R == 0]
  yAlive2 <- dAlive$Y[dAlive$R == 1]

  if(length(yAlive1) < 2 || length(yAlive2) < 2) {
    return(FALSE)
  }

  TRUE
}

adjustmentSpecification <- function(adjust) {
  if(is.null(adjust)) {
    return(NULL)
  }

  labels <- attr(stats::terms(adjust), "term.labels")
  if(length(labels) == 0) {
    return(NULL)
  }

  paste(labels, collapse = " + ")
}

isValid <- function(truncCompObj) {
  isTRUE(truncCompObj$success)
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

syncTruncComp2Aliases <- function(object) {
  object$mu_delta <- object$muDelta <- object$mu_delta
  object$mu_delta_ci <- object$muDeltaCI <- object$mu_delta_ci
  object$alpha_delta <- object$alphaDelta <- object$alpha_delta
  object$alpha_delta_ci <- object$alphaDeltaCI <- object$alpha_delta_ci
  object$delta <- object$Delta <- object$delta
  object$statistic <- object$W <- object$statistic
  object$conf_level <- object$conf.level <- object$conf_level
  object
}

newTruncComp2 <- function(muDelta = NULL, muDeltaCI = NULL,
                         alphaDelta = NULL, alphaDeltaCI = NULL,
                         Delta = NULL,
                         W = NULL, p = NULL,
                         method, conf.level, success,
                         error = "", init = NULL, data = NULL,
                         adjust = NULL, atom = NULL) {
  method <- truncCompMethod(method)

  out <- list(
    mu_delta = muDelta,
    mu_delta_ci = muDeltaCI,
    alpha_delta = alphaDelta,
    alpha_delta_ci = alphaDeltaCI,
    delta = Delta,
    statistic = W,
    p = p,
    method = method,
    conf_level = conf.level,
    success = success,
    error = error,
    init = init,
    data = data,
    adjust = adjust,
    atom = atom,
    muDelta = muDelta,
    muDeltaCI = muDeltaCI,
    alphaDelta = alphaDelta,
    alphaDeltaCI = alphaDeltaCI,
    Delta = Delta,
    W = W,
    conf.level = conf.level
  )
  class(out) <- c("trunc_comp_fit", "TruncComp2", "list")
  syncTruncComp2Aliases(out)
}

returnErrorData <- function(error, method, conf.level, init = NULL, data = NULL,
                            adjust = NULL, atom = NULL) {
  newTruncComp2(method = method,
               conf.level = conf.level,
               success = FALSE,
               error = error,
               init = init,
               data = data,
               adjust = adjust,
               atom = atom)
}
