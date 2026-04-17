#' Summarize a TruncComp2 fit
#'
#' Prints the estimated treatment contrasts, their confidence intervals, and the
#' joint likelihood-ratio result for a fitted `"TruncComp2"` object.
#'
#' @param object A `"TruncComp2"` object returned by [truncComp()].
#' @param ... Unused additional arguments.
#' @return The input object, returned invisibly.
#' @details
#' For successful unadjusted fits the summary also reports the fitted `Delta`
#' point estimate. Any confidence interval for `Delta` is computed on demand
#' through [confint.TruncComp2()].
#' @examples
#' library(TruncComp2)
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
#' summary(fit)
#' @rdname summary
#' @export
summary.TruncComp2 <- function(object, ...) {
  cat("Estimation method:", object$method, "\n")
  cat("Confidence level = ", object$conf.level * 100, "%\n\n", sep="")
  if(!is.null(object$adjust)) {
    cat("Adjusted for:", object$adjust, "\n\n")
  }

  if (isTRUE(object$success)) {
    cMat <- matrix(NA_real_, 2, 3)
    cMat[1,1] <- object$muDelta
    cMat[1, 2:3] <- object$muDeltaCI

    cMat[2, 1] <- object$alphaDelta
    cMat[2, 2:3] <- object$alphaDeltaCI

    if(is.finite(object$Delta)) {
      cMat <- rbind(cMat, c(object$Delta, NA_real_, NA_real_))
    }

    colnames(cMat) <- c("Estimate", "CI Lower", "CI Upper")
    rownames(cMat) <- c("Difference in means among the observed:",
                        "Odds ratio of being observed:",
                        if(is.finite(object$Delta)) "Combined-outcome mean difference (Delta):")

    cat("Treatment contrasts\n")
    print.default(cMat)

    cat("\nJoint test statistic: W =", object$W)
    cat("\np-value: p =", object$p, "\n")
  } else {
    cat("The estimation procedure failed.\n")
    cat("The error message was:", object$error, "\n")
  }
  cat("\n")

  invisible(object)
}

#' Print a TruncComp2 fit
#'
#' Delegates to [summary.TruncComp2()] so printing and summarizing use the same
#' package-specific display format.
#'
#' @param x A `"TruncComp2"` object returned by [truncComp()].
#' @param ... Unused additional arguments.
#' @return The input object, returned invisibly.
#' @examples
#' library(TruncComp2)
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' d <- simulateTruncatedData(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
#' print(fit)
#' @rdname print
#' @export
print.TruncComp2 <- function(x, ...) {
  summary(x, ...)
}
