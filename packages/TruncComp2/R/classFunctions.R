#' Summarize a TruncComp2 fit
#'
#' Prints the estimated treatment contrasts, their confidence intervals, and the
#' joint likelihood-ratio result for a fitted `"TruncComp2"` object.
#'
#' @param object A `"TruncComp2"` object returned by [truncComp()].
#' @param ... Unused additional arguments.
#' @return The input object, returned invisibly.
#' @details
#' For successful unadjusted fits the summary includes the stored
#' `DeltaMarginalCI` and `DeltaProfileCI`. The projected interval
#' `DeltaProjectedCI` remains available on demand through
#' [confint.TruncComp2()].
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
    cMat <- matrix(NA, 2, 3)
    cMat[1,1] <- object$muDelta
    cMat[1, 2:3] <- object$muDeltaCI

    cMat[2, 1] <- object$alphaDelta
    cMat[2, 2:3] <- object$alphaDeltaCI

    if(is.finite(object$Delta)) {
      delta_rows <- list()
      if(length(object$DeltaMarginalCI) >= 2 && all(is.finite(object$DeltaMarginalCI[1:2]))) {
        delta_rows[["Delta (marginal):"]] <- c(object$Delta, object$DeltaMarginalCI)
      }
      if(length(object$DeltaProjectedCI) >= 2 && all(is.finite(object$DeltaProjectedCI[1:2]))) {
        delta_rows[["Delta (projected):"]] <- c(object$Delta, object$DeltaProjectedCI)
      }
      if(length(object$DeltaProfileCI) >= 2 && all(is.finite(object$DeltaProfileCI[1:2]))) {
        delta_rows[["Delta (profile likelihood):"]] <- c(object$Delta, object$DeltaProfileCI)
      }

      if(length(delta_rows) > 0) {
        cMat <- rbind(cMat, do.call(rbind, delta_rows))
      }
    }

    colnames(cMat) <- c("Estimate", "CI Lower", "CI Upper")
    rownames(cMat) <- c("Difference in means among the observed:",
                        "Odds ratio of being observed:",
                        rownames(cMat)[-(1:2)])

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
