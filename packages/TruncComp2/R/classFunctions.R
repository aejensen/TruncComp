#' Summarize a truncated-comparison fit
#'
#' Prints the estimated treatment contrasts, their confidence intervals, and the
#' joint likelihood-ratio result for a fitted `"trunc_comp_fit"` object.
#'
#' @param object A `"trunc_comp_fit"` object returned by [truncComp()].
#' @param ... Unused additional arguments.
#' @return The input object, returned invisibly.
#' @details
#' For successful unadjusted fits the summary also reports the fitted `delta`
#' point estimate. Any confidence interval for `delta` is computed on demand
#' through [confint.trunc_comp_fit()].
#' @examples
#' library(TruncComp2)
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' d <- simulate_truncated_data(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "SPLRT")
#' summary(fit)
#' @rdname summary
#' @export
summary.trunc_comp_fit <- function(object, ...) {
  cat("Method:", object$method, "\n")
  cat("Confidence level = ", object$conf_level * 100, "%\n\n", sep = "")
  if(!is.null(object$adjust)) {
    cat("Adjusted for:", object$adjust, "\n\n")
  }

  if(isTRUE(object$success)) {
    estimates <- matrix(NA_real_, 2, 3)
    estimates[1, 1] <- object$mu_delta
    estimates[1, 2:3] <- object$mu_delta_ci
    estimates[2, 1] <- object$alpha_delta
    estimates[2, 2:3] <- object$alpha_delta_ci

    rownames(estimates) <- c(
      "Difference in means among the observed",
      "Odds ratio of being observed"
    )
    colnames(estimates) <- c("Estimate", "CI Lower", "CI Upper")

    if(is.finite(object$delta)) {
      estimates <- rbind(
        estimates,
        "Combined-outcome mean difference (delta)" = c(object$delta, NA_real_, NA_real_)
      )
    }

    cat("Treatment contrasts\n")
    print.default(estimates)
    cat("\nJoint test statistic:", object$statistic)
    cat("\np-value:", object$p, "\n")
  } else {
    cat("The estimation procedure failed.\n")
    cat("The error message was:", object$error, "\n")
  }
  cat("\n")

  invisible(object)
}

#' Print a truncated-comparison fit
#'
#' Prints a concise overview of a fitted `"trunc_comp_fit"` object.
#'
#' @param x A `"trunc_comp_fit"` object returned by [truncComp()].
#' @param ... Unused additional arguments.
#' @return The input object, returned invisibly.
#' @examples
#' library(TruncComp2)
#' f0 <- function(n) stats::rnorm(n, 3, 1)
#' f1 <- function(n) stats::rnorm(n, 3.5, 1)
#' d <- simulate_truncated_data(n = 40, f0 = f0, f1 = f1, pi0 = 0.6, pi1 = 0.5)
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
#' print(fit)
#' @rdname print
#' @export
print.trunc_comp_fit <- function(x, ...) {
  cat("<trunc_comp_fit>\n")
  cat("Method:", x$method, "\n")
  cat("Success:", isTRUE(x$success), "\n")

  if(isTRUE(x$success)) {
    estimates <- stats::setNames(
      c(x$mu_delta, x$alpha_delta, x$delta),
      c("mu_delta", "alpha_delta", "delta")
    )
    print.default(estimates)
    cat("\nUse summary() for confidence intervals and the joint test.\n")
  } else {
    cat("Error:", x$error, "\n")
  }

  invisible(x)
}

#' Extract fitted treatment-effect estimates
#'
#' Returns the primary fitted treatment-effect estimates from a
#' `"trunc_comp_fit"` object.
#'
#' @param object A `"trunc_comp_fit"` object returned by [truncComp()].
#' @param ... Unused additional arguments.
#' @return A named numeric vector with entries `mu_delta`, `alpha_delta`, and
#'   `delta`.
#' @examples
#' library(TruncComp2)
#' d <- load_trunc_comp2_example()
#' fit <- truncComp(Y ~ R, atom = 0, data = d, method = "LRT")
#' coef(fit)
#' @export
coef.trunc_comp_fit <- function(object, ...) {
  stats::setNames(
    c(object$mu_delta, object$alpha_delta, object$delta),
    c("mu_delta", "alpha_delta", "delta")
  )
}
