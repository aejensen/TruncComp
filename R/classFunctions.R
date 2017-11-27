summary.TruncComp <- function(object, ...) {
  cat("Composite Estimation Output\n\n")

  cat("Estimation method:", object$method, "\n")
  cat("Confidence level = ", object$conf.level * 100, "%\n\n", sep="")

  if (object$success == TRUE) {
    cMat <- matrix(NA, 2, 3)
    cMat[1,1] <- object$muDelta
    cMat[1, 2:3] <- object$muDeltaCI
    cMat[2, 1] <- object$alphaDelta
    cMat[2, 2:3] <- object$alphaDeltaCI
    colnames(cMat) <- c("Estimate", "CI Lower", "CI Upper")
    rownames(cMat) <- c("Mean difference among the observed:", "Odds ratio of being observed:")

    cat("Treatment contrasts:\n")
    print.default(cMat)

    cat("\nJoint LRT statistic: W =", object$W)
    cat("\np-value: p =", object$p, "\n")
  } else {
    cat("The estimation procedure failed.\n")
    cat("The error message was:", object$error, "\n")
  }
  cat("\n")

  invisible(object)
}

print.TruncComp <- function(x, ...) {
  summary(x)
}


