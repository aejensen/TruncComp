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

print.TruncComp2 <- function(x, ...) {
  summary(x, ...)
}
