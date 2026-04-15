library(TruncComp)
library(parallel)

simulatePowerData <- \(scenario, n) {
  R <- c(rep(0, n), rep(1, n))
  
  if(scenario == 1) {
    A <- sample(0:1, n*2, replace = TRUE, prob = c(0.65, 0.35))
    d <- data.frame(R = R, A = A, Y = NA)
    d$Y[d$A == 0] <- 0
    d$Y[d$A == 1 & d$R == 0] <- rnorm(sum(d$A == 1 & d$R == 0), 3, 1)
    d$Y[d$A == 1 & d$R == 1] <- rnorm(sum(d$A == 1 & d$R == 1), 3 + 0.5, 1)
  } else if(scenario == 2) {
    A <- c(sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5)),
           sample(0:1, n, replace = TRUE, prob = c(0.35, 0.65)))
    d <- data.frame(R = R, A = A, Y = NA)
    d$Y[d$A == 0] <- 0
    d$Y[d$A == 1] <- rnorm(sum(d$A == 1), 3.5, 1)
  } else if(scenario == 3) {
    A <- c(sample(0:1, n, replace = TRUE, prob = c(0.6, 0.4)),
           sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3)))
    d <- data.frame(R = R, A = A, Y = NA)
    d$Y[d$A == 0] <- 0 
    d$Y[d$A == 1 & d$R == 0] <- rnorm(sum(d$A == 1 & d$R == 0), 3, 1)
    d$Y[d$A == 1 & d$R == 1] <- rnorm(sum(d$A == 1 & d$R == 1), 3 + 0.5, 1)
  } else if(scenario == 4) {
    A <- sample(0:1, n*2, replace = TRUE, prob = c(0.65, 0.35))
    d <- data.frame(R = R, A = A, Y = NA)
    d$Y[d$A == 0] <- 0
    d$Y[d$A == 1 & d$R == 0] <- rbeta(sum(d$A == 1 & d$R == 0), 1, 1)
    d$Y[d$A == 1 & d$R == 1] <- rbeta(sum(d$A == 1 & d$R == 1), 1, 1 - 0.7)
  }else {
    stop("Invalid scenario")
  }
  
  list(wilcox = wilcox.test(d$Y ~ d$R)$p.value,
       ttest = t.test(d$Y ~ d$R)$p.value,
       LRT = truncComp(Y ~ R, atom = 0, data = d, method="LRT")$p,
       SPLRT = truncComp(Y ~ R, atom = 0, data = d, method="SPLRT")$p)
}

simulatePower <- \(scenario, n, R, mc.cores) {
  out <- mclapply(1:R, \(r) {
    if(r %% 1000 == 0) {
      cat(r, "\n")
    }
    simulatePowerData(scenario, n)
  }, mc.cores = mc.cores)
  c("Wilcoxon" = mean(sapply(out, \(q) q$wilcox < 0.05)),
    "T-test" = mean(sapply(out, \(q) q$ttest < 0.05)),
    "LRT" = mean(sapply(out, \(q) q$LRT < 0.05)),
    "SPLRT" = mean(sapply(out, \(q) q$SPLRT < 0.05)))
}

simulateNullData <- \(scenario, n, p) {
  R <- c(rep(0, n), rep(1, n))
  A <- c(sample(0:1, 2*n, replace = TRUE, prob = c(1 - p, p)))
  d <- data.frame(R = R, A = A, Y = NA)
  
  if(scenario == 1) {
    d$Y[d$A == 0] <- 0
    d$Y[d$A == 1] <- rnorm(sum(d$A == 1), 3, 1)
  } else if(scenario == 2) {
    d$Y[d$A == 0] <- 0
    d$Y[d$A == 1] <- rbeta(sum(d$A == 1), 1, 1)
  }

  list(wilcox = wilcox.test(d$Y ~ d$R)$p.value,
       ttest = t.test(d$Y ~ d$R)$p.value,
       LRT = truncComp(Y ~ R, atom = 0, data = d, method="LRT")$p,
       SPLRT = truncComp(Y ~ R, atom = 0, data = d, method="SPLRT")$p)
}

simulateTypeI <- \(scenario, n, p, R, mc.cores) {
  out <- mclapply(1:R, \(r) {
    if(r %% 1000 == 0) {
      cat(r, "\n")
    }
    simulateNullData(scenario, n, p)
  }, mc.cores = mc.cores)
  c("Wilcoxon" = mean(sapply(out, \(q) q$wilcox < 0.05)),
    "T-test" = mean(sapply(out, \(q) q$ttest < 0.05)),
    "LRT" = mean(sapply(out, \(q) q$LRT < 0.05)),
    "SPLRT" = mean(sapply(out, \(q) q$SPLRT < 0.05)))
}

##################################################################
# Power simulation
##################################################################
nSeq <- seq(50, 350, 25)

set.seed(12345)
power1 <- t(sapply(nSeq, \(n) {
  cat("****", n, "\n")
  simulatePower(1, n, 2*10^4, 128)
}))
rownames(power1) <- nSeq

set.seed(12345)
power2 <- t(sapply(nSeq, \(n) {
  cat("****", n, "\n")
  simulatePower(2, n, 2*10^4, 128)
}))
rownames(power2) <- nSeq

set.seed(12345)
power3 <- t(sapply(nSeq, \(n) {
  cat("****", n, "\n")
  simulatePower(3, n, 2*10^4, 128)
}))
rownames(power3) <- nSeq

set.seed(12345)
power4 <- t(sapply(nSeq, \(n) {
  cat("****", n, "\n")
  simulatePower(4, n, 2*10^4, 128)
}))
rownames(power4) <- nSeq


##################################################################
# Simulation of type I error
##################################################################
nSeqNull <- seq(100, 200, 400, 800, 1600)
pSeq <- c(0.2, 0.4, 0.6, 0.8)

set.seed(12345)
null1 <- lapply(pSeq, \(p) {
  cat("!!!!!", p, "\n")
  out <- t(sapply(nSeqNull, \(n) {
    cat("****", n, "\n")
    simulateTypeI(1, n, p, 10*10^4, 128)
  }))
  rownames(out) <- nSeqNull
  out
})

set.seed(12345)
null2 <- lapply(pSeq, \(p) {
  cat("!!!!!", p, "\n")
  out <- t(sapply(nSeqNull, \(n) {
    cat("****", n, "\n")
    simulateTypeI(2, n, p, 10*10^4, 128)
  }))
  rownames(out) <- nSeqNull
  out
})

##################################################################
# Save results
##################################################################
save(nSeq, power1, power2, power3, power4,
     nSeqNull, pSeq, null1, null2,
     file="simulationResults.RData")
