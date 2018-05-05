# TruncComp
Development version of the R package TruncComp for two-sample comparison of truncated continuous outcomes.

To install the development version of TruncComp run the following commands from within R

```{r}
library(devtools)
install_github('aejensen/TruncComp')
```

# Example
```{r}
library(TruncComp)

#Define the two distributions for the observed data
f0 <- function(n) stats::rnorm(n, 3, 1)
f1 <- function(n) stats::rnorm(n, 3.5, 1)

#Define probabilities of being observed
pi0 <- 0.35
pi1 <- 0.6

#Simulate data
d <- TruncComp::simulateTruncatedData(25, f0, f1, pi0, pi1)

#Estimate parameters using the semi-parametric method
model <- truncComp(Y ~ R, atom = 0, data = d, method="SPLRT")
summary(model)

#Get marginal confidence intervals
confint(model, type="marginal")

#Get simultaneous confidence region
confint(model, type="simultaneous", plot=TRUE, resolution = 10)
```
