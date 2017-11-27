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

d <- simTruncData(n = 100, mu0 = 3, mu1 = 3.5, pi0 = 0.6, pi1 = 0.5)

truncComp(d$Y, d$A, d$Z, method="LRT")

truncComp(d$Y, d$A, d$Z, method="SPLRT")
```
