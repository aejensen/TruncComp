rm(list=ls())
data <- simTruncData(50, 1, 1, 0.2, 0.2)

#
truncComp.default(data$Y, data$A, data$R, method="SPLRT")
truncComp.formula(Y ~ R, atom = 0, data = data, method="SPLRT")

#
truncComp.default(data$Y, data$A, data$R, method="LRT")
truncComp.formula(Y ~ R, atom = 0, data = data, method="LRT")


