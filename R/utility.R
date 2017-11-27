isDataOkay <- function(d) {
  dAlive <- subset(d, d$A == 1)

  yAlive1 <- dAlive$Y[dAlive$Z == 0]
  yAlive2 <- dAlive$Y[dAlive$Z == 1]

  if(length(yAlive1) < 2 | length(yAlive2) < 2) {
    return(FALSE)
  }

  TRUE
}

isValid <- function(truncComObj) {
  truncCombObj$success
}
