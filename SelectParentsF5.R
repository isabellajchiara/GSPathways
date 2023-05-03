
## select and define parents for next cycle

newCycleSelections = selectInd(F5, 5, use="ebv")

cross = matrix(nrow = 20, ncol=2)
cross[,1] <- c(1,1,1,1,2,2,2,3,3,4,2,3,4,5,3,4,5,4,5,5)
cross[,2] <- c(2,3,4,5,3,4,5,4,5,5,1,1,1,1,2,2,2,3,3,4)
