nCycle = 3
datalist = list()

#gather data from all consecutive cycles
for (cycle in 1:nCycle) {
  filename = paste("C", cycle, "_rrblup_random_F2_F5_F2_cors_snp_yield.csv", sep="")
  data = read.csv(filename)
  datalist[[cycle]] = data
}
varianceValues <- do.call(rbind, datalist)
varianceValues = varianceValues[ , colSums(is.na(varianceValues))==0] #remove columns(reps) containing an NA

values = as.matrix(varianceValues[,-1])
cumulativeVar = matrix(nrow = (nrow(values)-1), ncol = (ncol(values)))

#find the cumulative delta variance across each generation from the start of the sim
for (x in 1:(nrow(values)-1)) {
  gen1 = as.numeric((values[x+1,]))
  gen2 = as.numeric((values[1,]))
  change = gen1 - gen2
  cumulativeVar[x,] <- change
}

#find the mean cumulative delta variance across reps 
gens = as.data.frame(varianceValuesValues[-c(1,10,20),1])
cumulativeVar = as.data.frame(cumulativeVar)
meanVar = as.data.frame(rowMeans(cumulativeVar))
meanVar = meanVar[-c(10,20),]
resultsVars = cbind(gens, meanVar) # this DF has each variance value in consecutive order from the beginning of C1 to the end of C3


