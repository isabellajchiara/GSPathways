#collect results from all cycles into a DF

nCycle = 3
datalist = list()

for (cycle in 1:nCycle) {
  filename = paste("C", cycle, "_rrblup_random_F2_F5_F2_gvs_snp_yield.csv", sep="")
  data = read.csv(filename)
  datalist[[cycle]] = data
}
geneticValues <- do.call(rbind, datalist)
geneticValues <- geneticValues[-c(11,22,33),]

values = as.matrix(geneticValues[,-1])
cumulativeGain = matrix(nrow = (nrow(values)-1), ncol = (ncol(values)))

#find the cumulative gain across each generation from the start of the sim
for (x in 1:(nrow(values)-1)) {
  gen1 = as.numeric((values[x+1,]))
  gen2 = as.numeric((values[1,]))
  change = gen1 - gen2
  cumulativeGain[x,] <- change
}

#find the mean cumulative gain across reps 
gens = as.data.frame(geneticValues[-c(1,10,20),1])
cumulativeGain = as.data.frame(cumulativeGain)
meanGain = as.data.frame(rowMeans(cumulativeGain))
meanGain = meanGain[-c(10,20),]
resultsGVs = cbind(gens, meanGain) # this DF has each genetic value in consecutive order from the beginning of C1 to the end of C3

