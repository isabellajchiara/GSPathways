#collect results from all cycles into a DF

nCycle = 3
datalist = list()

for (cycle in 1:nCycle) {
  filename = paste("C", cycle, "_rrblup_random_trainAtF5_trainWithF2_F2Parents_gvs_snp_yield.csv", sep="")
  data = read.csv(filename)
  datalist[[cycle]] = data
}
geneticValues <- do.call(rbind, datalist)
geneticValues <- geneticValues[-c(11,12,22,23,33),]

values = as.matrix(geneticValues[,-1])
cumulativeGain = matrix(nrow = (nrow(values)-1), ncol = (ncol(values)))

#find the cumulative gain across each generation from the start of the sim
for (x in 0:(nrow(values)-1)) {
  gen1 = as.numeric((values[x+1,]))
  gen2 = as.numeric((values[1,]))
  change = gen1 - gen2
  cumulativeGain[x,] <- change
}

#find the mean cumulative gain across reps 
gens = as.data.frame(geneticValues[-1,1])
cumulativeGain = as.data.frame(cumulativeGain)
meanGain = as.data.frame(rowMeans(cumulativeGain))
resultsGVs = cbind(gens, meanGain) # this DF has each genetic value in consecutive order from the beginning of C1 to the end of C3

#find standard deviation of cumulative mean across reps (columns)
stdMat = matrix(nrow=nrow(cumulativeGain), ncol=1)
for (x in 1:nrow(cumulativeGain)) {
  row = cumulativeGain[x,]
  sd = sd(row)
  stdMat[x,1] = sd
}
std = as.data.frame(stdMat)

FINAL = cbind(resultsGVs, std)
write_xlsx(FINAL,"resultsGVs.xlsx")
