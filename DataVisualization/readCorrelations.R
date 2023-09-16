nCycle = 3
datalist = list()

#gather data from all consecutive cycles
for (cycle in 1:nCycle) {
  filename = paste("C", cycle, "_rrblup_random_F2_F5_F2_cors_snp_yield.csv", sep="")
  data = read.csv(filename)
  datalist[[cycle]] = data
}
correlationValues <- do.call(rbind, datalist)
correlationValues = correlationValues[ , colSums(is.na(correlationValues))==0] #remove columns(reps) containing an NA

# take the mean across reps 
values = as.data.frame(correlationValues[,-1])
means = rowMeans(values)
gens = as.data.frame(correlationValues[,1])
resultsCORs = cbind(gens,means) # this DF has each correlation value in consecutive order from the beginning of C1 to the end of C3

#find standard deviation of cumulative mean across reps (columns)
stdMat = matrix(nrow=nrow(cumulativeCor), ncol=1)
for (x in 1:nrow(cumulativeCor)) {
  row = cumulativeCor[x,]
  sd = sd(row)
  stdMat[x,1] = sd
}
std = as.data.frame(stdMat)

FINAL = cbind(resultsCors, std)
write_xlsx(FINAL,"resultsCors.xlsx")