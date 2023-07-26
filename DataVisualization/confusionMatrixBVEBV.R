library(caret)
library(writexl)

genList <- c("F2","F3","F4","F5","AYT","PYT")
datalist = list()
nreps = 15
x = 1
nCycle = 3

for (cycle in 1:nCycle) {
  for (gen in genList) {
    for (rep in 1:nreps){
    
  filename = paste("C", cycle, "_rrblup_random_F2_F5_F2_bvebv_snp_yield.rds", sep="") # read cycle file

  data <- readRDS(filename) 
  data <- data[[rep]] #take out one rep
  
  genData <- data[data$Gen==gen,] #take out one generation
  
  IDs <- as.data.frame(c(1:nrow(genData)))
  
  data1 <- cbind(genData,IDs) #order the data by bv
  colnames(data1) <- c("Gen","bv","ebv","IDs")
  data1 <- data1[order(data1$bv),]
  
  data2 <- cbind(genData,IDs) # order the data by ebv
  colnames(data2) <- c("Gen","bv","ebv","IDs")
  data2 <- data2[order(data2$ebv),]
  
  ## assign bv categories to each obs the data  
  cutoffBV =    data1[nrow(data1)*0.9,2]
  data1$bvCat <- as.factor(ifelse(data1$bv > cutoffBV, "1", "0"))
                
  ## put the data back in its original order
  data1Reordered <- data1[order(data1$IDs),]
  bvCat <- as.factor(data1Reordered$bvCat)

  # assign ebv catagories to each obs in the data
  cutoffEBV = data2[nrow(data2)*0.9,3]
  data2$ebvCat <- as.factor(ifelse(data2$ebv > cutoffEBV, "1", "0"))
  
  
  # put the data back in its original order
  data2Reordered <- data2[order(data2$IDs),]
  ebvCat <- as.factor(data2Reordered$ebvCat)

  # calculate CM
  CM <- confusionMatrix(data=bvCat, reference = ebvCat) # the CM is for ONE GENERATION OF ONE REP
  
  Acc <- as.data.frame(CM$overall) # pull accuracy
  Accuracy <- Acc[1,1]
  SensSpec <- as.data.frame(CM$byClass)
  Sensitivity <- SensSpec[1,1] #pull sensitivity
  Specificity <- SensSpec[2,1] #pull specificity
  
  df <-cbind(Accuracy, Sensitivity, Specificity) #bind values 
  rownames(df) = gen 
  datalist[[x]] = df #accumulate values 
  x=x+1 #move on to the next set of values 

  }
}  
}


cmValues <- do.call(rbind, datalist)

reps = 15
r = 1
cmMeans = list()
while (reps <= nrow(cmValues)) {
  set = as.matrix(cmValues[(reps-(reps-1)):reps,]) # pull out the first 15 obs (each rep for a given generation)
  repMean = colMeans(set) # find the means 
  cmMeans[[r]] = repMean # collect the means 
  r = r +1
  reps = 15 * r #change 15 to number of reps if number of reps is not 15. CMda
}

CMdata = do.call(rbind, cmMeans) # create DF
rownames(CMdata) = rep(genList,3) # name the generations 


#The results of the CM will tell us how well we were able to correctly identify the top 10% of individuals as top
#and the bottom 90% of individuals as bottom





