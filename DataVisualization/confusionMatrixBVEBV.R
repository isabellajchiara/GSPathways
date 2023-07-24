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
  cutoffBV =    data1[nrow(data1)*0.8,2]
  data1$bvCat <- as.factor(ifelse(data1$bv > cutoffBV, "1", "0"))
                
  ## put the data back in its original order
  data1Reordered <- data1[order(data1$IDs),]
  bvCat <- as.factor(data1Reordered$bvCat)

  # assign ebv catagories to each obs in the data
  cutoffEBV = data2[nrow(data2)*0.8,3]
  data2$ebvCat <- as.factor(ifelse(data2$ebv > cutoffEBV, "1", "0"))
  
  
  # put the data back in its original order
  data2Reordered <- data2[order(data2$IDs),]
  ebvCat <- as.factor(data2Reordered$ebvCat)

  # calculate CM
  CM <- confusionMatrix(data=bvCat, reference = ebvCat) # the CM is for ONE GENERATION OF ONE REP
  
  Acc <- as.data.frame(CM$overall)
  Accuracy <- Acc[1,1]
  SensSpec <- as.data.frame(CM$byClass)
  Sensitivity <- SensSpec[1,1]
  Specificity <- SensSpec[2,1]
  
  df <-cbind(Accuracy, Sensitivity, Specificity)
  rownames(df) = gen
  datalist[[x]] = df
  x=x+1

  }
}  
}


cmValues <- do.call(rbind, datalist)

C1 = cmValues[1:90,]
C2 = cmValues[91:180,]
C3 = cmValues[181:270,]

confusionMatrix = cbind(C1,C2,C3)





