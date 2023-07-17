library(caret)
library(writexl)

genList <- c("F2","F3","F4","F5")

for (gen in genList) {
  
  data <- readRDS("8C1rrblup_rd_bvebv_snp_yield.rds")
  data <- data[,1:3]
  
  data <- data[data$Gen==gen,]
  
  bv <- as.vector(range(data$bv))
  ebv <- as.vector(range(data$ebv))
  
  IDs <- as.data.frame(c(1:nrow(data)))
  
  data1 <- cbind(data,IDs)
  colnames(data1) <- c("Gen","bv","ebv","IDs")
  data1 <- data1[order(data1$bv),]
  
  data2 <- cbind(data,IDs)
  colnames(data2) <- c("Gen","bv","ebv","IDs")
  data2 <- data2[order(data2$ebv),]
  
  
  data1$bvCat <- as.factor(ifelse(data1$bv < data1[nrow(data)*0.1,2], "0",
                                  ifelse(data1$bv < data1[nrow(data)*0.2,2],"0",
                                         ifelse(data1$bv < data1[nrow(data)*0.3,2],"0",
                                                ifelse(data1$bv < data1[nrow(data)*0.4,2],"0",
                                                       ifelse(data1$bv < data1[nrow(data)*0.5,2],"0",
                                                              ifelse(data1$bv < data1[nrow(data)*0.6,2],"0",
                                                                     ifelse(data1$bv < data1[nrow(data)*0.7,2],"1","1"))))))))
  
  data1Reordered <- data1[order(data1$IDs),]
  bvCat <- as.factor(data1Reordered$bvCat)
  
  
  data2$ebvCat <- as.factor(ifelse(data2$ebv < data2[nrow(data)*0.1,3], "0",
                                   ifelse(data2$ebv < data2[nrow(data)*0.2,3],"0",
                                          ifelse(data2$ebv < data2[nrow(data)*0.3,3],"0",
                                                 ifelse(data2$ebv < data2[nrow(data)*0.4,3],"0",
                                                        ifelse(data2$ebv < data2[nrow(data)*0.5,3],"0",
                                                               ifelse(data2$ebv < data2[nrow(data)*0.6,3],"0",
                                                                      ifelse(data2$ebv < data1[nrow(data)*0.7,3],"1","1"))))))))
  
  
  
  data2Reordered <- data2[order(data2$IDs),]
  ebvCat <- as.factor(data2Reordered$ebvCat)
  
  CM <- confusionMatrix(data=bvCat, reference = ebvCat)
  
  Acc<- as.data.frame(CM$overall)
  Accuracy <- Acc[1,1]
  SensSpec <- as.data.frame(CM$byClass)
  Sensitivity <- SensSpec[1,1]
  Specificity <- SensSpec[2,1]
  
  df <- paste0("C1",gen)
  df <-cbind(Accuracy, Sensitivity, Specificity)
  assign(paste0("C1",gen), df)
  

  ###
  data <- readRDS("8C2rrblup_rd_bvebv_snp_yield.rds")
  data <- data[,1:3]
  data <- data[data$Gen==gen,]
  
  bv <- as.vector(range(data$bv))
  ebv <- as.vector(range(data$ebv))
  
  IDs <- as.data.frame(c(1:nrow(data)))
  
  data1 <- cbind(data,IDs)
  colnames(data1) <- c("Gen","bv","ebv","IDs")
  data1 <- data1[order(data1$bv),]
  
  data2 <- cbind(data,IDs)
  colnames(data2) <- c("Gen","bv","ebv","IDs")
  data2 <- data2[order(data2$ebv),]
  
  
  data1$bvCat <- as.factor(ifelse(data1$bv < data1[nrow(data)*0.1,2], "0",
                                  ifelse(data1$bv < data1[nrow(data)*0.2,2],"0",
                                         ifelse(data1$bv < data1[nrow(data)*0.3,2],"0",
                                                ifelse(data1$bv < data1[nrow(data)*0.4,2],"0",
                                                       ifelse(data1$bv < data1[nrow(data)*0.5,2],"0",
                                                              ifelse(data1$bv < data1[nrow(data)*0.6,2],"0",
                                                                     ifelse(data1$bv < data1[nrow(data)*0.7,2],"1","1"))))))))
  
  data1Reordered <- data1[order(data1$IDs),]
  bvCat <- as.factor(data1Reordered$bvCat)
  
  
  data2$ebvCat <- as.factor(ifelse(data2$ebv < data2[nrow(data)*0.1,3], "0",
                                   ifelse(data2$ebv < data2[nrow(data)*0.2,3],"0",
                                          ifelse(data2$ebv < data2[nrow(data)*0.3,3],"0",
                                                 ifelse(data2$ebv < data2[nrow(data)*0.4,3],"0",
                                                        ifelse(data2$ebv < data2[nrow(data)*0.5,3],"0",
                                                               ifelse(data2$ebv < data2[nrow(data)*0.6,3],"0",
                                                                      ifelse(data2$ebv < data1[nrow(data)*0.7,3],"1","1"))))))))
  
  
  
  data2Reordered <- data2[order(data2$IDs),]
  ebvCat <- as.factor(data2Reordered$ebvCat)
  
  CM <- confusionMatrix(data=bvCat, reference = ebvCat)
  
  Acc<- as.data.frame(CM$overall)
  Accuracy <- Acc[1,1]
  SensSpec <- as.data.frame(CM$byClass)
  Sensitivity <- SensSpec[1,1]
  Specificity <- SensSpec[2,1]
  
  df <- paste0("C2",gen)
  df <-cbind(Accuracy, Sensitivity, Specificity)
  assign(paste0("C2",gen), df)
  
  ###################
  
  data <- readRDS("8C3rrblup_rd_bvebv_snp_yield.rds")
  data <- data[,1:3]
  data <- data[data$Gen==gen,]
  
  
  bv <- as.vector(range(data$bv))
  ebv <- as.vector(range(data$ebv))
  
  IDs <- as.data.frame(c(1:nrow(data)))
  
  data1 <- cbind(data,IDs)
  colnames(data1) <- c("Gen","bv","ebv","IDs")
  data1 <- data1[order(data1$bv),]
  
  data2 <- cbind(data,IDs)
  colnames(data2) <- c("Gen","bv","ebv","IDs")
  data2 <- data2[order(data2$ebv),]
  
  
  data1$bvCat <- as.factor(ifelse(data1$bv < data1[nrow(data)*0.1,2], "0",
                                  ifelse(data1$bv < data1[nrow(data)*0.2,2],"0",
                                         ifelse(data1$bv < data1[nrow(data)*0.3,2],"0",
                                                ifelse(data1$bv < data1[nrow(data)*0.4,2],"0",
                                                       ifelse(data1$bv < data1[nrow(data)*0.5,2],"0",
                                                              ifelse(data1$bv < data1[nrow(data)*0.6,2],"0",
                                                                     ifelse(data1$bv < data1[nrow(data)*0.7,2],"1","1"))))))))
  
  data1Reordered <- data1[order(data1$IDs),]
  bvCat <- as.factor(data1Reordered$bvCat)
  
  
  data2$ebvCat <- as.factor(ifelse(data2$ebv < data2[nrow(data)*0.1,3], "0",
                                   ifelse(data2$ebv < data2[nrow(data)*0.2,3],"0",
                                          ifelse(data2$ebv < data2[nrow(data)*0.3,3],"0",
                                                 ifelse(data2$ebv < data2[nrow(data)*0.4,3],"0",
                                                        ifelse(data2$ebv < data2[nrow(data)*0.5,3],"0",
                                                               ifelse(data2$ebv < data2[nrow(data)*0.6,3],"0",
                                                                      ifelse(data2$ebv < data1[nrow(data)*0.7,3],"1","1"))))))))
  
  
  data2Reordered <- data2[order(data2$IDs),]
  ebvCat <- as.factor(data2Reordered$ebvCat)
  
  CM <- confusionMatrix(data=bvCat, reference = ebvCat)
  
  Acc<- as.data.frame(CM$overall)
  Accuracy <- Acc[1,1]
  SensSpec <- as.data.frame(CM$byClass)
  Sensitivity <- SensSpec[1,1]
  Specificity <- SensSpec[2,1]
  
  df <- paste0("C3",gen)
  df <-cbind(Accuracy, Sensitivity, Specificity)
  assign(paste0("C3",gen), df)
  
}
  


####
  
  C1 <- rbind(C1F2,C1F3,C1F4,C1F5)
  C2 <- rbind(C2F2,C2F3,C2F4,C2F5)
  C3 <- rbind(C3F2,C3F3,C3F4,C3F5)

results <- as.data.frame(cbind(C1,C2,C3))
write_xlsx(results,"confusionMatrix.xlsx")



