ncycles = 1
nreps=3
nSNP = 3547
genList=c("NP","F1","F2","F3","F4","F5","PYT","AYT","Variety")

freqList <- list()
datalist <- list()

for (cycle in 1:ncycles){
  for (i in 1:nreps) {
    for(gen in genList){
      
      genotypes <- readRDS((paste("1C",cycle,"_rrblup_rd_alleles_snp_yield.rds", sep=""))) #read in each cycle's SNP data
    
      repMat <- genotypes[[i]]
      genMat1 <- repMat[repMat$Gen==gen,]
      genMat2 <- genMat1[,-1]
      alleleSum <- as.matrix(colSums(genMat2))
      alleleFreq <- as.data.frame(alleleSum/nrow(genMat2))
    
      MAJMIN <- alleleSum
      MAJMIN[MAJMIN < (nrow(genMat2)/2)] <- 0
      MAJMIN[MAJMIN > (nrow(genMat2)/2)] <- 1
    
      majAllele <- cbind(MAJMIN, alleleFreq)
      x <- majAllele[,2]
    
      MajAlleleFreq <- as.matrix(ifelse(x<=0.5, x+1,x*(1)))
    
      freqList[[gen]] <- MajAlleleFreq
    }
    freqDF <- do.call(cbind, freqList)
    freqDF <- freqDF[,-1]
    datalist[[i]] <- freqList
  }
datalistDF = as.data.frame(do.call("rbind",datalist))
saveRDS(datalist,paste("datalistC",cycle,".rds"))
}

# datalist is a list of data frames
# each data frame represents one rep
# each column represents a generation
# each row represents the frequency of each SNP

#next we need the mean frequency across reps 

MeanFreq = matrix(nrow=nSNP, ncol=length(genList))
AllData = list()

for (cycle in 1:ncycles) {
  for (gen in (1:length(genList))){
    DF = readRDS(paste("datalistC",cycle,".rds"))
    DF = t(do.call(cbind,DF))
    colnames(DF) = c(1:length(genList))
    genDF <- DF[,gen]
    genDF1 <- do.call("cbind",genDF)
    genMeans <- as.matrix(rowMeans(genDF1))
    MeanFreq[,gen] = genMeans
  }
AllData[[gen]] <- MeanFreq
assign(paste0("meanFreqsC",cycle),MeanFreq) 

}  
      
    #bind 1-n column of each rep
    #calculate rowmeans
    #bind into one DF with means
    


###

AllFreq <- as.data.frame(cbind(C1MAF,C2MAF,C3MAF))

###
BOC <- as.data.frame(AllFreq[,1])
BOC$base <- rep(c("C1"), times = nrow(BOC))
colnames(BOC) <- c("freq","stage")

MOC <- as.data.frame(AllFreq[,14])
MOC$base <- rep(c("C2"), times = nrow(MOC))
colnames(MOC) <- c("freq","stage")

EOC <- as.data.frame(AllFreq[,27])
EOC$base <- rep(c("C3"), times = nrow(EOC))
colnames(EOC) <- c("freq","stage")

alldata <- rbind(MOC,EOC)


# Density plot with semi-transparent fill
p<-ggplot(alldata, aes(x=freq, fill=stage)) +
  geom_density(alpha=0.5)+
  xlim(0,0.15)+
  ylim(0,100)
p

# calcjulate overlap proportion
x <- list(C1=MOC$freq, C3=EOC$freq)

overlap(x, plot=FALSE, type="2")
