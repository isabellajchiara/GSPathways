ncycles = 3
nreps=3
nSNP = 3547
genList=c("NP","F1","F2","F3","F4","F5","PYT","AYT","Variety")

freqList <- list()
datalist <- list()

for (cycle in 1:ncycles){
  for (i in 1:nreps) {
    for(gen in genList){
      
      filename = paste("1C", cycle, "_svm_random_alleles_snp_yield.rds", sep="")
      genotypes <- readRDS(filename) #read in each cycle's SNP data
      
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
  saveRDS(datalist,paste("datalistC",cycle,".rds", sep=""))
}

# datalist is a list of data frames
# each data frame represents one rep
# each column represents a generation
# each row represents the frequency of each SNP

#next we need the mean frequency across reps 

MeanFreq = matrix(nrow=nSNP, ncol=length(genList)
AllData = list()

for (cycle in 1:ncycles) {
  for (gen in (1:length(genList))){
    datafile = paste("datalistC",cycle,".rds", sep="")
    DF = readRDS(datafile)
    DF = t(do.call(cbind,DF))
    colnames(DF) = c(1:length(genList))
    genDF <- DF[,gen]
    genDF1 <- do.call("cbind",genDF)
    genMeans <- as.matrix(rowMeans(genDF1))
    MeanFreq[,gen] = genMeans
  }
  assign(paste0("meanFreqsC",cycle),MeanFreq) 
  saveRDS(MeanFreq,paste("MeanFreqC",cycle,".rds",sep=""))
} 

FinalAlleleFreq =list()

for (x in 1:ncycles){
  data = paste0("meanFreqC",cycle, ".rds")
  data = as.data.frame(readRDS(data))
  FinalAlleleFreq[[x]] = data
}

Alldata = do.call("cbind",FinalAlleleFreq)
colnames(Alldata) = c(rep(genList, times=ncycles))
                      
###
BOC <- as.data.frame(Alldata[,1])
BOC$base <- rep(c("C1"), times = nrow(BOC))
colnames(BOC) <- c("freq","stage")
                      
MOC <- as.data.frame(Alldata[,ncol(Alldata*0.5)])
MOC$base <- rep(c("C2"), times = nrow(MOC))
colnames(MOC) <- c("freq","stage")
                      
EOC <- as.data.frame(Alldata[,ncol(Alldata)])
EOC$base <- rep(c("C3"), times = nrow(EOC))
colnames(EOC) <- c("freq","stage")
                      
alldata <- rbind(BOC,MOC,EOC)
                      
    
library(ggplot2)
library(lattice)

# Density plot with semi-transparent fill
p<-ggplot(alldata, aes(x=freq, fill=stage)) +
          geom_density(alpha=0.5)

p
                      
# calcjulate overlap proportion
x <- list(C1=MOC$freq, C3=EOC$freq)
                    

overlap(x, plot=FALSE, type="2")                      
                      
