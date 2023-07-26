# depending on ncycles and nreps, this can be a very slow script
# consider running in the cluster if your machine doesn't have a ton of memory

ncycles = 3
nreps=15
nSNP = 3547
genList=c("NP","F1","F2","F3","F4","F5","PYT","AYT","Variety")

freqList <- list()
datalist <- list()

for (cycle in 1:ncycles){
  for (i in 1:nreps) {
    for(gen in genList){
      
      filename = paste("C", cycle, "_rrblup_random_F2_F5_F2_alleles_snp_yield.rds", sep="")
      genotypes <- readRDS(filename) #read in each cycle's SNP data
      
      repMat <- genotypes[[i]] #pull out one rep
      genMat1 <- repMat[repMat$Gen==gen,] # pull out one generation
      genMat2 <- genMat1[,-1] # remove label column
      alleleSum <- as.matrix(colSums(genMat2)) # sum copy of alleles coded 1
      x <- as.matrix(alleleSum/nrow(genMat2)) # divide by number of alleles total
      
      Major <- as.matrix(ifelse(x<=0.5, 1-x, x*1)) #if the frequency is less than 0.5, then the other allele is major so we will calculate 1-x to make all Freqs the Major Freq
      
      freqList[[gen]] <- Major #add generation frequency to the freqList
      #cycle through all gens for a given rep
    }
    freqDF <- do.call(cbind, freqList) #turn list into DF
    datalist[[i]] <- freqList #add the freqList for 1 rep to the dataList
    cat("finished rep", i,"of", nreps, "for cycle", cycle,'\n')
  }
  datalistDF = as.data.frame(do.call("rbind",datalist)) # turn to DF
  saveRDS(datalist,paste("datalistC",cycle,".rds", sep="")) #we will have one dataList for each cycle
  cat("finished cycle", cycle, "of",ncycles,'\n')
}


# datalist is a list of data frames. there is one datalist for each cycle. 
# each data frame represents one rep
# each column in a data frame represents a generation
# each row represents the frequency of each SNP

#next we need the mean frequency across reps 

MeanFreq = matrix(nrow=nSNP, ncol=length(genList)) 
AllData = list()

for (cycle in 1:ncycles) {
    for (gen in (1:length(genList))){
    datafile = paste("datalistC",cycle,".rds", sep="") #read in the datalist file for each cycle
    DF = readRDS(datafile) 
    DF = t(do.call(cbind,DF)) # turn to DF with each generation (all reps) representing 1 column
    colnames(DF) = c(1:length(genList)) # number the columns
    genDF <- DF[,gen] #pull out one generation from the cycle
    genDF1 <- do.call("cbind",genDF) # turn to DF 
    genMeans <- as.matrix(rowMeans(genDF1)) #find the mean frequency across all reps 
    MeanFreq[,gen] = genMeans # add to the MeanFreq matrix
    }
  assign(paste0("meanFreqsC",cycle),MeanFreq) 
  saveRDS(MeanFreq,paste("MeanFreqC",cycle,".rds",sep=""))
  } 


FinalAlleleFreq =list()

for (x in 1:ncycles){
  data = paste0("meanFreqC",cycle, ".rds")
  data = as.data.frame(readRDS(data))# read in the mean freqs for each cycle 
  FinalAlleleFreq[[x]] = data # accumulate them in the FinalAlleleFreq list
}

Alldata = do.call("cbind",FinalAlleleFreq) # now we have a data frame with the mean freq for reps for each consecutive gen
colnames(Alldata) = c(rep(genList, times=ncycles))
                      
###
BOC <- as.data.frame(Alldata[,1]) # take the fist gen of the cycle
BOC$base <- rep(c("C1"), times = nrow(BOC))
colnames(BOC) <- c("freq","stage")
                      
MOC <- as.data.frame(Alldata[,ncol(Alldata*0.5)]) #take the middle of the 3 cycles 
MOC$base <- rep(c("C2"), times = nrow(MOC))
colnames(MOC) <- c("freq","stage")
                      
EOC <- as.data.frame(Alldata[,ncol(Alldata)]) # take the last gen in the cycle 
EOC$base <- rep(c("C3"), times = nrow(EOC))
colnames(EOC) <- c("freq","stage")
                      
alldata <- rbind(BOC,MOC,EOC) # we will look at just these 3 gens for simplicity
                      
    
library(ggplot2)
library(lattice)

# Density plot with semi-transparent fill
p<-ggplot(alldata, aes(x=freq, fill=stage)) + #plot the frequency and color by stage in the cycle 
          geom_density(alpha=0.5)

p
                      
# calculate overlap proportion
x <- list(C1=MOC$freq, C3=EOC$freq) # calculate the overlap between any 2stages 
                    

overlap(x, plot=FALSE, type="2")                      
                      
