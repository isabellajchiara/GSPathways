# to pull out mode genotypes across reps for visualiztion of genotypes 

cycle = 1 # one cycle at a time for computational purposes
nreps=15 # how manyr eps 
nSNP = 3547 #how many sites
gen = "PYT" #which gen are you viewing

allelelist <- list()
datalist <- list()

for (i in 1:nreps) {
      filename = paste("C", cycle, "_ann_random_trainAtF5_trainWithF2_F2Parents_alleles_snp_yield.rds", sep="")
      genotypes <- readRDS(filename) #read in each cycle's SNP data
      repMat <- genotypes[[i]] #pull out one rep
      genMat1 <- repMat[repMat$Gen==gen,] # pull out one generation
      genMat2 <- genMat1[,-1] # remove label column
      allelelist[[i]] <- genMat2 #add generation frequency to the allelelist
      cat("finished rep", i,"of", nreps, "for cycle", cycle,'\n')
}

repMode = matrix(nrow=120,ncol=nSNP)
nrow = 120
ncol = nSNP

 x = 1
  while (x < nrow){
    while (y < ncol){
      while (z < nreps){
        y = 1
        z = 1
        rep = allelelist[[z]] 
        cellVals = matrix(nrow=1, ncol=nreps)
        geno = as.numeric(rep[x,y])
        cellVals[,y] == geno
        mode = as.numeric(max(cellVals))
        repMode[x,y] = mode 
        y = y + 1
      }
      z= z+1
    }
    x = X+1
  }

        
        
        
        
