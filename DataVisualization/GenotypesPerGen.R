# depending on ncycles and nreps, this can be a very slow script
# consider running in the cluster if your machine doesn't have a ton of memory

cycle = 1
nreps=15
nSNP = 3547
gen = "PYT"

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
    x = x+1
  }

        
        
        
        
