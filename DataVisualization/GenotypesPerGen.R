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

  for (x in 1:nrow){ # for every inividual
    for (y in 1:ncol){ # for every locus
      for (z in 1:nreps){ #for every rep
        rep = allelelist[[z]] #pull one rep
        cellVals = matrix(nrow=1, ncol=nreps) # create empty matrix
        cellVals[[x,y]] = rep[[x,y]] #fill matrix with values for an individual at each locus
        mode = which.max(cellVals) #pull out the mode for that locus across reps 
        repMode[x,y] = mode } #assign genotype to final matrix
    } # finish one row, every locus and then move on to the next row 
  }

        
        
        
