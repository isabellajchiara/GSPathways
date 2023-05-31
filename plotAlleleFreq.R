genotypes <- readRDS("1C1rrblup_rd_alleles_snp_yield.rds")

C1MAF <- matrix(nrow=3547, ncol=9)

parentFrequencies <- matrix(nrow=3547,ncol=14)
F1Frequencies <- matrix(nrow=3547,ncol=14)
F2Frequencies <- matrix(nrow=3547,ncol=14)
F3Frequencies <- matrix(nrow=3547,ncol=14)
F4Frequencies <- matrix(nrow=3547,ncol=14)
F5Frequencies <- matrix(nrow=3547,ncol=14)
PYTFrequencies <- matrix(nrow=3547,ncol=14)
AYTFrequencies <- matrix(nrow=3547,ncol=14)
VarietyFrequencies <- matrix(nrow=3547,ncol=14)


nreps=14

for (i in 1:nreps) {
  parents <- genotypes[[i]]
  parents <- parents[parents$Gen=="NP",]
  parents <- parents[,-1]
  alleleSum <- as.matrix(colSums(parents))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))

  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  
  parentFrequencies[,i] <- MajAlleleFreq
}
parentFrequencies <-rowMeans(parentFrequencies)
C1MAF[,1] <-parentFrequencies



for (i in 1:nreps) {
  F1 <- genotypes[[i]]
  F1 <- F1[F1$Gen=="F1",]
  F1 <- F1[,-1]
  alleleSum <- as.matrix(colSums(F1))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F1Frequencies[,i] <- MajAlleleFreq
}

F1Frequencies <-rowMeans(F1Frequencies)
C1MAF[,2] <-F1Frequencies

for (i in 1:nreps) {
  F2 <- genotypes[[i]]
  F2 <- F2[F2$Gen=="F2",]
  F2 <- F2[,-1]
  alleleSum <- as.matrix(colSums(F2))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F2Frequencies[,i] <- MajAlleleFreq
}
F2Frequencies <-rowMeans(F2Frequencies)
C1MAF[,3] <-F2Frequencies

for (i in 1:nreps) {
  F3 <- genotypes[[i]]
  F3 <- F3[F3$Gen=="F3",]
  F3 <- F3[,-1]
  alleleSum <- as.matrix(colSums(F3))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F3Frequencies[,i] <- MajAlleleFreq
}
F3Frequencies <-rowMeans(F3Frequencies)
C1MAF[,4] <-F3Frequencies

for (i in 1:nreps) {
  F4 <- genotypes[[i]]
  F4 <- F4[F4$Gen=="F4",]
  F4 <- F4[,-1]
  alleleSum <- as.matrix(colSums(F4))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F4Frequencies[,i] <- MajAlleleFreq
}
F4Frequencies <-rowMeans(F4Frequencies)
C1MAF[,5] <-F4Frequencies

for (i in 1:nreps) {
  F5 <- genotypes[[i]]
  F5 <- F5[F5$Gen=="F5",]
  F5 <- F5[,-1]
  alleleSum <- as.matrix(colSums(F5))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F5Frequencies[,i] <- MajAlleleFreq
}
F5Frequencies <-rowMeans(F5Frequencies)
C1MAF[,6] <-F5Frequencies

for (i in 1:nreps) {
  PYT <- genotypes[[i]]
  PYT <- PYT[PYT$Gen=="PYT",]
  PYT <- PYT[,-1]
  alleleSum <- as.matrix(colSums(PYT))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  PYTFrequencies[,i] <- MajAlleleFreq
}
PYTFrequencies <-rowMeans(PYTFrequencies)
C1MAF[,7] <-PYTFrequencies

for (i in 1:nreps) {
  AYT <- genotypes[[i]]
  AYT <- AYT[AYT$Gen=="AYT",]
  AYT <- AYT[,-1]
  alleleSum <- as.matrix(colSums(AYT))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  AYTFrequencies[,i] <- MajAlleleFreq
}
AYTFrequencies <-rowMeans(AYTFrequencies)
C1MAF[,8] <-AYTFrequencies

for (i in 1:nreps) {
  Variety <- genotypes[[i]]
  Variety <- Variety[Variety$Gen=="Variety",]
  Variety <- Variety[,-1]
  alleleSum <- as.matrix(colSums(Variety))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  VarietyFrequencies[,i] <- MajAlleleFreq
  
}
VarietyFrequencies <-rowMeans(VarietyFrequencies)
C1MAF[,9] <-VarietyFrequencies

colnames(C1MAF) <- c("Base","F1","F2","F3","F4","F5","PYT","AYT","Variety")

###

genotypes <- readRDS("1C2rrblup_rd_alleles_snp_yield.rds")

C2MAF <- matrix(nrow=3547, ncol=9)

parentFrequencies <- matrix(nrow=3547,ncol=14)
F1Frequencies <- matrix(nrow=3547,ncol=14)
F2Frequencies <- matrix(nrow=3547,ncol=14)
F3Frequencies <- matrix(nrow=3547,ncol=14)
F4Frequencies <- matrix(nrow=3547,ncol=14)
F5Frequencies <- matrix(nrow=3547,ncol=14)
PYTFrequencies <- matrix(nrow=3547,ncol=14)
AYTFrequencies <- matrix(nrow=3547,ncol=14)
VarietyFrequencies <- matrix(nrow=3547,ncol=14)


nreps=14

for (i in 1:nreps) {
  parents <- genotypes[[i]]
  parents <- parents[parents$Gen=="NP",]
  parents <- parents[,-1]
  alleleSum <- as.matrix(colSums(parents))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  
  parentFrequencies[,i] <- MajAlleleFreq
}
parentFrequencies <-rowMeans(parentFrequencies)
C2MAF[,1] <-parentFrequencies



for (i in 1:nreps) {
  F1 <- genotypes[[i]]
  F1 <- F1[F1$Gen=="F1",]
  F1 <- F1[,-1]
  alleleSum <- as.matrix(colSums(F1))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F1Frequencies[,i] <- MajAlleleFreq
}

F1Frequencies <-rowMeans(F1Frequencies)
C2MAF[,2] <-F1Frequencies

for (i in 1:nreps) {
  F2 <- genotypes[[i]]
  F2 <- F2[F2$Gen=="F2",]
  F2 <- F2[,-1]
  alleleSum <- as.matrix(colSums(F2))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F2Frequencies[,i] <- MajAlleleFreq
}
F2Frequencies <-rowMeans(F2Frequencies)
C2MAF[,3] <-F2Frequencies

for (i in 1:nreps) {
  F3 <- genotypes[[i]]
  F3 <- F3[F3$Gen=="F3",]
  F3 <- F3[,-1]
  alleleSum <- as.matrix(colSums(F3))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F3Frequencies[,i] <- MajAlleleFreq
}
F3Frequencies <-rowMeans(F3Frequencies)
C2MAF[,4] <-F3Frequencies

for (i in 1:nreps) {
  F4 <- genotypes[[i]]
  F4 <- F4[F4$Gen=="F4",]
  F4 <- F4[,-1]
  alleleSum <- as.matrix(colSums(F4))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F4Frequencies[,i] <- MajAlleleFreq
}
F4Frequencies <-rowMeans(F4Frequencies)
C2MAF[,5] <-F4Frequencies

for (i in 1:nreps) {
  F5 <- genotypes[[i]]
  F5 <- F5[F5$Gen=="F5",]
  F5 <- F5[,-1]
  alleleSum <- as.matrix(colSums(F5))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F5Frequencies[,i] <- MajAlleleFreq
}
F5Frequencies <-rowMeans(F5Frequencies)
C2MAF[,6] <-F5Frequencies

for (i in 1:nreps) {
  PYT <- genotypes[[i]]
  PYT <- PYT[PYT$Gen=="PYT",]
  PYT <- PYT[,-1]
  alleleSum <- as.matrix(colSums(PYT))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  PYTFrequencies[,i] <- MajAlleleFreq
}
PYTFrequencies <-rowMeans(PYTFrequencies)
C2MAF[,7] <-PYTFrequencies

for (i in 1:nreps) {
  AYT <- genotypes[[i]]
  AYT <- AYT[AYT$Gen=="AYT",]
  AYT <- AYT[,-1]
  alleleSum <- as.matrix(colSums(AYT))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  AYTFrequencies[,i] <- MajAlleleFreq
}
AYTFrequencies <-rowMeans(AYTFrequencies)
C2MAF[,8] <-AYTFrequencies

for (i in 1:nreps) {
  Variety <- genotypes[[i]]
  Variety <- Variety[Variety$Gen=="Variety",]
  Variety <- Variety[,-1]
  alleleSum <- as.matrix(colSums(Variety))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  VarietyFrequencies[,i] <- MajAlleleFreq
  
}
VarietyFrequencies <-rowMeans(VarietyFrequencies)
C2MAF[,9] <-VarietyFrequencies

colnames(C2MAF) <- c("Base","F1","F2","F3","F4","F5","PYT","AYT","Variety")


###
genotypes <- readRDS("1C3rrblup_rd_alleles_snp_yield.rds")

C3MAF <- matrix(nrow=3547, ncol=9)

parentFrequencies <- matrix(nrow=3547,ncol=14)
F1Frequencies <- matrix(nrow=3547,ncol=14)
F2Frequencies <- matrix(nrow=3547,ncol=14)
F3Frequencies <- matrix(nrow=3547,ncol=14)
F4Frequencies <- matrix(nrow=3547,ncol=14)
F5Frequencies <- matrix(nrow=3547,ncol=14)
PYTFrequencies <- matrix(nrow=3547,ncol=14)
AYTFrequencies <- matrix(nrow=3547,ncol=14)
VarietyFrequencies <- matrix(nrow=3547,ncol=14)


nreps=14

for (i in 1:nreps) {
  parents <- genotypes[[i]]
  parents <- parents[parents$Gen=="NP",]
  parents <- parents[,-1]
  alleleSum <- as.matrix(colSums(parents))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  
  parentFrequencies[,i] <- MajAlleleFreq
}
parentFrequencies <-rowMeans(parentFrequencies)
C3MAF[,1] <-parentFrequencies



for (i in 1:nreps) {
  F1 <- genotypes[[i]]
  F1 <- F1[F1$Gen=="F1",]
  F1 <- F1[,-1]
  alleleSum <- as.matrix(colSums(F1))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F1Frequencies[,i] <- MajAlleleFreq
}

F1Frequencies <-rowMeans(F1Frequencies)
C3MAF[,2] <-F1Frequencies

for (i in 1:nreps) {
  F2 <- genotypes[[i]]
  F2 <- F2[F2$Gen=="F2",]
  F2 <- F2[,-1]
  alleleSum <- as.matrix(colSums(F2))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F2Frequencies[,i] <- MajAlleleFreq
}
F2Frequencies <-rowMeans(F2Frequencies)
C3MAF[,3] <-F2Frequencies

for (i in 1:nreps) {
  F3 <- genotypes[[i]]
  F3 <- F3[F3$Gen=="F3",]
  F3 <- F3[,-1]
  alleleSum <- as.matrix(colSums(F3))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F3Frequencies[,i] <- MajAlleleFreq
}
F3Frequencies <-rowMeans(F3Frequencies)
C3MAF[,4] <-F3Frequencies

for (i in 1:nreps) {
  F4 <- genotypes[[i]]
  F4 <- F4[F4$Gen=="F4",]
  F4 <- F4[,-1]
  alleleSum <- as.matrix(colSums(F4))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F4Frequencies[,i] <- MajAlleleFreq
}
F4Frequencies <-rowMeans(F4Frequencies)
C3MAF[,5] <-F4Frequencies

for (i in 1:nreps) {
  F5 <- genotypes[[i]]
  F5 <- F5[F5$Gen=="F5",]
  F5 <- F5[,-1]
  alleleSum <- as.matrix(colSums(F5))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  F5Frequencies[,i] <- MajAlleleFreq
}
F5Frequencies <-rowMeans(F5Frequencies)
C3MAF[,6] <-F5Frequencies

for (i in 1:nreps) {
  PYT <- genotypes[[i]]
  PYT <- PYT[PYT$Gen=="PYT",]
  PYT <- PYT[,-1]
  alleleSum <- as.matrix(colSums(PYT))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  PYTFrequencies[,i] <- MajAlleleFreq
}
PYTFrequencies <-rowMeans(PYTFrequencies)
C3MAF[,7] <-PYTFrequencies

for (i in 1:nreps) {
  AYT <- genotypes[[i]]
  AYT <- AYT[AYT$Gen=="AYT",]
  AYT <- AYT[,-1]
  alleleSum <- as.matrix(colSums(AYT))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  AYTFrequencies[,i] <- MajAlleleFreq
}
AYTFrequencies <-rowMeans(AYTFrequencies)
C3MAF[,8] <-AYTFrequencies

for (i in 1:nreps) {
  Variety <- genotypes[[i]]
  Variety <- Variety[Variety$Gen=="Variety",]
  Variety <- Variety[,-1]
  alleleSum <- as.matrix(colSums(Variety))
  alleleFreq <- as.data.frame(alleleSum/nrow(parents))
  
  MAJMIN <- alleleSum
  MAJMIN[MAJMIN < (nrow(parents)/2)] <- 0
  MAJMIN[MAJMIN > (nrow(parents)/2)] <- 1
  
  majAllele <- cbind(MAJMIN, alleleFreq)
  x <- as.vector(majAllele[,2])
  
  MajAlleleFreq <-- as.matrix(ifelse(x>0.5, x-1,x*(-1)))
  VarietyFrequencies[,i] <- MajAlleleFreq
  
}
VarietyFrequencies <-rowMeans(VarietyFrequencies)
C3MAF[,9] <-VarietyFrequencies

colnames(C3MAF) <- c("Base","F1","F2","F3","F4","F5","PYT","AYT","Variety")

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

alldata <- rbind(BOC,EOC)

library(overlapping)
library(lattice)

x <- list(C1=BOC$freq, C3=EOC$freq)

overlap(x, plot=TRUE, type="2")
