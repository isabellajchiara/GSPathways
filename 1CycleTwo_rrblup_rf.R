
gvMatC2 <- matrix(nrow=10, ncol=1)
corMatC2 <- matrix(nrow=7, ncol=1)
varMatC2 <- matrix(nrow=10, ncol=1)


F1 = makeCross(newCycleSelections, crossPlan = cross, nProgeny = 5)
varMatC2[3,] = varG(F1)
gvMatC2[3,] <- mean(gv(F1))


allelesMatF1 <- pullSegSiteHaplo(F1)
Gen <- as.data.frame(rep("F1", times=nInd(F1)))
colnames(Gen) <- "Gen"
allelesMatF1 <- cbind(Gen, allelesMatF1)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny = 30) 
varMatC2[4,] = varG(F2)
gvMatC2[4,] <- mean(gv(F2))

allelesMatF2 <- pullSegSiteHaplo(F2)
Gen <- as.data.frame(rep("F2", times=nInd(F2)))
colnames(Gen) <- "Gen"
allelesMatF2 <- cbind(Gen, allelesMatF2)

source("rrblup_F2data.R")

##set EBV using BLUP model##
M_F2 <-pullSegSiteGeno(F2)
G_F2 = M_F2-1
EBVF2 <- G_F2 %*% markerEffects

F2@ebv <- as.matrix(EBVF2)
corMatC2[2] = cor(bv(F2), ebv(F2))


rm(newCycleSelections)
source("SelectParentsF2.R")

## select top individuals from F2 bulk  to form F3 ##

TopFamF2 = selectFam(F2, 10, use="ebv", top=TRUE) 
SelectionsF2 = selectWithinFam(TopFamF2, 100, use="ebv", top=TRUE)

F3 = self(SelectionsF2)
F3 = setPheno(F3)
varMatC2[5,] = varG(F3)
gvMatC2[5,] <- mean(gv(F3))

allelesMatF3 <- pullSegSiteHaplo(F3)
Gen <- as.data.frame(rep("F3", times=nInd(F3)))
colnames(Gen) <- "Gen"
allelesMatF3 <- cbind(Gen, allelesMatF3)



##set EBV using BLUP model##
M_F3 <-pullSegSiteGeno(F3)
G_F3 = M_F3-1
EBVF3 <- G_F3 %*% markerEffects

F3@ebv <- as.matrix(EBVF3)
corMatC2[3,] = cor(bv(F3),ebv(F3))

##select top within familiy from F3 to form F4 ##

TopFamF3 = selectFam(F3,5,use="pheno", top=TRUE) 
SelectionsF3 = selectWithinFam(TopFamF3, 50, use="pheno", top=TRUE)

F4 = self(SelectionsF3)
F4 = setPheno(F4)
varMatC2[6,] = varG(F4)
gvMatC2[6,] <- mean(gv(F4))

allelesMatF4 <- pullSegSiteHaplo(F4)
Gen <- as.data.frame(rep("F4", times=nInd(F4)))
colnames(Gen) <- "Gen"
allelesMatF4 <- cbind(Gen, allelesMatF4)


##set EBV using BLUP model##
M_F4 <-pullSegSiteGeno(F4)
G_F4 = M_F4-1
EBVF4 <- G_F4 %*% markerEffects

F4@ebv <- as.matrix(EBVF4)
corMatC2[4,] = cor(bv(F4),ebv(F4))

## select top families from F4 to form F5 ##

SelectionsF4 = selectFam(F4, 4, use="ebv")
F5 = self(SelectionsF4)
varMatC2[7,]= varG(F5)
gvMatC2[7,] <- mean(gv(F5))

allelesMatF5 <- pullSegSiteHaplo(F5)
Gen <- as.data.frame(rep("F5", times=nInd(F5)))
colnames(Gen) <- "Gen"
allelesMatF5 <- cbind(Gen, allelesMatF5)


#continue pipeline

##set EBV using BLUP model##
M_F5 <-pullSegSiteGeno(F5)
G_F5 = M_F5-1
EBVF5 <- G_F5 %*% markerEffects2

F5@ebv <- as.matrix(EBVF5)
corMatC2[5,] = cor(bv(F5),ebv(F5))

## select top F5 families for preliminary yield trial ##

SelectionsF5 = selectFam(F5, 3, use="ebv") 
PYT = self(SelectionsF5)
varMatC2[8,] = varG(PYT)
gvMatC2[8,] <- mean(gv(PYT))

allelesMatPYT <- pullSegSiteHaplo(PYT)
Gen <- as.data.frame(rep("PYT", times=nInd(PYT)))
colnames(Gen) <- "Gen"
allelesMatPYT <- cbind(Gen, allelesMatPYT)



##set EBV using BLUP model##
M_PYT <-pullSegSiteGeno(PYT)
G_PYT = M_PYT-1
EBVPYT <- G_PYT %*% markerEffects2

PYT@ebv <- as.matrix(EBVPYT)
corMatC2[6,] = cor(bv(PYT),ebv(PYT))

## select top families from PYT for AYT ##

SelectionsPYT = selectFam(PYT,  1, use="ebv", reps=5, top=TRUE) 
AYT = self(SelectionsPYT)
varMatC2[9,] = varG(AYT)
gvMatC2[9,] <- mean(gv(AYT))

allelesMatAYT <- pullSegSiteHaplo(AYT)
Gen <- as.data.frame(rep("AYT", times=nInd(AYT)))
colnames(Gen) <- "Gen"
allelesMatAYT <- cbind(Gen, allelesMatAYT)



##set EBV using BLUP model##
M_AYT <-pullSegSiteGeno(AYT)
G_AYT = M_AYT-1
EBVAYT <- G_AYT %*% markerEffects2

AYT@ebv <- as.matrix(EBVAYT)
corMatC2[7,] = cor(bv(AYT),ebv(AYT))

## select top plants to form variety ##
VarietySel = selectInd(AYT, 1, use="ebv")
Variety = self(VarietySel)
varMatC2[10,] = varG(Variety)
gvMatC2[10,] <- mean(gv(Variety))

allelesMatVar <- pullSegSiteHaplo(Variety)
Gen <- as.data.frame(rep("Variety", times=nInd(Variety)))
colnames(Gen) <- "Gen"
allelesMatVar <- cbind(Gen, allelesMatVar)

allelesMatC2 <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)


###collect bvs and ebvs###

bvebv <- cbind(bv(newParents), ebv(newParents))
Gen <- as.data.frame(rep("NP", times=nInd(newParents)))
bvebv <- cbind(Gen, bvebv)

bvebv1 <- cbind(bv(F2), ebv(F2))
Gen <- as.data.frame(rep("F2", times=nInd(F2)))
bvebv1 <- cbind(Gen, bvebv1)

bvebv2 <- cbind(bv(F3), ebv(F3))
Gen <- as.data.frame(rep("F3", times=nInd(F3)))
bvebv2 <- cbind(Gen, bvebv2)

bvebv3 <- cbind(bv(F4), ebv(F4))
Gen <- as.data.frame(rep("F4", times=nInd(F4)))
bvebv3 <- cbind(Gen, bvebv3)

bvebv4 <- cbind(bv(F5), ebv(F5))
Gen <- as.data.frame(rep("F5", times=nInd(F5)))
bvebv4 <- cbind(Gen, bvebv4)

bvebv5 <- cbind(bv(PYT), ebv(PYT))
Gen <- as.data.frame(rep("PYT", times=nInd(PYT)))
bvebv5 <- cbind(Gen, bvebv5)

bvebv6 <- cbind(bv(AYT), ebv(AYT))
Gen <- as.data.frame(rep("AYT", times=nInd(AYT)))
bvebv6 <- cbind(Gen, bvebv6)


bv_ebvC2 <- rbind(bvebv,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6)
colnames(bv_ebv) <- c("Gen","bv","ebv")




#write files - naming convention: "model_trainingSet_descriptor_populationType_trait.csv"



source("1CycleThree_rrblup.R")
