F1 = makeCross(newCycleSelections, crossPlan = cross, nProgeny = 5)
var1 <-varG(F1)

F2 = self(F1, nProgeny = 5) ##nProgeny = number of progeny per cross## 
var2 <-varG(F2)


##set EBV using BLUP model##
M_F2 <-pullSegSiteGeno(F2)
G_F2 = M_F2-1
EBVF2 <- G_F2 %*% markerEffects
F2@ebv <- as.matrix(EBVF2)

cor2 = cor(gv(F2), ebv(F2))

source("SelectParentsF2.R")


## select top families to form F3 ##

F3Sel = selectWithinFam(F2, 6, use="ebv", top=TRUE) 
F3 = self(F3Sel, nProgeny = 6)

##set EBV using BLUP model##
M_F3 <-pullSegSiteGeno(F3)
G_F3 = M_F3-1
EBVF3 <- G_F3 %*% markerEffects
F3@ebv <- as.matrix(EBVF3)


cor3 = cor(gv(F3),ebv(F3))
var3 = varG(F3)


##select top families from F3 to form F4 ##

F4Sel = selectWithinFam(F3, 10, use="ebv", top=TRUE) 
F4 = self(F4Sel)

##set EBV using BLUP model##
M_F4 <-pullSegSiteGeno(F4)
G_F4 = M_F4-1
EBVF4 <- G_F4 %*% markerEffects
F4@ebv <- as.matrix(EBVF4)

cor4 = cor(gv(F4),ebv(F4))
var4 = varG(F4)

F5Sel = selectFam(F4, 10, use="ebv")
F5 = self(F5Sel, nProgeny=3)


source("RF_RD_RetrainF2F5.R")

#continue pipeline

##set EBV using BLUP model##
M = as.data.frame(pullSegSiteGeno(F5))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
F5ebv <- as.numeric(predict(rf_fit2, M))
F5@ebv <- as.matrix(F5ebv)

cor5 = cor(gv(F5),ebv(F5))
var5 = varG(F5)


## select top F5 families for preliminary yield trial ##

PYTSel = selectFam(F5, 4, use="ebv") 
PYT = self(PYTSel)
var6 = varG(PYT)

##set EBV using BLUP model##
M = as.data.frame(pullSegSiteGeno(PYT))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
PYTebv <- as.numeric(predict(rf_fit2, M))
PYT@ebv <- as.matrix(PYTebv)

cor6 = cor(gv(PYT),ebv(PYT))

## select top families from PYT for AYT ##

AYTSel = selectFam(PYT,  3, use="ebv", reps=5, top=TRUE) 
AYT = self(AYTSel)
var7 = varG(AYT)

##set EBV using BLUP model##
M = as.data.frame(pullSegSiteGeno(AYT))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
AYTebv <- as.numeric(predict(rf_fit2, M))
AYT@ebv <- as.matrix(AYTebv)

cor7 = cor(gv(AYT),ebv(AYT))

## select top plants to form variety ##
VarietySel = selectInd(AYT, 1, use="ebv")
Variety = self(VarietySel)
var8 = varG(Variety)


## pull genetic value for each generation and write results ##

F1gv2 <- as.vector(mean(gv(F1)))
F2gv2 <- as.vector(mean(gv(F2)))
F3gv2 <- as.vector(mean(gv(F3)))
F4gv2 <- as.vector(mean(gv(F4)))
F5gv2 <- as.vector(mean(gv(F5)))
PYTgv2 <- as.vector(mean(gv(PYT)))
AYTgv2 <- as.vector(mean(gv(AYT)))
Varietygv2 <- as.vector(mean(gv(Variety)))

cycle2gv <- matrix(nrow=8, ncol=1)
cycle2gv[1,] <- F1gv2
cycle2gv[2,] <- F2gv2
cycle2gv[3,] <- F3gv2
cycle2gv[4,] <- F4gv2
cycle2gv[5,] <- F5gv2
cycle2gv[6,] <- PYTgv2
cycle2gv[7,] <- AYTgv2
cycle2gv[8,] <- Varietygv2

cycle2gvs <- as.matrix(cycle2gv)

###list correlations to view model performacne ##
corMat <- matrix(nrow=6, ncol=1)
corMat[1,] <- cor2
corMat[2,] <- cor3
corMat[3,] <- cor4
corMat[4,] <- cor5
corMat[5,] <- cor6
corMat[6,] <- cor7

cycle2cors <- as.matrix(corMat)


varMat <- matrix(nrow=7, ncol=1)
varMat[1,] <- var1
varMat[2,] <- var2
varMat[3,] <- var3
varMat[4,] <- var4
varMat[5,] <- var5
varMat[6,] <- var6
varMat[7,] <- var7


cycle2vars <- as.matrix(varMat)


source("1CycleThreeF2.R")
