if (args$trainingData == "F2") 
    M = as.data.frame(pullSegSiteGeno(F2))
    y = as.data.frame(pheno(F2))

if (args$trainingData == "F5") 
    M = as.data.frame(pullSegSiteGeno(F5))
    y = as.data.frame(pheno(F5)) 

if (args$trainingData == "F2_and_F5") 
    F2M = as.data.frame(pullSegSiteGeno(F2))
    F2y = as.data.frame(pheno(F2))
    F5M = as.data.frame(pullSegSiteGeno(F5))
    F5y = as.data.frame(pheno(F5))
    M = rbind(F2M,F5M)
    y = rbind(F2y,F5y)

trainIndex <- as.matrix(sample(1:nrow(M)), 0.75*(nrow(M)))

phenoTrain <- y[trainIndex,]
genoTrain <- M[trainIndex,]

BV <- phenoTrain

EBVans <-mixed.solve(BV, Z=genoTrain, K=NULL, SE=FALSE, return.Hinv=FALSE)
markerEffects <- EBVans$u
markerEffects <- as.vector(markerEffects)
