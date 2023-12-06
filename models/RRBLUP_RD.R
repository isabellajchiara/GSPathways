

phenoTrain <- as.matrix(y)
genoTrain <- as.matrix(M)

EBVans <-mixed.solve(phenoTrain, Z=genoTrain, K=NULL, SE=FALSE, return.Hinv=FALSE)
markerEffects <- EBVans$u
markerEffects <- as.vector(markerEffects)