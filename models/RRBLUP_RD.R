

trainIndex <- as.matrix(sample(1:nrow(M)), 0.75*(nrow(M)))

phenoTrain <- y[trainIndex,]
genoTrain <- M[trainIndex,]

GM=tcrossprod(genoTrain)/dim(genoTrain)

BV <- phenoTrain

EBVans <-mixed.solve(BV, Z=GM, K=NULL, SE=FALSE, return.Hinv=FALSE)
markerEffects <- EBVans$u
markerEffects <- as.vector(markerEffects)
