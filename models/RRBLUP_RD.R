

phenoTrain <- y
genoTrain <- M

GM=tcrossprod(genoTrain)/dim(genoTrain)

BV <- phenoTrain

EBVans <-mixed.solve(BV, Z=GM, K=NULL, SE=FALSE, return.Hinv=FALSE)
markerEffects <- EBVans$u
markerEffects <- as.vector(markerEffects)
