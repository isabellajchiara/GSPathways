M <- as.data.frame(TrainingGeno)
y <- as.data.frame(TrainingPheno)

trainIndex <- as.matrix(sample(1:nInd(TP), 0.75*(nrow(M)))) 

phenoTrain <- y[trainIndex,]
genoTrain <- M[trainIndex,]

BV <- phenoTrain

EBVans <-mixed.solve(BV, Z=genoTrain, K=NULL, SE=FALSE, return.Hinv=FALSE)
markerEffects <- EBVans$u
markerEffects <- as.vector(markerEffects)
