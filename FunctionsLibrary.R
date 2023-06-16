# Defining trait parameters (AEG)
defineTraitAEG <- function(nQtl,mean,h2) {
  SP <<- SimParam$new(founderPop)
  SP$addTraitAEG(nQtl, mean=mean)
  SP$setVarE(h2=h2)
}

# Defining trait parameters (A)
defineTraitA <- function(nQtl,mean,h2) {
  SP <<- SimParam$new(founderPop)
  SP$addTraitA(nQtl, mean=mean)
  SP$setVarE(h2=h2)
}

# Selecting parents for the next cycle (pheno)
selectNewParents <- function(gen,nInd,criterion){
  selectInd(gen, nInd, use=criterion)
}

# Within Family Selections (pheno/ebv)
TopWithinFam <- function(gen,nFam,nIndPerFam,criterion){
  TopFam <-selectFam(gen,nFam, use=criterion, top=TRUE)
  Selections <<- selectWithinFam(TopFam, nIndPerFam,use=criterion, top=TRUE)
  self(Selections)
}

# Family Selections (pheno/ebv)
TopFamily <- function(gen,nFam,criterion){
  Top = selectFam(gen,nFam, use=criterion, top=TRUE)
  self(Top)
}

#Stratified Clusters

StratClusTRN <- function(y,M) { #y= matrix of training phenotypes M= matrix training genotypes
  
  newgeno <- M %>%  select(where(~ n_distinct(.) > 1))
  
  colnames(newgeno) =NULL
  
  PCAgeno <- prcomp(newgeno, center=TRUE, scale=TRUE) ##take out categorical columns##
  
  PCAselected = as.data.frame(-PCAgeno$x[,1:3])
  
  silhouette <- fviz_nbclust(PCAselected, kmeans, method = 'silhouette')
  kvalues <- silhouette$data ##largest value tells how many clusters are optimal ##
  kvalues <- kvalues[order(-kvalues$y),]
  
  k=as.numeric(kvalues[1,1])
  
  kmeans_geno = kmeans(PCAselected, centers = k, nstart = 50)
  clusters <- fviz_cluster(kmeans_geno, data = PCAselected)
  
  clusterData <- clusters$data
  
  clusterData <- clusterData[order(clusterData$cluster),]
  
  nclusters <- as.numeric(clusterData[as.numeric(nrow(clusterData)),as.numeric(ncol(clusterData))])
  
  datalist = vector("list", length = nclusters)
  
  for (x in 1:nclusters) {
    clustername <- paste0("cluster",x)
    clustername <- clusterData[clusterData$cluster==x,] 
    
    assign(paste0("cluster",x), clustername)
    
    trnname <- paste0("trn",x)
    trnname <- clustername[sample(0.75*nrow(clustername)),]
    datalist[[x]] <- trnname
    
  }
  
  TRN <- do.call(rbind, datalist)
  
  TRN <- TRN[,1]
  
  M <- as.data.frame(TrainingGeno)
  rownames(M) <- c(1:nrow(M))
  OptimGeno <<- as.matrix(M[TRN,])
  y <- as.data.frame(y)
  OptimPheno <<- y[TRN,]
  
}

# RRBLUP estimate ebvs
GetEBVrrblup <- function(gen){
  genMat <- pullSegSiteGeno(gen) 
  genMat <- genMat-1
  genMat %*% markerEffects
}

GetEBVrf <- function(gen){
  M = as.data.frame(pullSegSiteGeno(gen))
  colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
  as.numeric(predict(rf_fit, M))
}

# gets alleles matrix of genObj
getAllelesMat <- function(genObj, genName){
    allelesMat <- pullSegSiteHaplo(genObj)
    Gen <- as.data.frame(rep(genName, times=nInd(genObj)))
    colnames(Gen) <- "Gen"
    allelesMat <- cbind(Gen, allelesMat)
    allelesMat
}

# gets dataframe with bv and ebv
getBvEbv <- function(genObj, genName){
    bvebv <- cbind(bv(genObj), ebv(genObj))
    Gen <- as.data.frame(rep(genName, times=nInd(genObj)))
    bvebv <- cbind(Gen, bvebv)
    colnames(bvebv) <- c("Gen","bv","ebv")
    bvebv
}

## Functions to build dataframes

getAllGeneticValues <- function(geneticValues, lin1, lin2){
  geneticValues <- as.data.frame(geneticValues)
  colnames(geneticValues) <- 1:nReps
  gain <- as.data.frame(geneticValues[lin1,] - geneticValues[lin2,])
  colnames(gain) <- 1:nReps
  AllgeneticValues <- as.data.frame(rbind(geneticValues, gain))
  rownames(AllgeneticValues) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety","meanGV")
  colnames(AllgeneticValues) <- c(1:nReps)
  AllgeneticValues
}

getCorrelations <- function(correlations){
  correlations <- as.data.frame(correlations)
  rownames(correlations) <- c("NewParents","F2","F3","F4","F5","PYT","AYT")
  colnames(correlations) <- c(1:nReps)  
  correlations
}

getVariances <- function(variances){
  variances <- as.data.frame(variances)
  colnames(variances) <- c(1:nReps)
  rownames(variances) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT")
  variances
}
