# in Terminal - salloc session

library(dplyr)
library(cluster)
library(factoextra)

genList = c("PYT")
repList = list()
reps = 15

filename = paste("C3_ann_random_trainAtF2_trainWithF2_and_F5_F2Parents_alleles_snp_yield.rds")
genotypes <- readRDS(filename) #read in each cycle's SNP data

reps = 15
gen = "PYT"
repList = list()
for (rep in 1:reps){
  repMat = genotypes[[rep]]
  geno = repMat[repMat$Gen==gen,]
  M = geno[,-1]
  repList[[rep]] <- assign(paste0("M",rep),M)
}

ex = repList[[1]]
nInd = nrow(ex)
nLoci = ncol(ex)

indList = list()
for (n in 1:nInd){
  indMat= list()
  indList[[n]] = assign(paste0("ind",n),indMat)
}

for (rep in 1:length(repList)){
  for (x in 1:nInd){
    repMat = repList[[rep]] #pull on rep
    ind = as.matrix(repMat[x,]) #take x ind
    matrix = indList[[x]] #pull out the x ind matrix
    matrix[[rep]] = ind #assign rep values to the matrix
    indList[[x]] = matrix #replace matrix in IndList with new values
}
}

indListDF = list()
for (list in 1:length(indList)){
  ind = indList[[list]]
  indDF = as.data.frame(do.call(rbind,ind))
  indListDF[[list]] = indDF
}

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

finalModes = list()
for (DF in 1:length(indListDF)){
  indMode = matrix(nrow=1,ncol=nLoci)
  for (loc in 1:nLoci){
    ind = as.data.frame(indListDF[[DF]]) #take one individual, all reps(rows), all loci(columns)
    locus = as.data.frame((ind[,loc])) # take one locus
    mode = Modes(locus) # find mode state
    if (nrow(mode) > 1){
      mode = mode[1,1]
    }
    mode = as.matrix(mode)
    indMode[1,loc] = mode #add mode to new output and repeat for all loci
  }
  finalModes[[DF]] = indMode #collect modes for each ind 
}

M = do.call(rbind,finalModes)
      
#prepare DF
newgeno <- M %>%  select(where(~ n_distinct(.) > 1))
newgeno = scale(newgeno)
colnames(newgeno) = NULL
rownames(newgeno) = NULL
saveRDS(newgeno,paste("newgenoC",cycle,"mat.RDS"))

#visualize distances on genotype matrix
distance <- get_dist(newgeno)
saveRDS(distance,paste("distanceC",cycle,"mat.RDS"))


#
#
#
#
# in RStudio Console 

library(dplyr)
library(cluster)
library(factoextra)

distance = readRDS("distanceC3mat.RDS")
fviz_dist(distance,show_labels = FALSE ,gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# optimize K
newgeno = readRDS("newgenoC3mat.RDS")
PCAgeno <- prcomp(newgeno, center=TRUE, scale=TRUE) ##take out categorical columns##
PCAselected = as.data.frame(-PCAgeno$x[,1:3])
silhouette <- fviz_nbclust(PCAselected, kmeans, method = 'silhouette')
kvalues <- silhouette$data ##largest value tells how many clusters are optimal ##
kvalues <- kvalues[order(-kvalues$y),]
k=as.numeric(kvalues[1,1])

if (length(genList) > 1){
  main = paste(genList[[1]],genList[[2]])
}else{
  main=genList}

# apply kmeans function and visualize
k2 <- kmeans(newgeno, centers = k, nstart = 25)
fviz_cluster(k2, data = newgeno,geom="point",ggtheme=theme_minimal(), ellipse = TRUE,
             ellipse.type = "convex",
             ellipse.level = 0.95,
             ellipse.alpha = 0.2,
             main = paste0(main," genotypes"))

