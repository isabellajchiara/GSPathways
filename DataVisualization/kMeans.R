library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)
library(ggpubr)

cycle = 2
genList = c("F2","F5")
genoList = list()

for (gen in genList){
filename = paste("C", cycle, "_rrblup_random_trainAtF5_trainWithF2_F2Parents_alleles_snp_yield.rds", sep="")
genotypes <- readRDS(filename) #read in each cycle's SNP data
repMat = genotypes[[1]]
geno = repMat[repMat$Gen==gen,]
M = geno[,-1]

  if (gen == "F2"){
      M = M[1:240,]}

assign(paste0("M",gen),M)}

if (length(genList) > 1){
  M  = rbind(MF2,MF5)
}

#prepare DF
newgeno <- M %>%  select(where(~ n_distinct(.) > 1))
newgeno = scale(newgeno)
colnames(newgeno) =NULL
rownames(newgeno) = c(paste0("geno",1:nrow(newgeno)))

# optimize K
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



#visualize distances on genotype matrix
distance <- get_dist(newgeno)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
