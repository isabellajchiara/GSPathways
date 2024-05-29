
# in Terminal - salloc session

library(dplyr)
library(cluster)
library(factoextra)

cycle = 3
genList = c("PYT")
genoList = list()

for (i in 1:nreps){
    
    filename = paste("C",cycle,"_ann_random_trainAtF5_trainWithF2_F2Parents_alleles_snp_yield.rds", sep="")
    genotypes <- readRDS(filename) #read in each cycle's SNP data
    repMat = genotypes[[i]]
    geno = repMat[repMat$Gen==gen,]
    M = geno[,-1]
    genoList[i] = M
    }

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



