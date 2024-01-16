setwd("/Users/justin_folders/Desktop/Isabella_McGill_Files/sim_2024Jan08_114041 AM")
ANNpheno = read.csv("C1_ann_random_trainAtF2_trainWithF2_F5Parents_pheno_snp_yield.csv")
ANNpheno = ANNpheno[,-1]
ANNebv = readRDS("C1_ann_random_trainAtF2_trainWithF2_F5Parents_bvebv_snp_yield.rds")
setwd("/Users/justin_folders/Desktop/Isabella_McGill_Files")
RRBLUPpheno = read.csv("C1_rrblup_random_trainAtF2_trainWithF2_F5Parents_pheno_snp_yield.csv")
RRBLUPpheno = RRBLUPpheno[,-1]
RRBLUPebv = readRDS("C1_rrblup_random_trainAtF2_trainWithF2_F5Parents_bvebv_snp_yield.rds")

ANNebvs = list()
x = 1
for (x in 1:length(ANNebv)){
ANNebv1 = ANNebv[[x]]
GenNames = as.data.frame(ANNebv1[,1])
ANNebv1 = ANNebv1[,-c(1:2)]
ANNebv1 = cbind(GenNames,ANNebv1)
ANNebvs[[x]] = ANNebv1
x=x+1
}
ANNebvData = do.call(rbind, ANNebvs)

###EBVS ARE COLLECTED

##NOW COLLECT CORRESPONDING PHENOS

z = 1 
for (z in 1:nrow(ANNpheno)){
  check = grep("F1C1", ANNpheno[z,])
  if (length(check) > 0 ){
    ANNpheno[z,1:ncol(ANNpheno)] = NA
    ANNpheno[z,1:ncol(ANNpheno)] = NA
    z=z+1
  }
}
ANNpheno = as.data.frame(ANNpheno[-c(1:60),]) #remove prev pyts
ANNpheno =na.omit(ANNpheno)
ANNpheno = ANNpheno[-nrow(ANNpheno),] #remove last variety row, no ebv

GenNames = as.data.frame(ANNpheno[,2])

toDelete <- seq(0, length(ANNpheno), 2)
ANNPHENO = ANNpheno[,-c(toDelete)]

rep=1
phenodata = list()
for (col in 1:ncol(ANNPHENO)){
  column = ANNPHENO[,col]
  data = cbind(GenNames,column)
  phenodata[[rep]] = data
  rep = rep+1
}

ANNphenoData = do.call(rbind, phenodata)

dim(ANNebvData)
dim(ANNphenoData)

#### NOW RRBLUP
RRBLUPebvs = list()
x = 1
for (x in 1:length(RRBLUPebv)){
  RRBLUPebv1 = RRBLUPebv[[x]]
  GenNames = as.data.frame(RRBLUPebv1[,1])
  RRBLUPebv1 = RRBLUPebv1[,-c(1:2)]
  RRBLUPebv1 = cbind(GenNames,RRBLUPebv1)
  RRBLUPebvs[[x]] = RRBLUPebv1
  x=x+1
}
RRBLUPebvData = do.call(rbind, RRBLUPebvs)

###EBVS ARE COLLECTED

##NOW COLLECT CORRESPONDING BVs

RRBLUPebvs = list()
x = 1
for (x in 1:length(RRBLUPebv)){
  RRBLUPebv1 = RRBLUPebv[[x]]
  GenNames = as.data.frame(RRBLUPebv1[,1])
  RRBLUPebv1 = RRBLUPebv1[,-c(1,3)]
  RRBLUPebv1 = cbind(GenNames,RRBLUPebv1)
  RRBLUPebvs[[x]] = RRBLUPebv1
  x=x+1
}
RRBLUPebvData = do.call(rbind, RRBLUPebvs)



ANNebvData$model = rep(c("ANN"), times=nrow(ANNebvData))
ANNphenoData$model = rep(c("ANN"), times=nrow(ANNphenoData))
ANNDATA = cbind(ANNebvData,ANNphenoData)
ANNDATA = ANNDATA[,-c(4,6)]
colnames(ANNDATA) = c("gen","estimated","model","true")

RRBLUPebvData$model = rep(c("RRBLUP"), times=nrow(RRBLUPebvData))
RRBLUPphenoData$model = rep(c("RRBLUP"), times=nrow(RRBLUPphenoData))
RRBLUPDATA = cbind(RRBLUPebvData,RRBLUPphenoData)
RRBLUPDATA = RRBLUPDATA[,-c(4,6)]
colnames(RRBLUPDATA) = c("gen","estimated","model","true")



scale(RRBLUPDATA$estimated, center = TRUE, scale = TRUE)
scale(RRBLUPDATA$true, center = TRUE, scale = TRUE)

scale(ANNDATA$estimated, center = TRUE, scale = TRUE)
scale(ANNDATA$true, center = TRUE, scale = TRUE)


allData = rbind(RRBLUPDATA,ANNDATA)


plot <- ggplot(allData[allData$model=="RRBLUP",], aes(x=estimated, y=true, color=gen)) +
  geom_point()

plot + labs(title="RRBLUP Performance")