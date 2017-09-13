#This code generate differentially expressed gene list by cluster, when comparing ap2g+ to ap2g- in treated ap2gdd cells of NF54 cells
setwd("") #Set to the folder where all scripts and files are found
lapply(list("monocle","plyr","dplyr","Seurat","parallel","doParallel","ggplot2","magrittr","reshape2","pheatmap","scales","splines","zoo","RColorBrewer","reldist","entropy","grid","gridExtra","plotrix","stringr","png","igraph","scatterpie"),require,character.only = TRUE)
source("InternalFunctions.R") #load our DropSeq Functions
set.seed(123)
GeneInfo<-read.csv("GeneInfoTable_expanded_May2017_PlasmoDB32.csv",stringsAsFactors = F);
ggcols<-c(gg_color_hue(16)[1:11], "magenta1","deeppink1","deeppink2","hotpink","pink")


####First, create and load an "NF54 only" seurat object (modifying relevant lines in "Seurat_setup_code.R")
clus=9 ##Determine the clusters that are clearly asexual 
PfDS@data.info$treat[PfDS@data.info$strain=="NF54"]="ON"
PfDS@data.info$ap2g.treat=paste(PfDS@data["PF3D7_1222600",]>0,PfDS@data.info$treat,sep="_")
subset<-SubsetData(PfDS,cells.use=rownames(PfDS@data.info[PfDS@data.info$ap2g.treat!="TRUE_OFF" & PfDS@data.info$ap2g.treat!="FALSE_OFF",]))
subset=SubsetData(subset, cells.use = rownames(subset@data.info[subset@data.info$tp!="GC2",]))

ONvsOFF_dgeByCluster<-InternalDGE(
  SeuratObj = subset,
  compare.by=subset@data.info$ap2g.treat=="TRUE_ON",
  minPct=0.1,
  min.diff.pct=0.175,
  ths=0.4,  # previously was 0.4
  cellmin=12,
  posOnly = F)

deGenesAllClusNF54<-ONvsOFF_dgeByCluster[[2]]
deGenesAllClusNF54<-deGenesAllClusNF54[deGenesAllClusNF54$cluster%in% 1:clus & deGenesAllClusNF54$avg_diff >= 0.7,]

deGenesAllClusNF54=cbind(cluster=deGenesAllClusNF54$cluster,getInfo(as.character(deGenesAllClusNF54$geneID)),deGenesAllClusNF54[,c("pct.1","pct.2","pctDiff","avg_diff","p_val")])
table(deGenesAllClusNF54$cluster); length(unique(deGenesAllClusNF54$geneID)); length(unique(deGenesAllClusNF54$geneID[duplicated(deGenesAllClusNF54$geneID)]))

deGenesAllClusNF54$pctChange<-deGenesAllClusNF54$pctDiff/apply(deGenesAllClusNF54[,c("pct.1","pct.2")],1,max)
deGenesAllClusNF54$avg_diff<-log(exp(deGenesAllClusNF54$avg_diff),2) #log2 difference
deGenesAllClusNF54$motifs<-getInfo(deGenesAllClusNF54$geneID,"ap2gMotifs_3k") #log2 difference
deGenesAllClusNF54=deGenesAllClusNF54[deGenesAllClusNF54$geneID!="PF3D7_1222600",]

##Repeat for ap2g-dd cells alone
load("ap2g-dd_seurat.RData")
clus=11
PfDS@data.info$ap2g.treat=paste(PfDS@data["PF3D7_1222600",]>0,PfDS@data.info$treat,sep="_")
subset<-SubsetData(PfDS,cells.use=rownames(PfDS@data.info[PfDS@data.info$ap2g.treat!="TRUE_OFF" & PfDS@data.info$ap2g.treat!="FALSE_OFF",]))
subset=SubsetData(subset, cells.use = rownames(subset@data.info[subset@data.info$tp!="GC2",]))

ONvsOFF_dgeByCluster<-InternalDGE(
  SeuratObj = subset,
  compare.by=subset@data.info$ap2g.treat=="TRUE_ON",
  minPct=0.1,
  min.diff.pct=0.175,
  ths=0.4,  # previously was 0.4
  cellmin=12,
  posOnly = F)

deGenesAllClusAP2G<-ONvsOFF_dgeByCluster[[2]]
deGenesAllClusAP2G<-deGenesAllClusAP2G[deGenesAllClusAP2G$cluster%in% 1:clus & deGenesAllClusAP2G$avg_diff >= 0.7,]

deGenesAllClusAP2G=cbind(cluster=deGenesAllClusAP2G$cluster,getInfo(as.character(deGenesAllClusAP2G$geneID)),deGenesAllClusAP2G[,c("pct.1","pct.2","pctDiff","avg_diff","p_val")])
table(deGenesAllClusAP2G$cluster); length(unique(deGenesAllClusAP2G$geneID)); length(unique(deGenesAllClusAP2G$geneID[duplicated(deGenesAllClusAP2G$geneID)]))

deGenesAllClusAP2G$pctChange<-deGenesAllClusAP2G$pctDiff/apply(deGenesAllClusAP2G[,c("pct.1","pct.2")],1,max)
deGenesAllClusAP2G$avg_diff<-log(exp(deGenesAllClusAP2G$avg_diff),2) #log2 difference
deGenesAllClusAP2G$motifs<-getInfo(deGenesAllClusAP2G$geneID,"ap2gMotifs_3k") #log2 difference
deGenesAllClusAP2G=deGenesAllClusAP2G[deGenesAllClusAP2G$geneID!="PF3D7_1222600",]

save(deGenesAllClusAP2G,deGenesAllClusNF54,file="DE_EachStrainSeparately.RData")


