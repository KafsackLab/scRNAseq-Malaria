#This code generates the seurat object for all ap2g-dd cells, particularly for figure 2 and others.
setwd("") #Set to the folder where all scripts and files are found
lapply(list("Seurat","parallel","doParallel","ggplot2","magrittr","reshape2","pheatmap","plyr","dplyr","scales","splines","zoo","RColorBrewer","reldist","entropy"),require,character.only = TRUE)
source("InternalFunctions.R") #load our DropSeq Functions
set.seed(123)
ggcols<-c(gg_color_hue(16)[1:11], "magenta1","deeppink1","pink","hotpink","deeppink2")
GeneInfo<-read.csv("GeneInfoTable_expanded_May2017_PlasmoDB32.csv",stringsAsFactors = F);
load("All_merged_forWeb.RData") # made by RawMatrices_to_dataFrame.R

# cellID holds the meta date of the cells in all_merge
# filtering cells for AP2-G Induced only
cellID<-cellID[cellID$strain!="NF54",] #remove the NF54 cells
cellID<-cellID[cellID$tp %in% c(30,36,42,"GC2"),] #keep all time points
#cellID<-cellID[cellID$treat=="OFF",] #to remove untreated cells
data_merge<-all_merge[,cellID$newID]; rm(all_merge);

noPolyA<-readLines("notPolyA_PlasmoDB32.txt")
data_merge<-data_merge[!(rownames(data_merge) %in% noPolyA),]

### Following the Seurat version 1.4 strategy (will not work with version 2)
### http://satijalab.org/seurat/
PfDS <- new("seurat", raw.data = as.matrix(data_merge)) #new empty Seurat object, genes as columns, dont include sampletype 
PfDS <- Setup(PfDS, #data
							 min.cells = 3, # Gene filter: Remove gene expressing fewer than this number of cells
							 min.genes = 1, # Cell filter: Remove cells expressing fewer than this number of genes
							 do.logNormalize = T, # convert to frequency, then  take natural log
							 total.expr = 1e4, # Scale to frequency to this number (multiply by 10,000)
							 project = "PfDS" # doesn't matter
)
# Cell Filtering by number of UMIs
plot(density(log(colSums(PfDS@raw.data),10)))
PfDS <- SubsetData(PfDS, subset.name = "nUMI", accept.high = 5000, accept.low = 300)

print(paste("Filter low nUMI cells:",dim(PfDS@scale.data)[2],"of",dim(PfDS@raw.data)[2],"remaining"))
cellIDfiltered<-cellID[match(names(PfDS@ident),cellID$newID),]

PfDS@data.info$sample<-cellIDfiltered$sample
PfDS@data.info$strain<-cellIDfiltered$strain
PfDS@data.info$exp<-cellIDfiltered$exp
PfDS@data.info$tp<-cellIDfiltered$tp
PfDS@data.info$treat<-cellIDfiltered$treat
PfDS@data.info$strain.treat<-paste(cellIDfiltered$strain,cellIDfiltered$treat,sep=".")
PfDS@data.info$strain.tp.treat<-paste(cellIDfiltered$strain,cellIDfiltered$tp,cellIDfiltered$treat,sep=".")

# Violin Plot
ggplot(PfDS@data.info,aes(x=sample,y=nUMI,group=sample, fill=tp)) + 
	geom_violin()+
	geom_boxplot(width=0.25, fill="white") +
	xlab("")+scale_y_log10(breaks=c(100,200,500,1000,2000,5000))+theme(axis.text.x  = element_text(angle = 45, hjust = 1))
#ggsave("TranscriptsViolin.pdf")

ggplot(PfDS@data.info,aes(x=sample,y=nGene,group=sample, fill=tp)) + 
	geom_violin()+
	geom_boxplot(width=0.25, fill="white") +
	xlab("")+scale_y_log10(breaks=c(50,100,200,500,1000,2000,5000))+theme(axis.text.x  = element_text(angle = 45, hjust = 1))

PfDS <- RegressOut(PfDS, latent.vars = c("nUMI"))

PfDS <- MeanVarPlot(	PfDS ,
                     fxn.x = expMean, 		# creates bins
                     fxn.y = logVarDivMean,  # how to choose most dispersed genes in each bin
                     x.low.cutoff = 0.03, 	# lower normalized log expression limit
                     x.high.cutoff = 10, 		# upper normalized log expression limit, usally set to > max
                     y.cutoff = 0.9, 		# minimal dispersion cutoff
                     do.contour = F,
                     do.plot = T
)

PfDS <- PCA(PfDS, pc.genes = setdiff(PfDS@var.genes,"PF3D7_1222600"), do.print = TRUE, pcs.print = 5, genes.print = 5)# run PCA for these gene, excluding ap2g
PfDS <- ProjectPCA(PfDS) # see if we excluded genes that correlate strongly with identified PCs, then include those gene
PfDS<-FindClusters(PfDS, pc.use = 1:26, resolution = 1.4, print.output = 0) #Use with PCA without AP2G for clustering
PfDS <- PCA(PfDS, pc.genes = PfDS@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)# run PCA for all genes
PfDS <- ProjectPCA(PfDS) # see if we excluded genes that correlate strongly with identified PCs, then include those gene
PfDS<-RunTSNE(PfDS, dims.use = 1:26, do.fast = T, perplexity=100, max_iter=5000,dim_embed = 3)

###Uncomment and adapt the following lines to reorder clusters
# clusterOrder<-data.frame(old=c(0,2,3,5,11,7,12,13,9,10,6,14,1,4,15,8)+1,new=1:16)
# newident<-factor(clusterOrder$new[match(as.numeric(PfDS@ident),clusterOrder$old)])
# names(newident)<-names(PfDS@ident)
# PfDS@ident<-newident;rm(newident)
# PfDS@data.info$cluster<-PfDS@ident

###transition from 3D tSNE to 2D tSNE
tSNE_3D_coordinates=PfDS@tsne.rot
pmat=scatter3D(x=PfDS@tsne.rot[,1],y=PfDS@tsne.rot[,3],z=-PfDS@tsne.rot[,2],col = ggcols[as.numeric(PfDS@ident)],colvar = NULL,theta=295,phi=2,pch=19,alpha=0.4,bty="b")
project2D=t(ldply(trans3D(x=PfDS@tsne.rot[,1],y=PfDS@tsne.rot[,3],z=-PfDS@tsne.rot[,2],pmat),.id = NULL))
rownames(project2D)=rownames(PfDS@tsne.rot);colnames(project2D)=colnames(PfDS@tsne.rot)[1:2]
plot(project2D,col=ggcols[as.numeric(PfDS@ident)])
scatter3D(x=tSNE_3D_coordinates[,1],y=tSNE_3D_coordinates[,3],z=-tSNE_3D_coordinates[,2],col = ggcols[as.numeric(PfDS@ident)],colvar = NULL,theta=295,phi=2,pch=19,alpha=0.4,bty="b")
PfDS@tsne.rot=data.frame(project2D)
clusterOrder<-data.frame(old=c(0,2,3,5,11,7,12,13,9,10,6,14,1,4,15,8)+1,new=1:16)
newident<-factor(clusterOrder$new[match(as.numeric(PfDS@ident),clusterOrder$old)])
names(newident)<-names(PfDS@ident)
PfDS@ident<-newident;rm(newident)
PfDS@data.info$cluster<-PfDS@ident

TSNEPlot(PfDS,do.label = T)

###############################################################################################################################################
# timepoint assignment using RNA Seq
fpkm<-read.csv("Bartfai_RNASeqTC.pct",sep="\t",stringsAsFactors = F,row.names=1)
colnames(fpkm)<-c("tp6","tp11","tp16","tp21","tp26","tp31","tp36","tp46");fpkm<-fpkm[,c(2:8,1)];
windowRNA<-log(fpkm)

PfDS.cor.all.rs<-timepoint.cor(CorMe=PfDS@scale.data,CorWith=windowRNA)
PfDS@data.info$assignedTP.rs<-names(fpkm)[apply(PfDS.cor.all.rs,1,which.max)]
PfDS@data.info$maxTPcor.rs<-apply(PfDS.cor.all.rs,1,max)

###############################################################################################################################################
save(PfDS,PfDS.cor.all.rs,tSNE_3D_coordinates,file="ap2g-dd_seurat.RData")

