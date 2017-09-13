#This code, similarly to "Seurat_setup_code.R", demonstrates utilization of the Seurat package to cluster and visualize cells. Specifically, we choose cells from various biological and technical replicates at the different stages of the parasite life-cycle to create figures 1b,c
setwd("") #Set to the folder where all scripts and files are found
lapply(list("Seurat","parallel","doParallel","ggplot2","magrittr","reshape2","pheatmap","plyr","dplyr","scales","splines","zoo","RColorBrewer"),require,character.only = TRUE)
source("InternalFunctions.R")
GeneInfo<-read.csv("GeneInfoTable_expanded_May2017_PlasmoDB32.csv",stringsAsFactors = F);
ggcols<-c(gg_color_hue(16)[1:11], "magenta1","deeppink1","deeppink2","hotpink","pink")
set.seed(123)
#########################################################################################################
load("All_merged_forWeb.RData")
# select 36h Trophs, day 2 GC, Day 9 GC

cellID<-cellID[cellID$sample %in% c("DCJ.E.18.BSD.A","DCJ.E.18.BSD.B","AP2G.B.36.OFF.B","AP2G.B.36.OFF.A","AP2G.B.GC2.ON.A","AP2G.A.GC2.ON.A","NF54.C.GC9.Ind.A"),]
cellID<-cellID[-sample(which(cellID$sample=="AP2G.A.GC2.ON.A"),1000),] #downsample AP2G.A.GC2.ON.A
apply(cellID[,-1:-2],2,table)
data_merge<-all_merge[rowSums(all_merge[,cellID$newID])>0,cellID$newID]


PfDS <- new("seurat", raw.data = data_merge) #new empty Seurat object, genes as columns, dont include sampletype 
PfDS <- Setup(	PfDS, #data
							 min.cells = 3, # Gene filter: Remove gene expressing fewer than this number of cells
							 min.genes = 1, # Cell filter: Remove cells expressing fewer than this number of genes
							 do.logNormalize = T, # convert to frequency, then  take natural log
							 total.expr = 1e4, # Scale to frequency to this number (multiply by 10,000)
							 project = "PfDS" # doesn't matter
)
# Cell Filtering by number of UMIs
plot(density(log10(colSums(PfDS@raw.data))))
PfDS <- SubsetData(PfDS, subset.name = "nUMI", accept.high = 5000, accept.low = 100)

print(paste("Filter low nUMI cells:",dim(PfDS@scale.data)[2],"of",dim(PfDS@raw.data)[2],"remaining"))
cellIDfiltered<-cellID[match(names(PfDS@ident),cellID$newID),]

PfDS@data.info$sample<-factor(cellIDfiltered$sample)
PfDS@data.info$strain<-cellIDfiltered$strain
PfDS@data.info$exp<-cellIDfiltered$exp
PfDS@data.info$tp<-cellIDfiltered$tp
PfDS@data.info$treat<-cellIDfiltered$treat
PfDS@data.info$AP2G<-cellIDfiltered$AP2G=="ON"
PfDS@data.info$GDV1<-cellIDfiltered$GDV1=="ON"
PfDS@data.info$strain.treat<-paste(cellIDfiltered$strain,cellIDfiltered$treat,sep=".")
PfDS@data.info$strain.tp.treat<-paste(cellIDfiltered$strain,cellIDfiltered$tp,cellIDfiltered$treat,sep=".")
PfDS@data.info$strain.tp<-paste(cellIDfiltered$strain,cellIDfiltered$tp)
PfDS@data.info$both<-(PfDS@data.info$GDV1 & PfDS@data.info$AP2G)

#VlnPlot(PfDS, c("nGene", "nUMI"), nCol = 2, group.by="rename", size.use=0)
multiplot(
ggplot(PfDS@data.info,aes(x=sample,y=nUMI,group=sample, fill=tp)) + 
	geom_violin()+
	geom_boxplot(width=0.25, fill="white") +
	xlab("")+scale_y_log10(breaks=c(100,200,500,1000,2000,5000))+theme(axis.text.x  = element_text(angle = 45, hjust = 1))+

,

ggplot(PfDS@data.info,aes(x=sample,y=nGene,group=sample, fill=tp)) + 
	geom_violin()+
	geom_boxplot(width=0.25, fill="white") +
	xlab("")+scale_y_log10(breaks=c(50,100,200,500,1000,2000,5000))+theme(axis.text.x  = element_text(angle = 45, hjust = 1))+
,cols=1)

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
length(PfDS@var.genes)

PfDS <- PCA(PfDS, pc.genes = PfDS@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)# run PCA for these gene
PfDS <- ProjectPCA(PfDS) # see if we excluded genes that correlate strongly with identified PCs, then include those gene
PfDS<-FindClusters(PfDS, pc.use = 1:17, resolution = 0.9, print.output = 0) #Use with PCA without AP2G for clustering
PfDS<-RunTSNE(PfDS, dims.use = 1:17, do.fast = T, perplexity=50, max_iter=3000,dim_embed = 2)

save(cellID,cellIDfiltered,PfDS,GeneInfo,data_merge,file="Figure1_seuratObject.RData")
