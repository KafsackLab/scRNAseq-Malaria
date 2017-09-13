### This code demonstrates how various calculations were performed using monocle, through the example of pseudo-time ordering ap2g-dd treated cells 
setwd("") #Set to the folder where all scripts and files are found
lapply(list("monocle","Seurat","parallel","doParallel","ggplot2","magrittr","reshape2","pheatmap","plyr","dplyr","scales","splines","zoo","RColorBrewer","reldist","entropy","grid","gridExtra","plotrix","stringr","png","igraph"),require,character.only = TRUE)
source("InternalFunctions.R") #load our DropSeq Functions
set.seed(123)
GeneInfo<-read.csv("GeneInfoTable_expanded_May2017_PlasmoDB32.csv",stringsAsFactors = F);
ggcols<-c(gg_color_hue(16)[1:11], "magenta1","deeppink1","deeppink2","hotpink","pink")

##monocle - pick what group of cells to order (ap2gdd treat only, nf54 only, both, etc) 
#for example, ap2gdd treated cells only
load("ap2g-dd_seurat.RData")  #see "Seurat_setup_code.R"

set.seed(10027)
subset=SubsetData(PfDS, cells.use = rownames(PfDS@data.info[PfDS@data.info$tp!="GC2",]))
subset=SubsetData(subset, cells.use = rownames(subset@data.info[subset@data.info$cluster %in% 1:11,]))
subset=SubsetData(subset, cells.use = rownames(subset@data.info[subset@data.info$treat=="OFF",]))
pd=new("AnnotatedDataFrame", data = data.frame(row.names=names(subset@ident),cluster=subset@data.info$cluster,tp=subset@data.info$tp,treat=subset@data.info$treat,ap2g=subset@data["PF3D7_1222600",]>0))
HSMM <- newCellDataSet(as.matrix(subset@raw.data[,names(subset@ident)]),
                       phenoData = pd,
                       lowerDetectionLimit=1,
                       expressionFamily=negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- names(which(rowSums(subset@raw.data[,names(subset@ident)]!=0)>20))

### to order by var genes
disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.01 &
                           dispersion_empirical >= 1.2 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM, reverse=T)
plot_cell_trajectory(HSMM, color_by="cluster",show_branch_points = F,show_backbone = F,show_tree = F)
subset@data.info$pt=max(pData(HSMM)$Pseudotime)-pData(HSMM)$Pseudotime

save(subset,HSMM,file="ap2g-dd_treatedAsexuals_monocle.RData")

