#This file contains functions utilized throughout the other scripts for various purposes, including visualization, differential gene expression, etc
# allows grid plotting of multiple stored plots
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
	require(grid)
	plots <- c(list(...), plotlist)
	numPlots = length(plots)
	if (is.null(layout)) layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
	if (numPlots == 1) {print(plots[[1]])} 
	else {
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		for (i in 1:numPlots) {
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		}
	}
}
########################################################################################################################################

ParallelPcaRes<-function(seuratobj,PCAS,RES,cores=32,perplexity=50,max_iter = 3000,dim=2){
	registerDoParallel(cores = cores) # or registerDoParallel() if you want to use all the cores
	
	pr<-cbind(pca=rep(PCAS,times=length(RES)),res=rep(RES,each=length(PCAS)))
	pcas<-foreach(opts=iter(pr,by='row'))  %dopar%  
		Seurat::FindClusters(seuratobj, pc.use = 1:opts[1], resolution = opts[2], print.output = 0, save.SNN = F)@ident
	
	names(pcas)=paste(rep(PCAS,times=length(RES)),res=rep(RES,each=length(PCAS)),sep="_")
	#tsnes<-foreach(opts=PCAS)  %dopar%  Seurat::RunTSNE(seuratobj, dims.use = 1:opts, do.fast = T,perplexity=perplexity,max_iter=max_iter)@tsne.rot
	tsnes<-foreach(opts=PCAS)  %dopar%  Seurat::RunTSNE(seuratobj, dims.use = 1:opts, do.fast = T,perplexity=perplexity,dim_embed = dim)@tsne.rot
	tsneReps<-rep(tsnes,times=length(RES))
	names(tsneReps)=paste(rep(PCAS,times=length(RES)),res=rep(RES,each=length(PCAS)),sep="_")
	
	return(list(PCAS=PCAS,RES=RES,pcas=pcas,tsneReps=tsneReps,pr=pr))
}
##


########################################################################################################################################
# functions for comparing Differential Gene Expression by cluster

########################################################################################################################################
# functions for comparing Differential Gene Expression by cluster
InternalDGE<-function(SeuratObj,compare.by,subset=NULL,minPct=0.1,min.diff.pct=-Inf,ths=0.25,cellmin=10,posOnly=F){ # compare.by is a logical by which to compare gene expression
	
	SeuratObj@data.info$diff<-compare.by # logical vector for comparisons
	if(!is.null(subset)) {SeuratObj=SubsetData(SeuratObj,cells.use = subset)}
	
	print(table(SeuratObj@data.info$cluster,SeuratObj@data.info$treat))
	
	dges<-NULL #dge object
	
	ftest<-data.frame(matrix(ncol=8,nrow=length(unique(SeuratObj@data.info$cluster)))) #make dataframe with 8 columns & one row for each cluster
	colnames(ftest)<-c("Cluster","OFF","ON","pval","OR","lowCI","highCI","sig?")
	rownames(ftest)<-sort(unique(SeuratObj@data.info$cluster))
	
	for (clusNum in sort(unique(SeuratObj@data.info$cluster))){ #for each cluster
		
		InClus.cells=rownames(SeuratObj@data.info[SeuratObj@data.info$cluster==clusNum,]) # cell IDs for cells in that cluster
		print(paste("Cluster",clusNum,":",length(InClus.cells),"Cells"))
		
		# Test if cluster composition is skewed with respect to compare.by vector
		ct<-table(SeuratObj@data.info$diff,SeuratObj@data.info$cluster==clusNum); colnames(ct)<-c("InCluster","OutCluster"); rownames(ct)<-c("TRUE","FALSE"); # make contingency table
		ft<-fisher.test(ct,conf.level = (1-0.05/length(unique(SeuratObj@data.info$cluster)))) # fisher test p.val adjusted by cluster number
		ftest[clusNum,]<-c(clusNum,ct[1:2,1],scientific(ft$p.value,2),round(ft$estimate,2),round(ft$conf.int[1],2),round(ft$conf.int[2],2),sig=(ft$p.value<0.05/length(unique(SeuratObj@data.info$cluster))))
		# #ct<-cbind(rbind(ct,Total=colSums(ct)),Total=rowSums(rbind(ct,Total=colSums(ct))))
		
		newIdent=rep("IRR",dim(SeuratObj@data)[2]); names(newIdent)=rownames(SeuratObj@data.info); # make new ident vector of IRR
		DIFF=SeuratObj@data.info[match(InClus.cells,rownames(SeuratObj@data.info)),]$diff # get the comparison vector for this cluster
		
		if(min(table(DIFF))>=cellmin & length(table(DIFF))>1){  # makes sure both compare.by values are in the cluster and present in at least cellmin # of cells
			print(clusNum)
			newIdent[match(InClus.cells,colnames(SeuratObj@data))]=DIFF # for cells in the cluster assign DIFF as new identity
			SeuratObj@ident=newIdent
			
			cluster.markers <- FindMarkers(SeuratObj, ident.1 = "TRUE",ident.2 = "FALSE",only.pos = posOnly, min.pct = minPct,thresh.use = ths,min.diff.pct=min.diff.pct) #find markers
			ordered=cluster.markers[order(cluster.markers$avg_diff,decreasing = T),] #order by average difference
			ord<-cbind(geneID=rownames(ordered),cluster=rep(clusNum,dim(cluster.markers)[1]),ordered)
			dges<-rbind(dges,ord) #add to dge object
		}
	}
	rownames(dges)<-NULL
	dges$pctDiff<-dges$pct.1-dges$pct.2
	dges$pctChng<-round(abs(dges$pctDiff)/apply(dges[,c("pct.1","pct.2")],1,max),3)
	return(list(ftest,dges))
}
#########################################################################################################################################
timepoint.cor<-function(CorMe,CorWith,genelist=NULL,method="pearson"){
	if(is.null(genelist)){
		genelist<-intersect(rownames(CorMe),rownames(CorWith))
	}else{
		genelist<-intersect(rownames(CorMe[genelist]),rownames(CorWith[genelist,]))
	}
	cor.value<-cor(CorMe[genelist,],CorWith[genelist,])
	return(cor.value)
}
########################################################################################################################################
MyMerge <- function(x, y){
	df <- merge(x, y, by="row.names", all=TRUE)
	rownames(df)=df$Row.names; df=df[,-1];df[is.na(df)]=0
	return(df)
}
########################################################################################################################################
coX<-function(genePair,exprMat,thresh=0,corMethod = "spearman"){
	require("psych",quietly=F)
	exprMat<-as.matrix(exprMat)
	geneA<-as.character(genePair[1])
	geneB<-as.character(genePair[2])
	
	# if genes not above threshold or expressed in all cells, function returns a NULL row that won't show up in a data frame being generated when applied
	nullreturn<-data.frame(geneA=character(),geneB=character(),phi=numeric(),OddsRatio=numeric(),p.value=numeric(),cor=numeric(),BonInAon=numeric(),BonInAoff=numeric(),meanDGE=numeric())
	if(mean(exprMat[geneA,]>0)<=thresh | mean(exprMat[geneA,]>0)==1){return(nullreturn)}
	if(mean(exprMat[geneB,]>0)<=thresh | mean(exprMat[geneB,]>0)==1){return(nullreturn)}
	
	#print(paste(geneA,geneB,"by=",unique(exprMat["by",])))
				
	conTable<-table(exprMat[geneA,]!=0,exprMat[geneB,]!=0) #create contingency table
	#print(conTable)
	ft<-fisher.test(conTable) # calculate fishers exact test
	coreff<-round(cor(exprMat[geneA,],exprMat[geneB,],method=corMethod),3) #calculate correlation coefficient
	BonInAon<-round(conTable[2,2]/sum(conTable[2,]),3) #percentage in ON
	BonInAoff<-round(conTable[1,2]/sum(conTable[1,]),3) #percentage in OFF
	meanDGE<-round(log(mean(exprMat[geneB,exprMat[geneA,]>0])/mean(exprMat[geneB,exprMat[geneA,]==0]),2),3) #mean difference in gene expression
	
	# return GeneA GeneB Phi OddsRatio p.value and correlation coefficient
	return(data.frame(geneA=geneA,geneB=geneB,phi=round(phi(conTable),3),OddsRatio=round(as.numeric(ft$estimate),3),p.value=ft$p.value,cor=coreff,BonInAon=BonInAon,BonInAoff=BonInAoff,meanDGE=meanDGE))
}	
########################################################################################################################################
multicoX<-function(geneListA,geneListB=NULL,exprMat,thresh=0,by=NULL,p.adjust="fdr",parallel=FALSE, cores=29,corMethod ="spearman"){
	require(plyr)
  if(is.null(geneListB)){genePairs<-data.frame(t(combn(geneListA,2)))} # if no geneB do all combinations of GeneA
	else{ #if geneB is set 
		geneListB<-geneListB[!(geneListB %in% geneListA)] # Eliminate duplicate genes
		if(geneListB[1]=="all"){genePairs<-expand.grid(geneListA,rownames(exprMat))} #if geneB is "all", compare all the genes
		else{genePairs<-expand.grid(geneListA,geneListB)}
	}
	genePairs %>% dplyr::mutate_if(is.factor, as.character) -> genePairs # return factors columns to character
	genePairs<-genePairs[genePairs[,1]!=genePairs[,2],] # remove self-correlations
	genePairs=as.list(data.frame(t(genePairs)))
	
	if(parallel){
		registerDoParallel(cores=20) # or registerDoParallel() if you want to use all the cores
	}
	if(is.null(by)){
		cx<-plyr::ldply(.data = genePairs,.fun = coX,exprMat=exprMat,thresh=thresh,corMethod,.id=NULL,.parallel = parallel)
		cx$p.value<-p.adjust(cx$p.value,method=p.adjust)		
	}else{
		by<-factor(by)
		funx<-function(fd,genePairs,thresh,corMethod = "spearman"){  # then for each separate transposed expression matrix
			df<-t(fd)  # retranspose
			# if we can remove these first they dont need to run
			cxx<-plyr::ldply(.data = genePairs,.fun = coX,exprMat=df,thresh=thresh,corMethod,.id=NULL,.parallel = parallel)
			cxx$p.value<-p.adjust(cxx$p.value,method=p.adjust) #adjust p.value
			return(cxx) # return to ddply
		}
		taMrpxe<-cbind(by,t(as.matrix(exprMat)))
		rownames(taMrpxe)<-colnames(exprMat)
		cx<-ddply(data.frame(taMrpxe),.(by),.fun=funx, genePairs,thresh=thresh,corMethod)
	}
	return(cx[order(cx$phi),])
}
########################################################################################################################################
getInfo<-function(genes, cols=c("geneID","name","description"),file=NULL,sep="\t"){  # need failure mode for genes not in the list
	if(is.null(file)){
		missing<-genes[!(genes %in% GeneInfo$geneID)]
	
		if(cols[1]=="all"){return(GeneInfo[sapply(genes,grep,GeneInfo$geneID),])}
		else{return(GeneInfo[sapply(genes,grep,GeneInfo$geneID),cols])}
	}
	else{
		tab<-read.csv(file,sep="\t",stringsAsFactors =FALSE)
		if(cols[1]=="all"){return(tab[sapply(genes,grep,GeneInfo$geneID),])}
		else{return(tab[sapply(genes,grep,GeneInfo$geneID),cols])}
	}
}
########################################################################################################################################
RunSeurat<-function(exprMat,cells=NULL,id="all"){  #cells is vector of cell names
	if(is.null(cells)){clusterInClus.cells<- colnames(exprMat)}else {clusterInClus.cells<- cells}
	data_merge<-all_merge[,cells]
	SO <- new("seurat", raw.data = data_merge)
	SO <- Seurat::Setup(	SO,min.cells = 3,min.genes = 1,do.logNormalize = T,total.expr = 1e4, project = as.character(id))
	SO <- Seurat::SubsetData(SO, subset.name = "nUMI", accept.high = 5000, accept.low = 300)
	SO <- Seurat::RegressOut(SO, latent.vars = c("nUMI"))
	SO <- Seurat::MeanVarPlot(	SO ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.03, x.high.cutoff = 10, y.cutoff = 0.9,do.contour = F,do.plot = F)
	SO <- Seurat::PCA(SO, pc.genes = SO@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
	SO <- Seurat::ProjectPCA(SO)
	SO <- Seurat::FindClusters(SO, pc.use = 1:24, resolution = 1.0, print.output = 0, save.SNN = F)
	SO <- Seurat::RunTSNE(SO, dims.use = 1:24, do.fast = T)
	return(SO)
}
########################################################################################################################################
evenSample<-function(plotme,vect,minCount=NULL){
	if(is.null(minCount)) minCount<-min(table(vect))
	plotme$vect<-vect
	newplotme<-data.frame()
	for(i in unique(plotme$vect))newplotme<-rbind(newplotme,dplyr::sample_n(plotme[plotme$vect==i,],size=minCount))
	newplotme$vect<-NULL
	return(newplotme)
}
########################################################################################################################################
gg_color_hue <- function(n) {
	hues = seq(15, 375, length = n + 1)
	hcl(h = hues, l = 65, c = 100)[1:n]
}
########################################################################################################################################
coxPlotme<-function(SO,geneA,geneB,cells=NULL){
	if(!is.null(cells))SO<-SubsetData(SO,cells.use=cells)
	pltme<-cbind(
		SO@data.info,
		SO@tsne.rot,
		exprA=(exp(SO@data[geneA,])-1),
		onA=SO@data[geneA,]>0,
		exprB=(exp(SO@data[geneB,])-1),
		onB=SO@data[geneB,]>0
	)
	pltme$cox<-factor(as.numeric(pltme$onA)*1+as.numeric(pltme$onB)*2)
	dict<-data.frame(a=1:4,b=factor(c("Neither",as.character(paste(geneA,"alone")),as.character(paste(geneB,"alone")),"Co-Expr."),levels=c("Neither",as.character(paste(geneA,"alone")),as.character(paste(geneB,"alone")),"Co-Expr.")))
	pltme$cox<-dict[as.numeric(pltme$cox),"b"]
	return(pltme)
}
########################################################################################################################################
# returns top expressed genes for a cluster
topExpressed<-function(SO,cluster,sortByPctPos=T,top=25){
	if(sortByPctPos){
		exprMat<-SO@data[,rownames(SO@data.info[SO@data.info$cluster==cluster,])]
		pctPos<-round(sort(rowMeans(as.matrix(exprMat)>0),decreasing = T)[1:top],3)*100
		meanExpr<-round(exp(rowMeans(as.matrix(exprMat[names(pctPos),])))-1,1)
		returnme<-data.frame(row.names=NULL,cluster=rep(cluster,length(pctPos)),geneID=names(pctPos),pctPos=pctPos,meanExpr=meanExpr)
	}
	else{
		exprMat<-SO@data[,rownames(SO@data.info[SO@data.info$cluster==cluster,])]
		meanExpr<-round(sort(exp(rowMeans(as.matrix(exprMat)),decreasing = T)-1)[1:top],3)
		pctPos<-round(rowMeans(as.matrix(exprMat[names(meanExpr),])>0),1)*100
		returnme<-data.frame(row.names=NULL,cluster=rep(cluster,length(meanExpr)),geneID=names(meanExpr),meanExpr=meanExpr,pctPos=pctPos)
	}
	return(returnme)
}

#####################################################################################################################################
#defines how many breaks to have in ggplot axis
equal_breaks <- function(n = 4, ...){
  function(x){
    step=round(floor(max(x))/(n-1),1)
    floor(c(0,1:(n-1)*step))
  }
}