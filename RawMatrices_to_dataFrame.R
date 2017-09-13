#This code compiles all the single-cell gene expression matrices into one matrix, and creates a metadata object used for cell filtering in subsequent scripts
setwd("") #Set to the folder where all dge matrices are found
lapply(list("reshape2","plyr","dplyr"),require,character.only = TRUE)
source("codesForWeb/InternalFunctions.R")

########################################################################################################
#This script is used to load the matrices which are the DropSeq Toolbox output. The cells are placed in one matrix and a metadata structure is created
#label structure: AP2G.A.42.OFF.A (Strain,Batch,time,treatment,tRep)

# peg4 tdT NF54 Schizonts (1700 cells)
dge1<-read.table("./out_gene_exon_tagged.dge_NF54re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)

#correct exogenous genes
data_merge=dge1
# Change cellIDs and sampleType
orgID<-colnames(dge1)
newID<-paste("NF54.A.42.Ind.A",1:dim(dge1)[2],sep="_")
cellID<-cbind(orgID,newID,ldply(strsplit(ldply(strsplit(newID,"_"))[,1],".",fixed=T)))
colnames(data_merge)<-newID
colnames(cellID)<-c("orgID","newID","strain","exp","tp","treat","tRep")
cellID<-cbind(cellID,AP2G=rep("ON",dim(cellID)[1]),GDV1=rep("ON",dim(cellID)[1]));
cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
NF54_schiz_merge<-data_merge
NF54_schiz_ID<-cellID
NF54_schiz_ID$barcode=NF54_schiz_ID$orgID
rm(dge1,data_merge,orgID,newID,cellID)
#########################################################################################################
# AP2-G-ddFKBP Schizonts

data_CN2_1=read.table("./out_gene_exon_tagged.dge_CN2InducedON1re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_CN2_2=read.table("./out_gene_exon_tagged.dge_CN2InducedON2re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_CN3_1=read.table("./out_gene_exon_tagged.dge_CN3InducedOFF1re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_CN3_2=read.table("./out_gene_exon_tagged.dge_CN3InducedOFF2re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_CN4=read.table("./out_gene_exon_tagged.dge_CN4GCD2re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)

data_merge=Reduce(MyMerge, list(data_CN2_1,data_CN2_2,data_CN3_1,data_CN3_2,data_CN4))
barcodes=colnames(data_merge)
colnames(data_merge)=c(paste0("CN21_",1:ncol(data_CN2_1)),paste0("CN22_",1:ncol(data_CN2_2)),paste0("CN31_",1:ncol(data_CN3_1)),paste0("CN32_",1:ncol(data_CN3_2)),paste0("CN4_",1:ncol(data_CN4)))

# Change cellIDs and sampleType
pairing<-cbind(c("CN21","CN22","CN31","CN32","CN4"),c("AP2G.A.42.ON.A","AP2G.A.42.ON.B","AP2G.A.42.OFF.A","AP2G.A.42.OFF.B","AP2G.A.GC2.ON.A"))

orgID<-colnames(data_merge)
newID<-paste(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],ldply(strsplit(orgID,"_"))[,2],sep="_")
cellID<-cbind(orgID,newID,ldply(strsplit(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],".",fixed=T)))
colnames(data_merge)<-newID
colnames(cellID)<-c("orgID","newID","strain","exp","tp","treat","tRep")
cellID<-cbind(cellID,AP2G=cellID$treat,GDV1=rep("ON",dim(cellID)[1])); 
cellID$AP2G[cellID$AP2G=="UION"]<-"ON"
cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
exp52_merge<-data_merge
exp52_ID<-cellID
exp52_ID$barcode=barcodes
rm(data_merge,pairing,orgID,newID,cellID,barcodes)
########################################################################################################
# DCJ BSD
pairing<-cbind(c("BSD1","BSD2"),c("DCJ.E.18.BSD.A","DCJ.E.18.BSD.B"))

dge1<-read.csv("./out_gene_exon_tagged.dge_BSD1re.txt.gz",sep="\t",stringsAsFactors = F); rownames(dge1)<-dge1[,1]; dge1<-dge1[,-1]; #colnames(dge1)<-paste(rep(pairing[1,2],dim(dge1)[2]),1:dim(dge1)[2],sep="_")
dge2<-read.csv("./out_gene_exon_tagged.dge_BSD2re.txt.gz",sep="\t",stringsAsFactors = F); rownames(dge2)<-dge2[,1]; dge2<-dge2[,-1]; #colnames(dge2)<-paste(rep(pairing[2,2],dim(dge1)[2]),1:dim(dge1)[2],sep="_")
data_merge<-Reduce(MyMerge, list(dge1,dge2))
barcodes=colnames(data_merge)
colnames(data_merge)=c(paste0("BSD1_",1:ncol(dge1)),paste0("BSD2_",1:ncol(dge2)))
# Change cellIDs and sampleType
orgID<-colnames(data_merge)
newID<-c(paste(rep(pairing[1,2],dim(dge1)[2]),1:dim(dge1)[2],sep="_"),paste(rep(pairing[2,2],dim(dge2)[2]),1:dim(dge2)[2],sep="_"))
cellID<-cbind(orgID,newID,ldply(strsplit(ldply(strsplit(newID,"_"))[,1],".",fixed=T)))
colnames(data_merge)<-newID
colnames(cellID)<-c("orgID","newID","strain","exp","tp","treat","tRep")
cellID<-cbind(cellID,AP2G=rep("ON",dim(cellID)[1]),GDV1=rep("ON",dim(cellID)[1])); 
cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
exp54_Rings_merge<-data_merge
exp54_Rings_ID<-cellID
exp54_Rings_ID$barcode=barcodes
rm(dge1,dge2,data_merge,pairing,orgID,newID,cellID,barcodes)
#########################################################################################################
# AP2-G-ddFKBP timecourse
data_30Plus=read.table("./out_gene_exon_tagged.dge_30Pre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_30Minus=read.table("./out_gene_exon_tagged.dge_30Mre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_36PlusA=read.table("./out_gene_exon_tagged.dge_36PAre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_36PlusB=read.table("./out_gene_exon_tagged.dge_36PBre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_36MinusA=read.table("./out_gene_exon_tagged.dge_36MAre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_36MinusB=read.table("./out_gene_exon_tagged.dge_36MBre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42PlusA=read.table("./out_gene_exon_tagged.dge_42PAre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42PlusB=read.table("./out_gene_exon_tagged.dge_42PBre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42MinusA=read.table("./out_gene_exon_tagged.dge_42MAre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42MinusB=read.table("./out_gene_exon_tagged.dge_42MBre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42D2A=read.table("./out_gene_exon_tagged.dge_42D2Are.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42D2B=read.table("./out_gene_exon_tagged.dge_42D2Bre.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)

data_merge=Reduce(MyMerge, list(data_30Plus,data_30Minus,data_36MinusA,data_36MinusB,data_36PlusA,data_36PlusB,data_42MinusA,data_42MinusB,data_42PlusA,data_42PlusB,data_42D2A,data_42D2B))
barcodes=colnames(data_merge)
colnames(data_merge)=c(paste0("30P_",1:ncol(data_30Plus)),paste0("30M_",1:ncol(data_30Minus)),paste0("36MA_",1:ncol(data_36MinusA)),paste0("36MB_",1:ncol(data_36MinusB)),paste0("36PA_",1:ncol(data_36PlusA)),paste0("36PB_",1:ncol(data_36PlusB)),paste0("42MA_",1:ncol(data_42MinusA)),paste0("42MB_",1:ncol(data_42MinusB)),paste0("42PA_",1:ncol(data_42PlusA)),paste0("42PB_",1:ncol(data_42PlusB)),paste0("42D2A_",1:ncol(data_42D2A)),paste0("42D2B_",1:ncol(data_42D2B)))

#correct exogenous genes
pairing<-cbind(
  c("30P","30M","36MA","36MB","36PA","36PB","42MA","42MB","42PA","42PB","42D2A","42D2B"),
  c("AP2G.B.30.ON.A","AP2G.B.30.OFF.A","AP2G.B.36.OFF.A","AP2G.B.36.OFF.B","AP2G.B.36.ON.A","AP2G.B.36.ON.B","AP2G.B.42.OFF.A","AP2G.B.42.OFF.B","AP2G.B.42.ON.A","AP2G.B.42.ON.B","AP2G.B.GC2.ON.A","AP2G.B.GC2.ON.B")
)
orgID<-colnames(data_merge)
newID<-paste(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],ldply(strsplit(orgID,"_"))[,2],sep="_")
cellID<-cbind(orgID,newID,ldply(strsplit(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],".",fixed=T)))
colnames(data_merge)<-newID
colnames(cellID)<-c("orgID","newID","strain","exp","tp","treat","tRep")
cellID<-cbind(cellID,AP2G=cellID$treat,GDV1=rep("ON",dim(cellID)[1])); 
cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
exp54_merge<-data_merge
exp54_ID<-cellID
exp54_ID$barcode=barcodes
rm(data_merge,pairing,orgID,newID,cellID,barcodes)

#########################################################################################################
# GC stage IV (day 9) timecourse

data_NFGC9=read.table("./out_gene_exon_tagged.dge_GCD9re.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)

data_merge=data_NFGC9
barcodes=colnames(data_merge)
colnames(data_merge)=paste0("GCD9_",1:ncol(data_NFGC9))

pairing<-cbind(c("GCD9"),c("NF54.C.GC9.Ind.A"))
orgID<-colnames(data_merge)
newID<-paste(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),1],ldply(strsplit(orgID,"_"))[,2],sep="_")
cellID<-cbind(orgID,newID,ldply(strsplit(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],".",fixed=T)))
colnames(data_merge)<-newID
colnames(cellID)<-c("orgID","newID","strain","exp","tp","treat","tRep")
cellID<-cbind(cellID,AP2G=rep("ON",dim(cellID)[1]),GDV1=cellID$treat);
cellID$GDV1[1:500]<-"ON"
cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
gdv1_cellID<-cellID
gdv1_merge<-data_merge
gdv1_cellID$barcode=barcodes
rm(data_merge,pairing,orgID,newID,cellID,barcodes)

#########################################################################################################
# NF54 Time Course
data_42A=read.table("out_gene_exon_tagged.dge_NF42A.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_42B=read.table("out_gene_exon_tagged.dge_NF42B.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_48A=read.table("out_gene_exon_tagged.dge_NF48A.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)
data_48B=read.table("out_gene_exon_tagged.dge_NF48B.txt.gz", header=T,row.names = "GENE",stringsAsFactors=F)

data_merge=Reduce(MyMerge, list(data_42A,data_42B,data_48A,data_48B))
barcodes=colnames(data_merge)
colnames(data_merge)=c(paste0("NF42A_",1:ncol(data_42A)),paste0("NF42B_",1:ncol(data_42B)),paste0("NF48A_",1:ncol(data_48A)),paste0("NF48B_",1:ncol(data_48B)))

pairing<-cbind(
  c("NF42A","NF42B","NF48A","NF48B"), 
  c("NF54.tc.42.Ind.A","NF54.tc.42.Ind.B","NF54.tc.48.Ind.A","NF54.tc.48.Ind.B")
)
orgID<-colnames(data_merge)
newID<-paste(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],ldply(strsplit(orgID,"_"))[,2],sep="_")
cellID<-cbind(orgID,newID,ldply(strsplit(pairing[match(ldply(strsplit(colnames(data_merge),"_"))[,1],pairing[,1]),2],".",fixed=T)))
colnames(data_merge)<-newID
colnames(cellID)<-c("orgID","newID","strain","exp","tp","treat","tRep")
cellID<-cbind(cellID,AP2G=rep("ON",dim(cellID)[1]),GDV1=rep("ON",dim(cellID)[1]));
cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
NFtc_cellID<-cellID
NFtc_merge<-data_merge
NFtc_cellID$barcode=barcodes
rm(data_merge,pairing,orgID,newID,cellID,barcodes)

#########################################################################################################
# Combine the experiments
cellID<-rbind(NF54_schiz_ID,exp52_ID,exp54_Rings_ID,exp54_ID,gdv1_cellID,NFtc_cellID)
cellID$orgID<-as.character(cellID$orgID);cellID$newID<-as.character(cellID$newID);cellID$AP2G<-as.character(cellID$AP2G);cellID$GDV1<-as.character(cellID$GDV1);

#cellID$sample<-apply(cellID[,3:7],1,paste,collapse=".")
apply(cellID[,-1:-2],2,table)
all_merge=Reduce(MyMerge, list(NF54_schiz_merge,exp52_merge,exp54_Rings_merge,exp54_merge,gdv1_merge,NFtc_merge))
identical(colnames(all_merge),cellID$newID)

save(all_merge,cellID,file="All_merged_forWeb.RData")


