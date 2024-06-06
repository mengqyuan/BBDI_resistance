#########################figure 5a###########################
#########################figure 5a###########################
library("Seurat")
library(cluster)
load("SUM159.Rdata")
used=subset(SUM159,idents="SUM159RDMSO")
count=used@assays$RNA@counts
rawcount=CreateSeuratObject(counts=count[which(! duplicated(rownames(count))),],min.cells=3,min.features=200)
rawcount=NormalizeData(rawcount, normalization.method = "LogNormalize", scale.factor = 10000)
rawcount=FindVariableFeatures(rawcount, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rawcount)
scaledata <- ScaleData(rawcount,features = all.genes,verbose=FALSE)
scaledata <- RunPCA(scaledata, features=VariableFeatures(object=scaledata))
scaledata <- FindNeighbors(scaledata, reduction = "pca", dims = 1:10)
emb=Embeddings(scaledata@reductions$pca)
dist.matrix=dist(emb[,1:10])
score=c()

for(i in seq(0.1,1.5,0.1)){
scaledata<-scaledata%>% FindClusters(resolution = i)
#DimPlot(scaledata)
clusters=scaledata@meta.data$seurat_clusters
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)

score=rbind(score,c(mean(sil[,3]),sd(sil[,3])))
}
colnames(score)=c("mean_sihouette","sd")
score=as.data.frame(score)
score$res=seq(0.1,1.5,0.1)
score1=score[1:10,]

scaledata <-scaledata%>% FindClusters(resolution = 0.5)
scaledata <- RunUMAP(scaledata, reduction = "pca", dims = 1:10)
DimPlot(scaledata,cols=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F"),pt.size=1.5,label=F)
save(scaledata,file="SUM159_subtype.Rdata")
pdf("fig5a.pdf",width=5,height=4)
DimPlot(scaledata,cols=c("#C4A35C","#ED7EA4","#E37444","#61BBD4","#91C79F"),pt.size=1,label=F)
dev.off()

#########################figure 5b###########################
#########################figure 5b###########################

library("Seurat")
gene=read.table("volcaplot.gene",sep="\t")
genes=gene[,1]
up_index=c(1:8)
down_index=c(9:20)
load("SUM159_subtype.Rdata")
used=scaledata@assays$RNA@data[genes,]
score=colSums(used[up_index,])-colSums(used[down_index,])
re=as.data.frame(cbind(score,as.vector(Idents(scaledata))))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
SUM159_subtype=boxplot(score~Group,data=re,col=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F"),outline=F)

pdf("fig5b.pdf")
boxplot(score~Group,data=re,col=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F"),outline=F)
dev.off()
#########################figure 5c###########################
#########################figure 5c###########################
load("SUM159.Rdata")
write.table(SUM159@assays$RNA@data,file="SUM159_exp.tsv",sep="\t",quote=F)
##############R
library("Seurat")
data=read.csv("log2_IC50.pred.csv")
jq1=data[,c(1,which(colnames(data)=="X163"))]
rownames(jq1)=jq1[,1]
rownames(jq1)=gsub(".","-",rownames(jq1),fixed=T)
load("SUM159_subtype.Rdata")
RDMSO=jq1[is.element(rownames(jq1),colnames(scaledata)),]
identical(rownames(RDMSO),colnames(scaledata))#TRUE
RDMSO_subtype=as.data.frame(cbind(RDMSO[,2],as.vector(Idents(scaledata))))

colnames(RDMSO_subtype)=c("log2_IC50","Group")
RDMSO_subtype[,1]=as.numeric(RDMSO_subtype[,1])
pdf("fig5c.pdf")
boxplot(log2_IC50~Group,data=RDMSO_subtype,outline=F,col=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F"))
dev.off()
#########################figure 5d###########################
#########################figure 5d###########################
library(monocle)
library(Seurat)
#library(SeuratWrappers)
library(ggplot2)
require(dplyr)
library("openxlsx")
load("SUM159_subtype.Rdata")
subset_data=scaledata
new.metadata <- merge(subset_data@meta.data,data.frame(Idents(subset_data)),by = "row.names",sort = FALSE)
rownames(new.metadata)<-new.metadata[,1]
head(subset_data@meta.data)
subset_data@meta.data<-new.metadata
expression_matrix <- as(as.matrix(subset_data@assays$RNA@counts), 'sparseMatrix')
cell.type<-as.vector(Idents(subset_data))
cell_metadata <- new('AnnotatedDataFrame',data=cbind(subset_data@meta.data,cell.type))# 132510
gene_annotation <- new('AnnotatedDataFrame',data=data.frame(gene_short_name = row.names(subset_data), row.names = row.names(subset_data)))
monocle_cds <- newCellDataSet(expression_matrix,phenoData = cell_metadata,featureData = gene_annotation,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())		      

cds <- monocle_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)  ## Removing 99 outliers
diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~ cell.type")
### inference the pseudotrajectory########################################################
# step1: select genes for orderding setOrderingFilter() #
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))#
cds <- setOrderingFilter(cds, ordering_genes)  
plot_ordering_genes(cds) #mark the genes ###
# step2: dimension reduction=> reduceDimension()  DDRTree #
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
# step3: ordering the cells=> orderCells()
cds <- orderCells(cds)

pdf("fig5d-1.pdf")
plot_cell_trajectory(cds, color_by = "cell.type") + scale_color_manual(breaks = c("0","1", "2", "3","4"), values=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F"))
dev.off()


cds <- orderCells(cds,root_state=2)
pdf("fig5d-2.pdf")
plot_cell_trajectory(cds, color_by = "Pseudotime")+scale_color_gradientn(colours=c("#440154","#3C508B","#228E8D","#53C569","#EEE51C")) 
dev.off()

#########################figure 5e###########################
#########################figure 5e###########################
library(Seurat)
load("SUM159_subtype.Rdata")
exp=as.matrix(scaledata@assays$RNA@counts)
pheno=as.vector(Idents(scaledata))
library(CytoTRACE)
re=CytoTRACE(mat=exp)
plotCytoTRACE(re,phenotype=pheno,outputDir="SUM159")
#########################figure 5f###########################
#########################figure 5f###########################
library("Seurat")
load("SUM159_subtype.Rdata")
gene=read.table("volcaplot.gene",sep="\t")
genes=gene[,1]

p=DotPlot(scaledata,features=genes,cols = c("lightgrey", "#F07C53"),dot.scale=3)+scale_x_discrete("")+scale_y_discrete("")+theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.4))+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

pdf("fig5f.pdf",width=4,height=4)
p
dev.off()


#########################figure 5g###########################
#########################figure 5g###########################
load("SUM159_subtype.Rdata")
color=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F")
VlnPlot(scaledata,features=c("FTH1"),pt.size=0,stack=T,fill.by="ident",cols=rev(color))+guides(fill=F)+theme(axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 0),
        axis.text.x = element_text(size = 0,color ='black',angle=0),
        axis.text.y = element_text(size = 10,color ='black',angle=0),axis.ticks=element_blank(),panel.border = element_rect(colour = "black"),axis.line=element_blank())
pdf("fig5g-1.pdf",width=3,height=3)		
RidgePlot(scaledata,features=c("FTH1"),cols=color)
dev.off()
pdf("fig5g-2.pdf",width=3,height=2)
RidgePlot(scaledata,features=c("FTL"),cols=color)
dev.off()
pdf("fig5g-3.pdf",width=3,height=3)
RidgePlot(scaledata,features=c("GPX4"),cols=color)
dev.off()
#########################figure 5h###########################
#########################figure 5h###########################
library(GSVA)
library(GSEABase)
library(msigdbr)
library(Seurat)
load("SUM159_subtype.Rdata")
exp=as.matrix(scaledata@assays$RNA@data)
msgdC2 = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
keggSet = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_description)
res_ssGSEA <- gsva(as.matrix(exp),keggSet,method = "ssgsea",kcdf = "Gaussian",min.sz= 5,max.sz=1000,ssgsea.norm= T,parallel.sz= 1L)

Idents(scaledata)
out0=c()
for(i in 1:dim(res_ssGSEA)[1]){
temp=wilcox.test(res_ssGSEA[i,which(Idents(scaledata)==0)],res_ssGSEA[i,which(Idents(scaledata)!=0)],alternative="greater")
out0=rbind(out0,c(rownames(res_ssGSEA)[i],temp$p.value))
}


out1=c()
for(i in 1:dim(res_ssGSEA)[1]){
temp=wilcox.test(res_ssGSEA[i,which(Idents(scaledata)==1)],res_ssGSEA[i,which(Idents(scaledata)!=1)],alternative="greater")
out1=rbind(out1,c(rownames(res_ssGSEA)[i],temp$p.value))
}

out2=c()
for(i in 1:dim(res_ssGSEA)[1]){
temp=wilcox.test(res_ssGSEA[i,which(Idents(scaledata)==2)],res_ssGSEA[i,which(Idents(scaledata)!=2)],alternative="greater")
out2=rbind(out2,c(rownames(res_ssGSEA)[i],temp$p.value))
}



out3=c()
for(i in 1:dim(res_ssGSEA)[1]){
temp=wilcox.test(res_ssGSEA[i,which(Idents(scaledata)==3)],res_ssGSEA[i,which(Idents(scaledata)!=3)],alternative="greater")
out3=rbind(out3,c(rownames(res_ssGSEA)[i],temp$p.value))
}


out4=c()
for(i in 1:dim(res_ssGSEA)[1]){
temp=wilcox.test(res_ssGSEA[i,which(Idents(scaledata)==4)],res_ssGSEA[i,which(Idents(scaledata)!=4)],alternative="greater")
out4=rbind(out4,c(rownames(res_ssGSEA)[i],temp$p.value))
}


all=cbind(out0[,2],out1[,2],out2[,2],out3[,2],out4[,2])
all=apply(all,2,as.numeric)
rownames(all)=rownames(res_ssGSEA)
plotpath=unique(c(out0[order(out0[,2]),][1:10,1],out1[order(out1[,2]),][1:10,1],out2[order(out2[,2]),][1:10,1],out3[order(out3[,2]),][1:10,1],out4[order(out4[,2]),][1:10,1]))

filter=c()
for(i in 1:dim(all)[1]){
	if(length(which(all[i,]<0.05))==1){
		filter=c(filter,i)
	}
}
filter_path=all[filter,]

t1=names(filter_path[order(filter_path[,1]),1][1:5])
t2=names(filter_path[order(filter_path[,2]),2][1:4])
t3=names(filter_path[order(filter_path[,3]),3][1:5])
t4=names(filter_path[order(filter_path[,4]),4][1:5])
t5=names(filter_path[order(filter_path[,5]),5][1:5])
final=filter_path[c(t1,t2,t3,t4,t5),]
final[which(final>0.05)]=1
final[which(final<0.05)]=0

colnames(final)=c("R1","R2","R3","R4","R5")
library("pheatmap")
p1=pheatmap(final,scale="none",show_colnames=T,show_rownames=T,cluster_cols=F,cluster_rows=F,color = c("#EA5415","#E1E3EC"))
pdf("ssgsea.heatmap.teyi.pdf")
p1
dev.off()
