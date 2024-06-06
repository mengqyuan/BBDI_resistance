#########################figure 4a###########################
#########################figure 4a###########################

library(dplyr)
library(cowplot)
library(Seurat)
#InferCNV is used to explore tumor single cell RNA-Seq data to identify evidence for large-scale chromosomal copy number variations, #such as gains or deletions of entire chromosomes or large segments of chromosomes. This is done by exploring expression intensity of #genes across positions of the genome in comparison to the average or a set of reference ‘normal’ cells. A heatmap is generated #illustrating the relative expression intensities across each chromosome, and it becomes readily apparent as to which regions of the #genome are over-abundant or less-abundant as compared to normal cells (or the average, if reference normal cells are not provided).
library("infercnv")
infercnvsample<-function(matrix,anno,gene_order,out_dir,ref_group){
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=matrix,
                                    annotations_file=anno,
                                    delim="\t",
                                    gene_order_file=gene_order,
                                    ref_group_names=ref_group,
				    chr_exclude = c("chrM"))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=TRUE, 
                             plot_steps=F,
                             mask_nonDE_genes = F,
		             denoise=TRUE,
		             HMM=TRUE,
					 analysis_mode="subclusters"
#                           include.spike=T  # used for final scaling to fit range (0,2) centered at 1.
                             )
return(infercnv_obj)
}

ref=as.matrix(Read10X(data.dir = "hg19"))
length(intersect(which(colSums(ref)>=15500),which(colSums(ref)<=16000))) #105
index=intersect(which(colSums(ref)>=15500),which(colSums(ref)<=16000))
ref_100=ref[,index]
ref_100=as.data.frame(ref_100)
ref_100$symbol=rownames(ref_100)
load("SUM159.Rdata")
library("Seurat")
scaledata=SUM159
used=subset(scaledata,idents=c("SUM159DMSO","SUM159RDMSO"))
SUM159count=used@assays$RNA@counts
datause=as.data.frame(SUM159count)
datause$symbol=rownames(SUM159count)
data_merge=merge(ref_100,datause,all=TRUE,by="symbol")
rownames(data_merge)=data_merge$symbol
data_merge <- subset( data_merge, select = -c(symbol))
data_merge[is.na(data_merge)] <- 0
ref_100 <- subset( ref_100, select = -c(symbol))
ref_anno=cbind(colnames(ref_100),rep("ref",105))
write.table(ref_anno,"anno159.txt",quote=F,col.names=F,row.names=F,sep="\t")
ob_anno=cbind(colnames(used),as.character(Idents(used)))
write.table(ob_anno,"anno159.txt",quote=F,col.names=F,row.names=F,sep="\t",append=T)
ref_group="ref"
outdir="SUM159DMSO_RDMSO_chrX_chrY"
anno="SUM159DMSO_RDMSO_chrX_chrY/anno159.txt"
gene_order="gencode_v19_gene_pos.txt"

infercnvobj=infercnvsample(data_merge,anno,gene_order,outdir,ref_group)


library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library(scales)

mat=read.table("infercnv.observations.txt",header=T,row.names=1,sep=" ",stringsAsFactors=F)
gn <- rownames(mat)
geneFile <- read.table("gencode_v19_gene_pos.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(geneFile$V1,gn),]
mat=mat[intersect(geneFile$V1,gn),]

top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = c(1:22,"X"),labels_gp = gpar(cex = 1.5)))

gene=read.table("volcaplot.gene",header = F,sep = "\t",stringsAsFactors = F)
key_gene=intersect(rownames(mat),gene[,1])
ha=columnAnnotation(foo=anno_mark(at=which(rownames(mat) %in% key_gene),labels = rownames(mat)[which(rownames(mat) %in% key_gene)],which = "column",side = "bottom"))
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = c(1:22,"X"),labels_gp = gpar(cex = 1.5)))

forCNVanno=cbind(c(rep("DMSO",length(grep(".5",colnames(mat)))),rep("RDMSO",length(grep(".7",colnames(mat))))),colnames(mat))
colnames(forCNVanno)=c("subgroup","cell")
rownames(forCNVanno)=forCNVanno[,2]
forCNVanno=as.data.frame(forCNVanno)
forCNVanno$subgroup=factor(forCNVanno$subgroup,levels = sort(unique(forCNVanno$subgroup)))
color_c=RColorBrewer::brewer.pal(12, "Paired")[5:6]
color_c=c("#8491B4","#FDC086")
names(color_c)=names(table(forCNVanno$subgroup))
forCNVanno$cell=NULL
left_anno <- rowAnnotation(df = forCNVanno,col=list(subgroup=color_c),border = F)
pdf("fig4a-1.pdf",width = 15,height = 11)
ht_tp = Heatmap(t(mat),
                col = colorRamp2(c(0.6,1,1.4), c("#54AEB0","#F0F0F0","#FCDA14")),
                cluster_rows = F,cluster_columns = F,
                show_column_names = F,show_row_names = F,
                column_split = factor(sub_geneFile$V2, paste("chr",c(1:22,"X"),sep = "")),
                column_gap = unit(2, "mm"),
                top_annotation = top_anno,
                heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
                bottom_annotation = ha,
                left_annotation = left_anno,
                row_title = NULL,
)

draw(ht_tp, heatmap_legend_side = "right")
dev.off()


mat=read.table("infercnv.references.txt",header=T,row.names=1,sep=" ",stringsAsFactors=F)

gn <- rownames(mat)
geneFile <- read.table("gencode_v19_gene_pos.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(geneFile$V1,gn),]
mat=mat[intersect(geneFile$V1,gn),]

top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = c(1:22,"X"),labels_gp = gpar(cex = 1.5)))
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = c(1:22,"X"),labels_gp = gpar(cex = 1.5)))


pdf("fig4a-2.pdf",width = 15,height = 11)
ht_tp = Heatmap(t(mat),
                col = colorRamp2(c(0.6,1,1.4), c("#54AEB0","#F0F0F0","#FCDA14")),
                cluster_rows = F,cluster_columns = F,
                show_column_names = F,show_row_names = F,
                column_split = factor(sub_geneFile$V2, paste("chr",c(1:22,"X"),sep = "")),
                column_gap = unit(2, "mm"),
                top_annotation = top_anno,
                heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
                row_title = NULL,
)

draw(ht_tp, heatmap_legend_side = "right")
dev.off()



#########################figure 4b###########################
#########################figure 4b###########################
library(Seurat)
load("SUM159.Rdata")
data=subset(SUM159,idents=c("SUM159DMSO","SUM159RDMSO"))
data_use=data@assays$RNA@counts
write.csv(t(as.matrix(data_use)),file="SUM159_res_scenic_input.csv")
#########python########
import os, sys
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("SUM159_res_scenic_input.csv");#R中导出的表达矩阵
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("SUM159_res.loom",x.X.transpose(),row_attrs,col_attrs)
#########shell########
pyscenic grn SUM159_res.loom hs_hgnc_tfs.txt -o SUM159_res_scenic_input.adj.tsv
pyscenic ctx SUM159_res_scenic_input.adj.tsv hg19-500bp-upstream-7species.mc9nr.feather -o SUM159_res_scenic_input.motif.tsv --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname SUM159_res.loom
pyscenic aucell SUM159_res.loom SUM159_res_scenic_input.motif.tsv  --output SUM159_res_SCENIC.loom --num_workers 3


##########R########################################
library(SCENIC)
packageVersion("SCENIC")
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
#library(scRNAseq)
rm(list=ls())

loom <- open_loom('SUM159_res_SCENIC.loom')      
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons_incidMat[1:4,1:4]  
regulons <- regulonsToGeneLists(regulons_incidMat)    
geneall=read.table("volcaplot.gene",header = F,sep = "\t",stringsAsFactors = F)
up=geneall[1:8,]
down=geneall[9:20,]
net=c()
fer_tf=c()
fer_tf_new=c()
for(i in 1:length(regulons)){
	tf=names(regulons)[i]
	tf=gsub("[(+)]","",tf)
	temp=tf
	for(j in 1:length(regulons[i][[1]])){
		gene=regulons[i][[1]][j]
		if(is.element(gene,geneall[,1])){
			temp=paste(temp,gene,sep="_")
			fer_tf=c(fer_tf,tf)
			net=rbind(net,c(tf,gene))
		}
		
		
	}
	if(temp!=tf){
	fer_tf_new=c(fer_tf_new,temp)}
}


regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
tfname=names(regulonAUC)
tfname_use=gsub("[(+)]","",tfname)
heatmap_gene=intersect(tfname_use,unique(fer_tf)) 
index=is.element(tfname_use,heatmap_gene)
aucvalue=(regulonAUC@assays@data@listData$AUC)
plot_value=aucvalue[index,]
allpvaluegreater=c()
allpvalueless=c()
retain=c()
up_tf=c()
down_tf=c()
for(i in 1:dim(plot_value)[1]){
	p_greater=wilcox.test(plot_value[i,1:751],plot_value[i,752:1881],alternative ="greater")
	p_less=wilcox.test(plot_value[i,1:751],plot_value[i,752:1881],alternative ="less")
	allpvaluegreater=c(allpvaluegreater,p_greater$p.value)
	allpvalueless=c(allpvalueless,p_less$p.value)
	if(p_less$p.value<0.01){
		retain=c(retain,i)
		up_tf=c(up_tf,rownames(plot_value)[i])
	}
	if(p_greater$p.value<0.01){
		retain=c(retain,i)
		down_tf=c(down_tf,rownames(plot_value)[i])
	}
}
pvalue_out=cbind(rownames(plot_value),allpvaluegreater,allpvalueless) #ATF3(+) 0.102872547324011
up_tf_use=gsub("[(+)]","",up_tf)
down_tf_use=gsub("[(+)]","",down_tf)
print(length(retain)) 
plot_value_re=plot_value[retain,]
net_filter_temp=gsub("[(+)]","",rownames(plot_value_re))
net_filter=net[is.element(net[,1],net_filter_temp),]
net_filter[,2]=gsub("ATF4","ATF4_1",net_filter[,2])
net_filter=cbind(net_filter,1)
net_filter=net_filter[,c(2,1,3)]
colnames(net_filter)=c("term","genes","logFC")
library("pheatmap")
load("SUM159.Rdata")
data=subset(SUM159,idents=c("SUM159DMSO","SUM159RDMSO"))
identical(colnames(data),colnames(regulonAUC))
Idents(data)
table(Idents(data))
annotation_col = data.frame(DataType = Idents(data))
rownames(annotation_col)<-colnames(plot_value_re)
ann_colors = list(DataType= c(SUM159DMSO="#8491B4",SUM159RDMSO="#FDC086"))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
result = plot_value_re[p$tree_row$order,]
name=gsub("[(+)]","",rownames(result))

result = pvalue_out[p3$tree_row$order,][c(1:11,14:15,42:43,20,24,30,12:13,16:19,21:23,25:29,31:41,44:46),1]
p4=pheatmap(plot_value[result,],scale="row",show_colnames=F,cluster_cols=F,cluster_rows=F,color = c(colorRampPalette(colors = c("#234B44","#234B44","#234B44","#F0F0F0","white"))(length(bk)/2),colorRampPalette(colors =c("white","#F0F0F0","#FBAC5F","#FBAC5F","#FBAC5F"))(length(bk)/2)),legend_breaks=seq(-2,2,1),breaks=bk,annotation_col=annotation_col,annotation_colors=ann_colors)
pdf("Fig4b.pdf")
p4
dev.off()
#########################figure 4c###########################
#########################figure 4c###########################

library(ggraph)
library(ggcor)
library(tidygraph)
library(netET)
library("igraph")
library(vegan)
data("varespec")
data("varechem")
geneall=read.table("volcaplot.gene",header = F,sep = "\t",stringsAsFactors = F)
geneall[10,1]="ATF4_1"
fer_up=geneall[1:8,1]
fer_up=fer_up[-c(2,8)]
fer_down=geneall[9:20,1]

up_tf=read.table("pyscenic_up_tf.txt")
down_tf=read.table("pyscenic_down_tf.txt")

df <- read.csv("pyscenic_tf_gene_network_fer.txt")
head(df)
df[,3]=0
for(i in 1:dim(df)[1]){
	if(is.element(df[i,1],fer_up) & is.element(df[i,2],up_tf[,1])){df[i,3]=1}
	if(is.element(df[i,1],fer_down) & is.element(df[i,2],down_tf[,1])){df[i,3]=1}
}
#nodes <- data.frame(name = unique(union(df[,1], df[,2])))
nodes=data.frame(name = c(fer_up,up_tf[,1],fer_down,down_tf[,1]),up=c(rep(0,length(c(fer_up,up_tf[,1]))),rep(1,length(c(fer_down,down_tf[,1])))))

inode=c(names(sort(table(df[,1])[fer_up])),names(sort(table(df[,1])[fer_down],decreasing=T)))
onode=c(names(sort(table(df[,2])[up_tf[,1]])),names(sort(table(df[,2])[down_tf[,1]],decreasing=T)))
#onode=c(names(sort(table(df[,2])[up_tf[,1]])),names(sort(table(df[,2])[down_tf[,1]])))
p <- tbl_graph(nodes=nodes,edges=df[,1:3]) %>%
  mutate(Degree = centrality_degree(mode = "all")) %>%
  as_bipartite_circular(start=0,end=180,inner_nodes=inode,outer_nodes =onode )
D <- ggraph(p, layout_bipartite_circular(p)) +
  annotate_arc_rect(0, 180, 
                    fill = "#e0eee8", 
                    r0 = 0.55, 
                    r1 = 1.05) +
  geom_edge_circular(aes(colour = logFC > 0), edge_width = 0.5, edge_alpha = 0.8) +
  geom_node_point(aes(size = Degree, colour = up == 0)) +
  geom_node_text_circular(expand = 0.08,size = 3) +
  scale_size(range = c(1,3)) +
  scale_colour_manual(values = c("TRUE" = "#FF6A48","FALSE" = "#234B44"),
                      guide = "none") +
  scale_edge_colour_manual(values = c("TRUE" = "#F5A396", "FALSE" = "#8491B4"),
                           labels = c("In-Consistent", "consistent")) +
  coord_fixed(clip = "off", xlim = c(-1.2, 1.2), ylim = c(0, 1.1)) +
  theme(panel.background = element_blank()) +
  guides(edge_colour = guide_legend(override.aes = list(edge_width = 1))) +
  labs(edge_colour = "Spearman's r")
pdf("fig4c.pdf", width = 6, height = 6)
D
dev.off()
#########################figure 4d###########################
#########################figure 4d###########################

library("Seurat")
gene=read.table("volcaplot.gene",sep="\t")
genes=gene[,1]
up_index=c(1:8)
down_index=c(9:20)
setwd("/mnt/yuanmengqin/Single_cell_resistance/Step2/up_down_risk_score")
load("/mnt/yuanmengqin/Single_cell_resistance/Step1/SUM159.Rdata")
ab=subset(SUM159,idents=c("SUM159DMSO","SUM159RDMSO"))
used=ab@assays$RNA@data[genes,]
rawcount=ab@assays$RNA@data
rawcount=CreateSeuratObject(counts=rawcount,min.cells=0,min.features=0)
rawcount=NormalizeData(rawcount, normalization.method = "LogNormalize", scale.factor = 10000)
rawcount=FindVariableFeatures(rawcount, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rawcount)
scaledata <- ScaleData(rawcount,features = all.genes,verbose=FALSE)
scaledata <- RunPCA(scaledata, features=VariableFeatures(object=scaledata))
scaledata <- RunUMAP(scaledata, reduction = "pca", dims = 1:10)
Idents(scaledata)=as.vector(Idents(ab))
score=colSums(used[up_index,])-colSums(used[down_index,])
re=as.data.frame(cbind(score,as.vector(Idents(ab))))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
scaledata$fer20=re[,1]
library(RColorBrewer)
pdf("fig4d.pdf")
FeaturePlot(scaledata,features="fer20",cols=rev(brewer.pal(11,name="RdBu")))
dev.off()


#########################figure 4e###########################
#########################figure 4e###########################
load("SUM149.Rdata")
ab=subset(SUM149,idents=c("SUM149DMSO","SUM149RDMSO"))
ge=intersect(genes,rownames(SUM149))
up_index=c(1:8)
down_index=c(9:19)
used=ab@assays$RNA@data[ge,]
rawcount=ab@assays$RNA@data
rawcount=CreateSeuratObject(counts=rawcount,min.cells=0,min.features=0)
rawcount=NormalizeData(rawcount, normalization.method = "LogNormalize", scale.factor = 10000)
rawcount=FindVariableFeatures(rawcount, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rawcount)
scaledata <- ScaleData(rawcount,features = all.genes,verbose=FALSE)
scaledata <- RunPCA(scaledata, features=VariableFeatures(object=scaledata))
scaledata <- RunUMAP(scaledata, reduction = "pca", dims = 1:10)
Idents(scaledata)=as.vector(Idents(ab))
score=colSums(used[up_index,])-colSums(used[down_index,])
re=as.data.frame(cbind(score,as.vector(Idents(ab))))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
scaledata$fer20=re[,1]
save(scaledata,file="SUM149_fer20.Rdata")
library(RColorBrewer)
col=c("#5163AD","#5163AD","#3581B8","#54AEB0","#B2DFA0","#E3F392","#FCDA14","#E65747","#B92546","#B92546")
pdf("fig4e.pdf",width=5,height=5)
FeaturePlot(scaledata,features="fer20",cols=col)
dev.off()

#########################figure 4f###########################
#########################figure 4f###########################
library(tidyverse)
library(data.table)
tmp=read.table("gene_length.txt",sep="\t",header=T)
tmp=tmp[,-1]
tmp=tmp[!duplicated(tmp[,1]),]
count0=read.table("GSM3763465_SUM159_DMSO_12h_R1.counts.txt",sep="\t")
count1=read.table("GSM3763466_SUM159_DMSO_12h_R2.counts.txt",sep="\t")
count2=read.table("GSM3763467_SUM159_DMSO_24h_R1.counts.txt",sep="\t")
count3=read.table("GSM3763468_SUM159_DMSO_24h_R2.counts.txt",sep="\t")
count4=read.table("GSM3763477_SUM159_JQ1R_DMSO_24h_R1.counts.txt",sep="\t")
count5=read.table("GSM3763478_SUM159_JQ1R_DMSO_24h_R2.counts.txt",sep="\t")
count6=read.table("GSM3763479_SUM159_JQ1R_DMSO_3h_R1.counts.txt",sep="\t")
count7=read.table("GSM3763480_SUM159_JQ1R_DMSO_3h_R2.counts.txt",sep="\t")

count=cbind(count0,count1[,2],count2[,2],count3[,2],count4[,2],count5[,2],count6[,2],count7[,2])
colnames(count)=c("gene_symbol","DMSO1","DMSO2","DMSO3","DMSO4","RDMSO1","RDMSO2","RDMSO3","RDMSO4")
tmp <- left_join(count,tmp,by = 'gene_symbol')
tmp=tmp[!is.na(tmp[,3]),]
fpkm <- data.frame(row.names = tmp$gene_symbol)
for (i in 2:(dim(tmp)[2]-1)){
  col <- as.numeric(tmp[[i]])
  N <- sum(col) # 计算每个样本的mapped reads数
  FPKMi <- (col*1e9)/(N*as.numeric(tmp[[dim(tmp)[2]]]))# 计算FPKM值
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}

gene=read.table("volcaplot.gene")
genes=gene[,1]
used=fpkm[genes,]
#used=log2(used+1)
score=colSums(used[1:8,])-colSums(used[9:20,])
re=as.data.frame(cbind(score,c("Control","Control","Control","Control","Resistant","Resistant","Resistant","Resistant")))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
pdf("fig4f-1.pdf",height=5,width=3)
boxplot(score~Group,data=re,col=c("#FDC086","#99DAA9"),outline=F)
dev.off()

library(tidyverse)
library(data.table)
tmp=read.table("gene_length.txt",sep="\t",header=T)
tmp=tmp[,-1]
tmp=tmp[!duplicated(tmp[,1]),]
count0=read.table("GSM3763447_SUM149_DMSO_12h_R1.counts.txt",sep="\t")
count1=read.table("GSM3763448_SUM149_DMSO_12h_R2.counts.txt",sep="\t")
count2=read.table("GSM3763449_SUM149_DMSO_12h_R3.counts.txt",sep="\t")
count3=read.table("GSM3763456_SUM149R_DMSO_12h_R1.counts.txt",sep="\t")
count4=read.table("GSM3763457_SUM149R_DMSO_12h_R2.counts.txt",sep="\t")
count5=read.table("GSM3763458_SUM149R_DMSO_12h_R3.counts.txt",sep="\t")
count=cbind(count0,count1[,2],count2[,2],count3[,2],count4[,2],count5[,2])
colnames(count)=c("gene_symbol","DMSO1","DMSO2","DMSO3","RDMSO1","RDMSO2","RDMSO3")
tmp <- left_join(count,tmp,by = 'gene_symbol')
tmp=tmp[!is.na(tmp[,3]),]
fpkm <- data.frame(row.names = tmp$gene_symbol)
for (i in 2:(dim(tmp)[2]-1)){
  col <- as.numeric(tmp[[i]])
  N <- sum(col) 
  FPKMi <- (col*1e9)/(N*as.numeric(tmp[[dim(tmp)[2]]]))
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() 
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}

gene=read.table("volcaplot.gene")
genes=gene[,1]
used=fpkm[genes,]
score=colSums(used[1:8,])-colSums(used[9:20,])
#re=as.data.frame(cbind(score,c("Resistant","Resistant","Resistant","Control","Control","Control")))
re=as.data.frame(cbind(score,c("Control","Control","Control","Resistant","Resistant","Resistant")))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
pdf("fig4f-2.pdf",height=5,width=3)
boxplot(score~Group,data=re,col=c("#FDC086","#99DAA9"),outline=F)
dev.off()

#########################figure 4g###########################
#########################figure 4g###########################
load("GSE164813_exp.Rdata")
gene=read.table("volcaplot.gene")
genes=gene[,1]
used=GSE_expr[genes,]
used[is.na(used)]=0
score=colSums(used[1:8,])-colSums(used[9:20,])
re=as.data.frame(cbind(score,c("Resistant","Resistant","Resistant","Control","Control","Control","Control Treated","Control Treated","Control Treated","Resistant","Resistant","Resistant","Control","Control","Control","Control Treated","Control Treated","Control Treated")))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
use=re[c(1:6,10:15),]
use=cbind(use,c(rep("H23",6),rep("H1975",6)))
colnames(use)[3]="cell line"
library(ggpubr)
use[4:6,2]="Sensitive"
use[10:12,2]="Sensitive"
p=ggbarplot(use, x = "cell line", y = "score", 
          add = c("mean_se"),
          color = "Group",fill = "Group", palette = c("#FDC086","#8491B4"),position = position_dodge(0.8))+
stat_compare_means(aes(group = Group), label = "p.format",method = "t.test") +
  #rotate_x_text(angle = 45) +
  ylab("Score") + #
  theme(legend.position="right") +  
  theme(legend.key.size = unit(0.3, "cm"))
  
pdf("fig4g.pdf",width=5,height=5)
p+NoLegend()
dev.off()

#########################figure 4h###########################
#########################figure 4h###########################
RV1=read.table("GSM2753412_22RV1.genes.fpkm_tracking",fill=TRUE,sep="\t",header=T)
RV1=RV1[,c(1,10)]
RV1_1=read.table("GSM2753413_22RV1JQ1RPool1.genes.fpkm_tracking",fill=TRUE,sep="\t",header=T)
RV1_2=read.table("GSM2753414_22RV1JQ1RPool2.genes.fpkm_tracking",fill=TRUE,sep="\t",header=T)
RV1_1=RV1_1[,c(1,10)]
RV1_2=RV1_2[,c(1,10)]
LNCaP=read.table("GSM2753415_LNCaP.genes.fpkm_tracking",fill=TRUE,sep="\t",header=T)
LNCaP_1=read.table("GSM2753416_LNCaPJQ1RPool1.genes.fpkm_tracking",fill=TRUE,sep="\t",header=T)
LNCaP_2=read.table("GSM2753417_LNCaPJQ1RPool2.genes.fpkm_tracking",fill=TRUE,sep="\t",header=T)
LNCaP=LNCaP[,c(1,10)]
LNCaP_1=LNCaP_1[,c(1,10)]
LNCaP_2=LNCaP_2[,c(1,10)]
identical(LNCaP_1[,1],LNCaP[,1])
identical(LNCaP_1[,1],RV1[,1])
exp=cbind(RV1,RV1_1[,2],RV1_2[,2],LNCaP[,2],LNCaP_1[,2],LNCaP_2[,2])
rownames(exp)=exp[,1]
exp=exp[,-1]
gene=read.table("volcaplot.gene")
genes=gene[,1]
used=exp[genes,]
used[is.na(used)]=0
score=colSums(used[1:8,])-colSums(used[9:20,])
re=as.data.frame(cbind(score,c("Control","Resistant","Resistant","Control","Resistant","Resistant")))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
use=cbind(re,c(rep("RV1",3),rep("LNCaP",3)))
colnames(use)[3]="cell line"
use[1,2]="Sensitive"
use[4,2]="Sensitive"
p=ggbarplot(use, x = "cell line", y = "score", 
          #add = c("mean_se"),
          color = "Group",fill = "Group", palette = c("#FDC086","#8491B4"),position = position_dodge(0.8))+
stat_compare_means(aes(group = Group), label = "p.format",method = "t.test") +
  #rotate_x_text(angle = 45) +
  ylab("Score") + 
  theme(legend.position="right") +  
  theme(legend.key.size = unit(0.3, "cm"))
  
pdf("fig4h.pdf",width=5,height=5)
p+NoLegend()
dev.off()

#########################figure 4i###########################
#########################figure 4i###########################

drug_response=read.csv("GCT-1.3_AUC.gct",sep="\t")
jq_cell=drug_response[which(rownames(drug_response)=="JQ_1"),]
jq_cell=jq_cell[-which(is.na(jq_cell)==TRUE)]#length(jq_cell)

data_1=read.table("CCLE_rpkm.csv",sep="\t")
data_used=data_1[is.element(data_1[,2],genes),]
dim(data_used)#20
length(unique(rownames(data_used[,1:2])))#20
rownames(data_used)=data_used[,2]
exp=data_used[,c(-1,-2)] 
colnames(exp)=unlist(lapply(strsplit(colnames(exp), split = "_"), function(x) x[1]))#
exp=exp[genes,]
common_sample=intersect(names(jq_cell),colnames(exp))
exp_used=exp[,common_sample]
jq_cell_used=jq_cell[common_sample]
score_1=colSums(exp_used[up_index,])-colSums(exp_used[down_index,])
yuzhi=quantile(score_1)
up=names(which(score_1>=yuzhi[4]))
down=names(which(score_1<yuzhi[2]))

jq_cell_used=jq_cell_used[c(up,down)]
class=c(rep("High",length(up)),rep("Low",length(down)))
re=as.data.frame(cbind(as.numeric(jq_cell_used),class))
colnames(re)=c("AUC","Group")
re[,1]=as.numeric(re[,1])

CCLE_score=boxplot(AUC~Group,data=re,col=c("#FDC086","#99DAA9"),outline=F)
pdf("fig4i-1.pdf",height=5,width=3)
boxplot(AUC~Group,data=re,col=c("#FDC086","#99DAA9"),outline=F)
dev.off()
library(ggpubr)
compare_means(AUC~Group,data=re) # 0.055

##########GDSC
data=readRDS("GDSC_cell_exp.RDS")
drug_response=read.csv("GCT-1.3_AUC.gct",sep="\t")
jq_cell=drug_response[which(rownames(drug_response)=="JQ_1"),]
jq_cell=jq_cell[colnames(data)]
jq_cell=jq_cell[-which(is.na(jq_cell)==TRUE)]
comm=intersect(names(jq_cell),colnames(data))
jq_cell_used=jq_cell[comm]
data_used=data[,comm]
data_used=data_used[genes,]
data_used[is.na(data_used)]=0
score_1=colSums(data_used[up_index,])-colSums(data_used[down_index,])
yuzhi=quantile(score_1)
up=names(which(score_1>=yuzhi[4]))
down=names(which(score_1<yuzhi[2]))

jq_cell_used=jq_cell_used[c(up,down)]
class=c(rep("High",length(up)),rep("Low",length(down)))
re=as.data.frame(cbind(as.numeric(jq_cell_used),class))
colnames(re)=c("AUC","Group")
re[,1]=as.numeric(re[,1])

pdf("fig4i-2.pdf",height=5,width=3)
boxplot(AUC~Group,data=re,col=c("#FDC086","#99DAA9"),outline=F)
dev.off()


