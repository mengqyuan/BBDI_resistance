#########################figure 2a###########################
#########################figure 2a###########################
library(dplyr)
library(Seurat)
library(patchwork)
load("GSE131135_BRD4JQ1_sc10x_raw_counts.RData")
rawCount=rawCounts
rownames(rawCounts)=rawCounts.ann[,2] #ENSG change to gene symbol
single_data_process<-function(rawcount,idents,umap_or_tsne)
{
	rawcount=CreateSeuratObject(counts=rawcount[which(! duplicated(rownames(rawcount))),],min.cells=3,min.features=200)
	dim(rawcount@assays$RNA@counts)
	Idents(rawcount) <- idents
	rawcount[["percent.mt"]]<-PercentageFeatureSet(rawcount, pattern = "^MT-")
	rawcount=subset(rawcount, subset = nFeature_RNA > 200 & percent.mt < 20)
	rawcount=NormalizeData(rawcount, normalization.method = "LogNormalize", scale.factor = 10000)
	rawcount=FindVariableFeatures(rawcount, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(rawcount)
	scaledata <- ScaleData(rawcount,features = all.genes,verbose=FALSE)
	scaledata <- RunPCA(scaledata, features=VariableFeatures(object=scaledata))
	scaledata <- RunUMAP(scaledata, reduction = "pca", dims = 1:10)
	if(umap_or_tsne=="tsne")
	{
	        scaledata <- RunTSNE(scaledata, reduction = "pca", dims = 1:10)
	}
	return(scaledata)
}


idents=as.vector(rawCounts.grp[c(3657:8111)])
SUM159=single_data_process(rawCounts[,3657:8111],idents,"umap")
save(SUM159,file="SUM159.Rdata")
pdf("fig2a.pdf")
DimPlot(SUM159, cols=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4"))
dev.off()

#########################figure 2b###########################
#########################figure 2b###########################
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
 
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}
SUM159DMSO_marker=FindMarkers(SUM159, ident.1="SUM159DMSO")
SUM159JQ1_marker=FindMarkers(SUM159, ident.1="SUM159JQ1")
SUM159RDMSO_marker=FindMarkers(SUM159, ident.1="SUM159RDMSO")
SUM159RJQ1_marker=FindMarkers(SUM159, ident.1="SUM159RJQ1")
library("dplyr")
library("Seurat")
top20_1_1=rownames(SUM159DMSO_marker %>% top_n(n = 20, wt = avg_log2FC))
top20_2_1=rownames(SUM159JQ1_marker %>% top_n(n = 20, wt = avg_log2FC))
top20_3_1=rownames(SUM159RDMSO_marker %>% top_n(n = 20, wt = avg_log2FC))
top20_4_1=rownames(SUM159RJQ1_marker %>% top_n(n = 20, wt = avg_log2FC))
select_gene_1=unique(c(top20_1_1,top20_2_1,top20_3_1,top20_4_1))
DoHeatmap(SUM159,features=select_gene_1,angle=30,size=3)

data=SUM159@assays$RNA@scale.data[select_gene_1,]
library("pheatmap")
celltype=gsub("SUM159","",Idents(SUM159))
annotation_col= data.frame(celltype=factor(celltype))
rownames(annotation_col) = colnames(SUM159)
color=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4")
ann_colors = list(celltype = c(DMSO=color[1],JQ1=color[2],RDMSO=color[3],RJQ1=color[4]))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
colors=colorRampPalette(c('#66CC66','#B3CAC8','#ffffff','#F1C341',"#EA8405"))(length(bk))
p=pheatmap(data,cluster_rows=F,cluster_cols = F,color = colors,annotation_col = annotation_col,annotation_colors=ann_colors,show_colnames = F,legend_breaks=seq(-2,2,1),breaks=bk,fontsize=6)

gene_name=c("CCND1", "FLNA", "TFDP1", "TUBA1B", "ZFP36L2", "ID3","ASNS","FTH1", "FTL", "GPX4", "HMOX1", "SAT1", "SLC3A2", "GSTP1","SMS")
pdf("Fig2b.pdf")
add.flag(p,kept.labels = gene_name,repel.degree = 0)
dev.off()
add.flag(p,kept.labels = gene_name,repel.degree = 0.2)


#########################figure 2c###########################
#########################figure 2c###########################
mfuzz_plot_defined=function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
    time.points, ylim.set = c(0, 0), xlab = "Time", ylab = "Expression changes", 
    x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", 
    col.lab = "black", col.main = "black", col.sub = "black", 
    col = "black", centre = FALSE, centre.col = "black", centre.lwd = 2, 
    Xwidth = 5, Xheight = 5, single = FALSE, line.lwd=1,...) 
{
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < min.mem] <- -1
    colorindex <- integer(dim(exprs(eset))[[1]])
    if (missing(colo)) {
        colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800", 
            "#FF6000", "#FF7800", "#FF8F00", "#FFA700", "#FFBF00", 
            "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", 
            "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", 
            "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28", 
            "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F", 
            "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", 
            "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", 
            "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF", 
            "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF", 
            "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", 
            "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", 
            "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018")
    }
    else {
        if (colo == "fancy") {
		fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                c(255:150))
            colo <- rgb(b = fancy.blue/255, g = fancy.green/255, 
                r = fancy.red/255)
        }
    }
    colorseq <- seq(0, 1, length = length(colo))
    for (j in 1:dim(cl[[1]])[[1]]) {
        if (single) 
            j <- single
        tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
            if (x11) 
                X11(width = Xwidth, height = Xheight)
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            }
            else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            if (sum(ylim.set == c(0, 0)) == 2) {
                ylim <- c(ymin, ymax)
            }
            else {
                ylim <- ylim.set
				}
            if (!is.na(sum(mfrow))) {
                par(mfrow = mfrow, bg = bg, col.axis = col.axis, 
                  col.lab = col.lab, col.main = col.main, col.sub = col.sub, 
                  col = col)
            }
            else {
                par(bg = bg, col.axis = col.axis, col.lab = col.lab, 
                  col.main = col.main, col.sub = col.sub, col = col)
            }
            xlim.tmp <- c(1, dim(exprs(eset))[[2]])
            if (!(missing(time.points))) 
                xlim.tmp <- c(min(time.points), max(time.points))
            plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                xlab = xlab, ylab = ylab, main = paste("Cluster", 
                  j), axes = FALSE, ...)
            if (missing(time.labels) && missing(time.points)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
                  col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.labels) && !(missing(time.points))) {
                axis(1, time.points, 1:length(time.points), time.points, 
                  col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.points) & !(missing(time.labels))) {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels,
				col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (!(missing(time.points)) & !(missing(time.labels))) {
                axis(1, time.points, time.labels, col = ax.col, 
                  ...)
                axis(2, col = ax.col, ...)
            }
        }
        else {
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            }
            else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            if (sum(ylim.set == c(0, 0)) == 2) {
                ylim <- c(ymin, ymax)
            }
            else {
                ylim <- ylim.set
            }
            xlim.tmp <- c(1, dim(exprs(eset))[[2]])
            if (!(missing(time.points)))
xlim.tmp <- c(min(time.points), max(time.points))
            plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                xlab = xlab, ylab = ylab, main = paste("Cluster", 
                  j), axes = FALSE, ...)
            if (missing(time.labels) && missing(time.points)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
                  col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.labels) && !(missing(time.points))) {
                axis(1, time.points, 1:length(time.points), time.points, 
                  col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.points) & !(missing(time.labels))) {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
                  col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (!(missing(time.points)) & !(missing(time.labels))) {
                axis(1, time.points, time.labels, col = ax.col, 
                  ...)
                axis(2, col = ax.col, ...)
            }
        }
        if (length(tmpmem) > 0) {
            for (jj in 1:(length(colorseq) - 1)) {
tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                  colorseq[jj + 1])
                if (sum(tmpcol) > 0) {
                  tmpind <- which(tmpcol)
                  for (k in 1:length(tmpind)) {
                    if (missing(time.points)) {
                      lines(tmp[tmpind[k], ], col = colo[jj],lwd = line.lwd)
                    }
                    else lines(time.points, tmp[tmpind[k], ], 
                      col = colo[jj],lwd = line.lwd)
                  }
                }
            }
        }
        if (centre) {
            lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
        }
        if (single) 
            return()
    }
}



library(dplyr)
library(Seurat)
library(patchwork)
load("SUM159.Rdata")
marker_exp=AverageExpression(object =SUM159)
diff=FindMarkers(SUM159, ident.1="SUM159DMSO")
diff1=FindMarkers(SUM159, ident.1="SUM159JQ1")
diff2=FindMarkers(SUM159, ident.1="SUM159RDMSO")
diff3=FindMarkers(SUM159, ident.1="SUM159RJQ1")
#################|logFC|>0.25 p<0.01
all=rbind(diff,diff1,diff2,diff3)
temp=all[which(as.numeric(all[,2])>=0.25),] 
up=temp[which(as.numeric(temp[,5])<= 0.01),]
temp=all[which(as.numeric(all[,2])<= -0.25),] 
down=temp[which(as.numeric(temp[,5])<= 0.01),]
gene=unique(c(rownames(up),rownames(down)))
data_4=as.matrix(marker_exp$RNA[gene,])
eset <- new("ExpressionSet",exprs = data_4)
tmp <- filter.std(eset,min.std=0) 
data.s <- standardise(tmp) 
m1 <- mestimate(data.s)
set.seed(2022101702)
cl_8 <- mfuzz(data.s,c=6,m=m1)##fuzzy c-means 
gene_filter=c()
cluster_filter=c()
for(i in 1:dim(cl_8$membership)[1])
{
	if(max(cl_8$membership[i,])>0.6){
		gene_filter=c(gene_filter,rownames(cl_8$membership)[i])
		cluster_filter=c(cluster_filter,cl_8$cluster[i])
	}
}
peusodtime_gene=cbind(gene_filter,cluster_filter)#816
write.table(peusodtime_gene,"mfuzz_peusodtime_gene_159.txt",quote=F,sep="\t",row.names=F)
gene_cluster <- cbind(cl_8$cluster, data_4)
colnames(gene_cluster)[1] <- 'cluster'
plot_line=exprs(data.s)
plot_line=cbind(plot_line,cl_8[[3]])
colnames(plot_line)[5]="cluster"


color=c("#E3DCA3","#CADABD","#9AA193","#A1B2CC","#CB9C7A","#8A7B93","#AF777E","#F18772")
for(i in 1:6)
{
	outfile=paste("SUM159_6_",i,".pdf",sep="")
	color1=rep(color[i],64)
	pdf(outfile)
	mfuzz_plot_defined(data.s,cl=cl_8,mfrow=c(1,1),single=i,centre=T,colo=color1,line.lwd=2,time.labels=c("SUM159DMSO","SUM159JQ1","SUM159RDMSO","SUM159RJQ1"),x11 = FALSE,centre.lwd=3,min.mem =0.6)
	dev.off()
}
#########################figure 2d###########################
#########################figure 2d###########################
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggplot2")
outFile = function(temp,outFilePath,outFiletxt,outPDFPath){
	temp$function.gene.num = unlist(strsplit(temp[,"BgRatio"],"/"))[seq(1,length(unlist(strsplit(temp[,"BgRatio"],"/"))),2)]
	temp$GeneRatio = temp$Count / as.numeric(temp$function.gene.num)
	write.table(temp,outFiletxt,col.names=T,row.names=F,quote=F,sep="\t")
	if(dim(temp)[1]>=10){temp1 = temp[1:10,]}else{temp1 = temp[1:dim(temp)[1],]}
	temp1 = temp1[order(temp1$GeneRatio, decreasing=T),]
	temp1$Description = factor(temp1$Description,levels=rev(temp1$Description))
	gap = (max(temp1$GeneRatio) - min(temp1$GeneRatio))/5
	pdf(outPDFPath)
	p = ggplot(temp1,mapping = aes(x=GeneRatio,y=Description))+geom_point(aes(size=Count,color=p.adjust))+theme_bw()+ scale_color_gradient(low = "blue", high = "red")+scale_x_continuous(expand = c(gap, gap))+ylab(NULL)
	print(p)
	dev.off()
}
data=read.table("mfuzz_peusodtime_gene_159.txt",header=T)
rownames(data)=data[,1]
data1=read.table("mfuzz_peusodtime_gene_159_metascape_entrezid.txt",fill=T)
rownames(data1)=data1[,1]
used=merge(data,data1,all.x=TRUE,by="row.names")
for (i in 1:6){
filedir=paste(outdir,"/","cluster_",i,sep="")
dir.create(filedir,recursive=T)
setwd(filedir)
DEG.entrez=used[which(used[,3]==i),5]
ego.KEGG = enrichKEGG(gene=DEG.entrez, organism = "hsa",keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1)
temp = ego.KEGG[ego.KEGG$p.adjust<=1,]
dim(temp)	#	22  9
outxlsx=paste("cluster",i,"enrichKEGG.xlsx",sep="_")
outtxt=paste("cluster",i,"enrichKEGG.txt",sep="_")
outpdf=paste("cluster",i,"enrichKEGG.pdf",sep="_")
outFile(temp,outxlsx,outtxt,outpdf)

}
data=read.table("cluster_2_enrichKEGG.txt",sep="\t",header=T,quote ="",stringsAsFactors=F)
data[,11]=-log10(data[,5])
data1=data[1:10,]
p1=ggplot(data=data1, aes(x=reorder(Description,V11),y=V11)) + geom_bar(stat="identity", width=0.6,fill="#DEE8D6")+ coord_flip() +  xlab("KEGG Pathway") + ylab("-log10(p)") + theme_bw()
pdf("fig2d.pdf")
p1
dev.off()
#########################figure 2e###########################
#########################figure 2e###########################
library("Seurat")
load("SUM159.Rdata")
library("genefu")
str(pam50)
library("Seurat")
data=read.table("SUM159_id_symbol_sort.txt",sep="\t")
exp=SUM159@assays$RNA@counts[data[,1],] 
pam50genes=pam50$centroids.map[c(1,3)]
intersect(rownames(exp),rownames(pam50genes)) # 44个基因
exp<-t(exp)
colnames(data)=c("probe","Gene.Symbol","EntrezGene.ID")
data=as.data.frame(data)
SubPred_pam50_SUM159<-molecular.subtyping(sbt.model = "pam50",data = exp,annot = data,do.mapping = T)
table(SubPred_pam50_SUM159$subtype)
num=as.matrix(table(as.data.frame(cbind(as.vector(SubPred_pam50_SUM159$subtype),as.vector(Idents(SUM159))))))
anno=cbind(colnames(SUM159),as.vector(SubPred_pam50_SUM159$subtype))
write.table(anno,"pam50_SUM159_anno.txt",quote=F,sep="\t",row.names=F,col.names=F)
num.prop=prop.table(num,m=2)
library("reshape")
colnames(num.prop)=c("DMSO","JQ1","RDMSO","RJQ1")
cell.prop=num.prop
cell.prop.melt=melt(as.matrix(cell.prop))
cell.prop.melt=data.frame(subtype=cell.prop.melt[,1],condition =cell.prop.melt[,2],value=cell.prop.melt[,3])
col_vector=c("#D8BA51","#F8C3B4","#2D81AE","#109880","#B1C1D4")
library("ggplot2")
p1=ggplot(cell.prop.melt,mapping = aes(condition,value,fill=subtype))+geom_bar(stat='identity',position='fill') +labs(x = 'condition',y = 'Fraction') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))+theme(axis.text.x = element_text(angle=0, hjust = 0.5,vjust=1))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual(values=col_vector)

p1=p1+theme(panel.background = element_rect(fill="transparent",colour=NA))
ggsave("fig2e.pdf",p1)

#########################figure 2f###########################
#########################figure 2f###########################
load("SUM159.Rdata")
library("reshape")
library(Seurat)
library(ggplot2)
g2m_genes=cc.genes$g2m.genes
g2m_genes=CaseMatch(search=g2m_genes,match=rownames(SUM159))
s_genes=cc.genes$s.genes
s_genes=CaseMatch(search=s_genes,match=rownames(SUM159))
score=CellCycleScoring(SUM159,g2m.features=g2m_genes,s.features=s_genes)
cell=as.matrix(table(score$Phase,Idents(score)))
colnames(cell)=c("DMSO","JQ1","RDMSO","RJQ1")
cell.prop=prop.table(cell,m=2)
cell.prop.melt=melt(cell.prop)
cell.prop.melt=data.frame(cell_cycle=cell.prop.melt[,1],condition =cell.prop.melt[,2],value=cell.prop.melt[,3])
col_vector=c("#DB7B67","#82D8A9","#8064A2")
p1=ggplot(cell.prop.melt,mapping = aes(condition,value,fill=cell_cycle))+geom_bar(stat='identity',position='fill') +labs(x = 'condition',y = 'Fraction') +theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))+theme(axis.text.x = element_text(angle=0, hjust = 0.5,vjust=1))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual(values=col_vector)

p1=p1+theme(panel.background = element_rect(fill="transparent",colour=NA))
write.table(cell,"cell_cycle_propotion_SUM159.txt",quote=F,sep="\t")
ggsave("fig2f.pdf",p1)
Idents(SUM159)=factor(score$Phase,levels=c("G1","G2M","S"))

#########################figure 2g###########################
#########################figure 2g###########################
rm(list=ls())
library(monocle)
library(Seurat)
#library(SeuratWrappers)
library(ggplot2)
require(dplyr)
library("openxlsx")
load("SUM159.Rdata")
subset_data = SUM159
new.metadata <- merge(subset_data@meta.data,data.frame(Idents(subset_data)),by = "row.names",sort = FALSE)
rownames(new.metadata)<-new.metadata[,1]
head(subset_data@meta.data)
subset_data@meta.data<-new.metadata
expression_matrix <- as(as.matrix(subset_data@assays$RNA@counts), 'sparseMatrix')
cell.type<-Idents(subset_data)
cell_metadata <- new('AnnotatedDataFrame',data=cbind(subset_data@meta.data,cell.type))# 132510
gene_annotation <- new('AnnotatedDataFrame',data=data.frame(gene_short_name = row.names(subset_data), row.names = row.names(subset_data)))
monocle_cds <- newCellDataSet(expression_matrix,phenoData = cell_metadata,featureData = gene_annotation,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())		      
cds <- monocle_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)  ## Removing 99 outliers
diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~ cell.type")
### inference the pseudotrajectory########################################################
# step1: select genes for orderding setOrderingFilter() #
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)  
plot_ordering_genes(cds) #mark the genes ###
# step2: dimension reduction=> reduceDimension()  DDRTree #
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
# step3: ordering the cells=> orderCells()
cds <- orderCells(cds)

pdf("fig2g_1.pdf")
plot_cell_trajectory(cds, color_by = "cell.type",show_branch_points=F) + scale_color_manual(breaks = c("SUM159DMSO", "SUM159JQ1", "SUM159RDMSO","SUM159RJQ1"), values=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4"))   + NoLegend()
dev.off()
cds <- orderCells(cds,root_state=4)
pdf("fig2g-2.pdf")
plot_cell_trajectory(cds, color_by = "Pseudotime")+scale_color_gradientn(colours=c("#440154","#3C508B","#228E8D","#53C569","#EEE51C")) 
dev.off()
