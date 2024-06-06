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
p1=ggplot(score1,aes(res,mean_sihouette))+
  geom_point(size=4,color="#FF9900")+
  geom_line(position = position_dodge(0.1),cex=1.3,color="#3D7E9C")+
  scale_y_continuous(limits = c(0.1,0.25),expand = c(0,0))+scale_x_continuous(breaks = seq(0.1,1,0.1))+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())
pdf("figs5.pdf")
p1
dev.off()

