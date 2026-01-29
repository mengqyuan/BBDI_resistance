#########################figure s4c############################
#########################figure s4c###########################
mat=read.table("infercnv.observations.txt",header=T,row.names=1,sep=" ",stringsAsFactors=F)
gene=read.table("volcaplot.gene",sep="\t")
mat_1=mat[gene[,1],]
mat_1=cbind(rownames(mat_1),mat_1)
colnames(mat_1)[1]="Gene"
#length(grep(".5",colnames(mat))) 751
#length(grep(".7",colnames(mat))) 1130
out=c()
filter_gene=c()
for(i in 1:dim(mat_1)[1])
{
        p=wilcox.test(as.numeric(mat_1[i,1:751]),as.numeric(mat_1[i,752:1881]))
	p1=wilcox.test(as.numeric(mat_1[i,1:751]),as.numeric(mat_1[i,752:1881]),alternative="greater")
	p2=wilcox.test(as.numeric(mat_1[i,1:751]),as.numeric(mat_1[i,752:1881]),alternative="less")
        print(rownames(mat_1)[i])
        print(p$p.value)
        if(p$p.value<0.01){
        filter_gene=c(filter_gene,rownames(mat_1)[i])
        }
	out=rbind(out,c(rownames(mat_1)[i],p$p.value,p1$p.value,p2$p.value))
}
new=melt(mat_1)
new1=cbind(new,c(rep("DMSO",20*751),rep("RDMSO",20*1130)))
colnames(new1)[4]="group"
library(ggplot2)
new1[,1]=factor(new[,1],levels=gene[,1])
pdf("figs4c.pdf",width=15,height=6)
ggplot(new1,aes(x=Gene,y=value,fill=group))+geom_boxplot(aes(fill=group),width=0.3,outlier.size=1,position = position_dodge(width = 0.4))+
theme(axis.text.x = element_text(angle=45, hjust = 1,vjust=1))+
scale_fill_manual(values=c("#8491B4","#FDC086"))

dev.off()


#########################figure s4d############################
#########################figure s4d###########################
library(ggpubr)
library("Seurat")
gene=read.table("volcaplot.gene",sep="\t")
genes=gene[,1]
load("SUM149.Rdata")
ge=intersect(genes,rownames(SUM149))
used=SUM149@assays$RNA@data[ge,]
up_index=c(1:8)
down_index=c(9:19)
score=colSums(used[up_index,])-colSums(used[down_index,])
re=as.data.frame(cbind(score,as.vector(Idents(SUM149))))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
a=as.numeric(re[which(re[,2]=="SUM149DMSO"),1])
b=as.numeric(re[which(re[,2]=="SUM149RDMSO"),1])
wilcox.test(a,b,alternative="less") 
pdf("figs4d.pdf",height=5,width=5)
boxplot(score~Group,data=re,col=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4"),outline=F)
dev.off()
#########################figure s4e############################
#########################figure s4e###########################
gene=read.table("volcaplot.gene",sep="\t")
genes=gene[,1]
load("SUM159.Rdata")
used=SUM159@assays$RNA@data[genes,]
score=colSums(used[up_index,])-colSums(used[down_index,])
re=as.data.frame(cbind(score,as.vector(Idents(SUM159))))
colnames(re)=c("score","Group")
re[,1]=as.numeric(re[,1])
SUM159_plot=boxplot(score~Group,data=re,col=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4"),outline=F)
pdf("figs4e.pdf",height=5,width=5)
boxplot(score~Group,data=re,col=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4"),outline=F)
dev.off()