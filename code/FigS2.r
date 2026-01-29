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
data=read.table("cluster_4_enrichKEGG.txt",sep="\t",header=T,quote ="",stringsAsFactors=F)
data[,11]=-log10(data[,5])
data2=data[1:10,]
p=ggplot(data=data2, aes(x=reorder(Description,V11),y=V11)) + geom_bar(stat="identity", width=0.6,fill="#A1B2CC")+ coord_flip() +  xlab("KEGG Pathway") + ylab("-log10(p)") + theme_bw()
pdf("figs2.pdf")
p
dev.off()