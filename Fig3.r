#########################figure 3a###########################
#########################figure 3a###########################
library(ggplot2)
library(ggrepel)
load("SUM159.Rdata")
data=subset(SUM159,idents=c("SUM159DMSO","SUM159RDMSO"))
SUM159RDMSO_DMSO_marker=FindMarkers(SUM159, ident.1="SUM159RDMSO",ident.2="SUM159DMSO",min.pct=0,logfc.threshold=0)
data=SUM159RDMSO_DMSO_marker
data$symbol=rownames(data)
PValue=0.01
FC=0.25
data$sig[(-1*log10(data$p_val_adj) < -1*log10(PValue))|(data$avg_log2FC < FC)& data$avg_log2FC > -FC] <- "NotSig"
data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC >= FC] <- "Up"
data$sig[-1*log10(data$p_val_adj) >= -1*log10(PValue) & data$avg_log2FC <= -FC] <- "Down"
data[order(data[,2],decreasing=T),][1:10,]
data[order(data[,2],decreasing=T),][17367:17377,]
data=data[order(data[,2],decreasing=T),]
data$Marker =0
gene=read.table("volcaplot.gene",sep="\t")
data$Marker[is.element(data$symbol,gene[,1])]=1
data$label=ifelse(data$Marker == 1, as.character(data$symbol), '')
color=c("#234B44","#A6A873","#E28A64")
p=ggplot(data,aes(data$avg_log2FC,-1*log10(data$p_val_adj))) +   
  geom_point(aes(color = sig)) +                           
  labs(title="volcanoplot",                                
       x="log2(FC)", 
       y="-log10(p_val_adj)") + 
  scale_color_manual(values = color) + 
  geom_hline(yintercept=-log10(PValue),linetype=2)+        
  geom_vline(xintercept=c(-FC,FC),linetype=2)+ 
  theme(legend.background=element_rect(fill="transparent"),axis.line = element_line(color = "black"),panel.background=element_rect(fill="transparent"))+
  geom_text_repel(aes(x = data$avg_log2FC, y = -1*log10(data$p_val_adj), label=label), max.overlaps = 10000,size=3,
                  box.padding=unit(0.5,'lines'),           
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                   
                  show.legend=FALSE)
pdf("fig3a.pdf")
p
dev.off()