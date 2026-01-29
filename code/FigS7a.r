library(ggplot2)
library(reshape2)
library(ggbreak)
data=read.table("filgotinib_jq1_fc.txt",sep="\t",header=T,quote="")
colnames(data)=c("drug","up","down")

temp=data[1:8,][order(data[1:8,3],decreasing=F),]
temp1=data[9:20,][order(data[9:20,3],decreasing=F),]

out=rbind(temp1,temp)
a=cbind(out[,1:2],rep("Filgotinib",20))
colnames(a)=c("gene","value","group")
b=cbind(out[,c(1,3)],rep("JQ1",20))
colnames(b)=c("gene","value","group")
out1=rbind(a,b)
#out[,1]=factor(out[,1],levels = out[,1])
out1[,1]=factor(out1[,1],levels = unique(out1[,1]))


p9<-ggplot(out1,aes(x=gene,ymin=0,ymax=value,group=group)) +
 geom_linerange(position = position_dodge2(.5),linetype="dotdash",linewidth=0.5)+
geom_point(aes(x=gene, y=value,color=group),size=2,position = position_dodge(.5), fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=1) +
scale_color_manual(values=c("#2F9894","#D24625"))+
theme_light()+theme(panel.grid.major.x=element_blank(),
panel.border=element_blank(),
axis.ticks.x=element_blank())+
xlab("")+
ylab("log2FC")+
#scale_y_break(c(2.5,5),space=0.1,scales=0.2)+
coord_flip()
#p9+coord_flip()
p=p9+ theme_bw()+theme(axis.text.x = element_text(angle=30, hjust = 1,vjust=1))
pdf("fig7b.pdf",width=4,height=5)
print(p)
dev.off()