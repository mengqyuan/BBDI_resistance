library(ggpubr)
score=read.table("SUM159CytoTRACE_plot_table.txt",header=T,sep="\t")
p=ggboxplot(
  score,
  x = "Phenotype",
  y = "CytoTRACE",
  color = "black",
  fill = "Phenotype",
  palette=c("#8E288B","#60C5E1","#E38E40","#D86FA7","#56AD8F"),
  xlab = "Phenotype",
  ylab = "CytoTRACE Score",
  main = "CytoTRACE")  +theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))
pdf("figs6.pdf")
p
dev.off()