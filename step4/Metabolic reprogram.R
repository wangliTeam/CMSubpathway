library(Seurat)
library(patchwork)
library(dplyr)
library(data.table)
library(stringr)
library(data.table)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
gse246231<-readRDS("./step4/data/GSE246231_gene_log2.rds")
subpathway<-c("ADH1A","ADH1B","ADH1C","ADH4","ADH5","ADH6","ADH7",
              "ALDH1B1","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1")
geneset<-list()
geneset[["pathway"]]<-subpathway
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
cellscore<-gsva(as.matrix(gse246231), geneset,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
cellscore<-as.data.frame(t(cellscore))
cellscore$group<-c(rep("chow",6),rep("HFD",6),rep("Glucose",5))
save(cellscore,file = "./GSE246231_gsva.Rdata")

p<-ggplot(data = cellscore, aes(x = group, y = pathway, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("sub-pathway Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.7)
ggsave("./pathway_score.pdf",p,width = 4,height = 4)

my_comparisons = list(c("HFD","chow"),c("chow","Glucose"),c("Glucose","HFD"))
p<-ggplot(data = cellscore, aes(x = group, y = pathway, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("sub-pathway Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.7)
p<-ggboxplot(cellscore, x = "group", y = "pathway",
             color = "group", palette = "jco",add = "jitter",
             short.panel.labs = FALSE)+
  xlab("") +ylab("sub-pathway Score")+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.8)+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")
ggsave("./pathway_score_p.pdf",p,width = 4,height = 4)

geneexpr<-cbind(cellscore,t(na.omit(gse246231[subpathway,rownames(cellscore)])))
meltgene<-reshape2::melt(geneexpr[,5:12],id.var="group")
p<-ggplot(data = meltgene, aes(x = variable, y = value, fill = group))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.signif",label.y = 10)
ggsave("./pathwaygene_expression.pdf",p,width = 8,height = 4)

my_comparisons = list(c("HFD","chow"),c("chow","Glucose"),c("Glucose","HFD"))
colnames(meltgene)[2]<-"sub-pathway"
p<-ggboxplot(meltgene, x = "group", y = "value",
             color = "group", palette = "jco",add = "jitter",
             facet.by = "sub-pathway", short.panel.labs = FALSE)+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 5)+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")
ggsave("./pathwaygene_expression_pvalue.pdf",p,width = 10,height = 10)

food<-unique(cellscore$group)
plot.list<-list()
for(i in 1:length(food)){
  a<-cellscore[which(cellscore$group==food[i]),]
  p<-ggplot(data=a, aes(x=pathway, y=tumor))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab("Tumor Score")+ggtitle(food[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("./GSE246231_pathway_tumor_cor.pdf",width = 9,height = 3)
wrap_plots(plot.list,ncol = 3)
dev.off()
