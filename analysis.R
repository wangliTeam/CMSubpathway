setwd("../metabolic_pathway/")
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
load("../metabolic_pathway/scRNA/allcell/pbmc.harmony.rda")
Idents(pbmc.harmony)<-"celltype"
deg_allcell<-FindAllMarkers(pbmc.harmony,min.pct = 0.25, logfc.threshold = 0.25)
load("../metabolic_pathway/scRNA/epi/pbmc_epi.rda")
Idents(pbmc.harmony)<-"group"
deg_epi<-FindAllMarkers(pbmc.harmony,min.pct = 0.25, logfc.threshold = 0.25)
load("../metabolic_pathway/scRNA/t/pbmc.harmony.rda")
load("../metabolic_pathway/scRNA/t/score.rda")
scRNAsub<-subset(pbmc.harmony,cells = rownames(score))
scRNAsub@meta.data$cluster<-score[rownames(scRNAsub@meta.data),"cluster"]
Idents(scRNAsub)<-"cluster"
deg_CD8T<-FindAllMarkers(scRNAsub,min.pct = 0.25, logfc.threshold = 0.25)
save(deg_allcell,deg_CD8T,deg_epi,file = "scRNA/DEG.Rdata")

data<-GetAssayData(scRNAsub,layer = "data")
p<-DotPlot(scRNAsub,features = "TCF7",group.by = "cluster")
ggsave("scRNA/TCF1.pdf",p,width = 3,height = 5)

p<-DotPlot(scRNAsub,features = "TOX",group.by = "cluster")
score$TCF7<-data["TCF7",rownames(score)]
metadatasc<-score[which(score$TCF7>0),]
p<-ggplot(metadatasc,aes(x=cluster,y=TCF7))+
  geom_boxplot(width=0.4,outlier.shape = NA,aes(color=cluster),
               position = position_dodge(0.5))+
  theme_bw()+xlab("")+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1),
        panel.grid = element_blank())+
  stat_compare_means(aes(group=cluster))

score$PDCD1<-data["PDCD1",rownames(score)]
metadatasc<-score[which(score$PDCD1>0),]
p<-ggplot(metadatasc,aes(x=cluster,y=PDCD1))+
  geom_boxplot(width=0.4,outlier.shape = NA,aes(color=cluster),
               position = position_dodge(0.5))+
  theme_bw()+xlab("")+
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1),
        panel.grid = element_blank())+
  stat_compare_means(aes(group=cluster))
#top30<-deg_CD8T.sig%>%group_by(cluster)%>%top_n(n=30,wt = avg_log2FC)
deg_CD8T.sig<-deg_CD8T[which(deg_CD8T$p_val_adj<0.05&deg_CD8T$pct.1>0.25&deg_CD8T$avg_log2FC>0.25),]
deg_allcell.sig<-deg_allcell[which(deg_allcell$p_val_adj<0.05&deg_allcell$pct.1>0.25&deg_allcell$avg_log2FC>0.25),]
deg_epi.sig<-deg_epi[which(deg_epi$p_val_adj<0.05&deg_epi$pct.1>0.25&deg_epi$avg_log2FC>0.25),]
tumor_gene<-intersect(deg_epi.sig$gene[which(deg_epi.sig$cluster=="malig")],
                      deg_allcell.sig$gene[which(deg_allcell.sig$cluster=="Epithelial")])
Tcell0_gene<-intersect(deg_CD8T.sig$gene[which(deg_CD8T.sig$cluster==0)],
                       deg_allcell.sig$gene[which(deg_allcell.sig$cluster=="Tc")])
Tcell7_gene<-intersect(deg_CD8T.sig$gene[which(deg_CD8T.sig$cluster==7)],
                       deg_allcell.sig$gene[which(deg_allcell.sig$cluster=="Tc")])
gse246231<-readRDS("bulk/GSE246231_gene_log2.rds")
subpathway<-c("ADH1A","ADH1B","ADH1C","ADH4","ADH5","ADH6","ADH7",
              "ALDH1B1","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1")
geneset<-list()
geneset[["pathway"]]<-subpathway
geneset[["tumor"]]<-tumor_gene
geneset[["Tcell0"]]<-Tcell0_gene
geneset[["Tcell7"]]<-Tcell7_gene
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
cellscore<-gsva(as.matrix(gse246231), geneset,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
cellscore<-as.data.frame(t(cellscore))
cellscore$group<-c(rep("chow",6),rep("HFD",6),rep("Glucose",5))
save(cellscore,file = "bulk/GSE246231_gsva.Rdata")

p<-ggplot(data = cellscore, aes(x = group, y = pathway, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("sub-pathway Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.7)
ggsave("bulk/pathway_score.pdf",p,width = 4,height = 4)

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
ggsave("bulk/pathway_score_p.pdf",p,width = 4,height = 4)

p<-ggplot(data = cellscore, aes(x = group, y = Tcell0, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("Cell Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.7)
ggsave("bulk/Tcell0_score.pdf",p,width = 4,height = 4)
p<-ggboxplot(cellscore, x = "group", y = "Tcell0",
             color = "group", palette = "jco",add = "jitter",
             short.panel.labs = FALSE)+
  xlab("") +ylab("CD8 T Cell Score") +
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.7)+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")
ggsave("bulk/Tcell0_score_p.pdf",p,width = 4,height = 4)

p<-ggplot(data = cellscore, aes(x = group, y = Tcell7, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("Cell Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.7)
ggsave("bulk/Tcell7_score.pdf",p,width = 4,height = 4)
p<-ggplot(data = cellscore, aes(x = group, y = tumor, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("Cell Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 1.2)
ggsave("bulk/tumor_score.pdf",p,width = 4,height = 4)

my_comparisons = list(c("HFD","chow"),c("chow","Glucose"),c("Glucose","HFD"))
p<-ggplot(data = cellscore, aes(x = group, y = tumor, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  geom_jitter(shape=16, position = position_jitter(0.2))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("Cell Score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 1.4)
p<-ggboxplot(cellscore, x = "group", y = "tumor",
             color = "group", palette = "jco",add = "jitter",
             short.panel.labs = FALSE)+
  xlab("") +ylab("Tumor Cell Score") +
  stat_compare_means(aes(group=group),label = "p.format",label.y = 1.5)+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")
ggsave("bulk/tumor_score_p.pdf",p,width = 4,height = 4)

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
ggsave("bulk/pathwaygene_expression.pdf",p,width = 8,height = 4)

my_comparisons = list(c("HFD","chow"),c("chow","Glucose"),c("Glucose","HFD"))
colnames(meltgene)[2]<-"sub-pathway"
p<-ggboxplot(meltgene, x = "group", y = "value",
               color = "group", palette = "jco",add = "jitter",
               facet.by = "sub-pathway", short.panel.labs = FALSE)+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 5)+
  stat_compare_means(comparisons = my_comparisons,label = "p.format")
ggsave("bulk/pathwaygene_expression_pvalue.pdf",p,width = 10,height = 10)



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
pdf("bulk/GSE246231_pathway_tumor_cor.pdf",width = 9,height = 3)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(food)){
  a<-cellscore[which(cellscore$group==food[i]),]
  p<-ggplot(data=a, aes(x=pathway, y=Tcell0))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab("Tcell Score")+ggtitle(food[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/GSE246231_pathway_Tcell0_cor.pdf",width = 9,height = 3)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(food)){
  a<-cellscore[which(cellscore$group==food[i]),]
  p<-ggplot(data=a, aes(x=tumor, y=Tcell0))+
    geom_point(color="red")+xlab("Tumor score")+
    ylab("Tcell Score")+ggtitle(food[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/GSE246231_Tumor_Tcell0_cor.pdf",width = 9,height = 3)
wrap_plots(plot.list,ncol = 3)
dev.off()

p<-ggplot(data=a, aes(x=Tcell0, y=tumor))+
  geom_point(color="red")+xlab("Tcell score")+
  ylab("Tumor Score")+ggtitle(dataset[i])+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")


#CTC
load("../metabolic_pathway/bulk/CTC17dataset_exorssion.Rdata")
load("../metabolic_pathway/bulk/CTC17dataset_group.Rdata")
dataset<-c("GSE111842","GSE109761","GSE111065",
           "GSE51827","GSE55807","GSE67939",
           "GSE75367","GSE86978","GSE41245")
score.list<-list()
for(i in 1:length(dataset)){
  expr<-ctc.expression[[dataset[i]]]
  ctcgroup<-ctc.group[[dataset[i]]]
  gsvaobj<-ssgseaParam(as.matrix(expr),geneset,normalize = T)
  gsvascore<-gsva(gsvaobj)
  gsvascore<-as.data.frame(t(gsvascore))
  gsvascore$group<-NA
  gsvascore[ctcgroup$sample,"group"]<-as.character(ctcgroup$group)
  gsvascore<-cbind(gsvascore,t(na.omit(expr[subpathway,rownames(gsvascore)])))
  score.list[[dataset[i]]]<-gsvascore
}

scorep<-do.call(rbind,lapply(score.list,function(x){return(x[,1:5])}))
scorep$dataset<-str_split(rownames(scorep),"\\.",simplify = T)[,1]
p<-ggplot(data = scorep, aes(x = dataset, y = pathway, fill = group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
#  geom_jitter(shape=16, position = position_jitter(0.5))+
  scale_fill_manual(values = c("red","blue","orange"))+
  xlab("") +ylab("sub-pathway score") +theme_classic() + 
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.signif",label.y = 1)
ggsave("bulk/CTCdata_pathwayscore.pdf",p,width = 7,height = 3)

plot.list<-list()
for(i in 1:length(dataset)){
  a<-score.list[[dataset[i]]]
  a<-a[which(a$group=="cancer"),]
  p<-ggplot(data=a, aes(x=pathway, y=tumor))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab("Tumor Score")+ggtitle(dataset[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/CTCdata_pathwaytumor_cor.pdf",width = 10,height = 10)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(dataset)){
  a<-score.list[[dataset[i]]]
  a<-a[which(a$group=="cancer"),]
  p<-ggplot(data=a, aes(x=Tcell0, y=tumor))+
    geom_point(color="red")+xlab("Tcell score")+
    ylab("Tumor Score")+ggtitle(dataset[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/CTCdata_Tcell0tumor_cor.pdf",width = 10,height = 10)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(dataset)){
  a<-score.list[[dataset[i]]]
  a<-a[which(a$group=="cancer"),]
  p<-ggplot(data=a, aes(x=Tcell7, y=tumor))+
    geom_point(color="red")+xlab("Tcell score")+
    ylab("Tumor Score")+ggtitle(dataset[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/CTCdata_Tcell7tumor_cor.pdf",width = 10,height = 10)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(dataset)){
  a<-score.list[[dataset[i]]]
  a<-a[which(a$group=="cancer"),]
  p<-ggplot(data=a, aes(x=pathway, y=Tcell0))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab("Tcell Score")+ggtitle(dataset[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/CTCdata_pathwayTcell0_cor.pdf",width = 10,height = 10)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(dataset)){
  a<-score.list[[dataset[i]]]
  a<-a[which(a$group=="cancer"),]
  p<-ggplot(data=a, aes(x=pathway, y=Tcell7))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab("Tcell Score")+ggtitle(dataset[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/CTCdata_pathwayTcell7_cor.pdf",width = 10,height = 10)
wrap_plots(plot.list,ncol = 3)
dev.off()

plot.list<-list()
for(i in 1:length(dataset)){
  a<-score.list[[dataset[i]]]
  meltgene<-reshape2::melt(a[,5:17],id.var="group")
  p<-ggplot(data = meltgene, aes(x = variable, y = value))+
    geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
    scale_color_manual(values = c("blue","red"))+
    xlab("") +ylab("Expression") +ggtitle(dataset[i])+
    theme_classic() + 
    theme(panel.grid.major = element_blank()) +
    theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
    theme(legend.key=element_blank())+
    stat_compare_means(aes(group=group),label = "p.signif",label.y = (max(meltgene$value)-1))
  plot.list[[i]]<-p
}
pdf("bulk/CTCdata_pathwaygene_expression.pdf",width = 12,height = 15)
wrap_plots(plot.list,ncol = 2)
dev.off()

pvalue<-matrix(0,12,9,dimnames = list(subpathway,dataset))
for(i in 1:length(dataset)){
  a<-ctc.expression[[dataset[i]]]
  a1<-ctc.group[[dataset[i]]]
  a1<-cbind(a1,t(na.omit(a[subpathway,a1$sample])))
  p<-sapply(subpathway,function(x){
    return(wilcox.test(as.formula(paste0(x,"~group")),a1)$p.value)
  })
  p<-signif(p,digits = 2)
  pvalue[,i]<-p
}
pvalue<-as.data.frame(pvalue)
write.csv(pvalue,file = "bulk/CTCpvalue.csv")

#单细胞数据中刻画恶性细胞和正常细胞间，通路相关基因的差异表达
load("../metabolic_pathway/scRNA/epi/pbmc_epi.rda")
subpathway<-c("ADH1A","ADH1B","ADH1C","ADH4","ADH5","ADH6","ADH7",
              "ALDH1B1","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1")
metadata<-pbmc.harmony@meta.data
scdata<-GetAssayData(pbmc.harmony,layer = "data")
metadata<-cbind(metadata,t(na.omit(scdata[subpathway[subpathway%in%rownames(scdata)],rownames(metadata)])))
genemelt<-reshape2::melt(metadata[,c(7,10:19)],id.var="group")
p<-ggplot(data = genemelt, aes(x = variable, y = value,color=group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  scale_color_manual(values = c("blue","red"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(0,2))+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.signif",label.y = 1.8)
ggsave("scRNA/subpathway_gene.pdf",p,width = 8,height = 3)
p<-ggplot(data = genemelt, aes(x = variable, y = value,color=group))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5))+
  scale_color_manual(values = c("blue","red"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(0,2))+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 1.8)
ggsave("scRNA/subpathway_gene_pvalue.pdf",p,width = 10,height = 3)

genedata<-as.matrix(scdata[subpathway[subpathway%in%rownames(scdata)],])
metadata<-metadata[order(metadata$group),]
annodata<-data.frame(row.names = rownames(metadata),Group=metadata$group)
genedata<-genedata[,rownames(annodata)]
library(pheatmap)
pheatmap(genedata,cluster_rows = T,cluster_cols = F,
         show_rownames = T,show_colnames = F,annotation_col = annodata,
         color=colorRampPalette(c("white","#EC1B21"))(100),
         filename = "scRNA/subpathway_expression_pheatmap.pdf",width = 8,height = 3)
n<-t(scale(t(genedata)))
n<-na.omit(n)
pheatmap(n,cluster_rows = T,cluster_cols = F,
         show_rownames = T,show_colnames = F,annotation_col = annodata,
         breaks = seq(-1,1,length.out=100),
         color=colorRampPalette(c("#4169E1","white","#EC1B21"))(100),
         filename = "scRNA/subpathway_expression_pheatmap2.pdf",width = 8,height = 3)

#TCGA和GEO
load("../metabolic_pathway/bulk/TCGA_exp.rda")
subpathway<-c("ADH1A","ADH1B","ADH1C","ADH4","ADH5","ADH6","ADH7",
              "ALDH1B1","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1")
genemeta<-exp[subpathway,]
genemeta<-as.data.frame(t(genemeta))
genemeta$group<-ifelse(as.numeric(str_sub(rownames(genemeta),14,15))<10,"cancer","normal")
genemelt<-reshape2::melt(genemeta,id.var="group")
p<-ggplot(data = genemelt, aes(x = variable, y = value))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(0,20))+ggtitle("TCGA")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.signif",label.y = 18)
ggsave("bulk/subpathway_gene_TCGA.pdf",p,width = 10,height = 4)
p<-ggplot(data = genemelt, aes(x = variable, y = value))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(0,20))+ggtitle("TCGA")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 18)
ggsave("bulk/subpathway_gene_TCGA_pvalue.pdf",p,width = 12,height = 4)

library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
gsvaobj<-ssgseaParam(as.matrix(exp),geneset,normalize = T)
gsvascore<-gsva(gsvaobj)
gsvascore<-as.data.frame(t(gsvascore))
gsvascore$group<-ifelse(as.numeric(str_sub(rownames(gsvascore),14,15))<10,"cancer","normal")
save(gsvascore,file = "bulk/TCGA_pathway_score.Rdata")
p<-ggplot(data = gsvascore, aes(x = group, y = pathway))+
  geom_boxplot(width=0.8,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  geom_jitter(shape=16, position = position_jitter(0.3),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("sub-pathway score") +theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ggtitle("TCGA")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.8)
ggsave("bulk/subpathway_score_TCGA.pdf",p,width = 4,height = 4)

a<-gsvascore[which(gsvascore$group=="cancer"),]
p<-ggplot(data=a, aes(x=pathway, y=tumor))+
  geom_point(color="red")+xlab("key sub-pathway score")+
  ylab("Tumor Score")+ggtitle("TCGA")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")
ggsave("bulk/subpathway_tumor_cor_TCGA.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=Tcell0, y=pathway))+
  geom_point(color="red")+xlab("CD8 T cell subtype 0")+
  ylab("key sub-pathway score")+ggtitle("TCGA")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")
ggsave("bulk/subpathway_Tcell0_cor_TCGA.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=tumor, y=Tcell0))+
  geom_point(color="red")+xlab("Tumor score")+
  ylab("CD8 T cell subtype 0")+ggtitle("TCGA")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")
ggsave("bulk/Tumor_Tcell0_cor_TCGA.pdf",p,width = 4,height = 4)

gsvascore<-cbind(gsvascore,genemeta[rownames(gsvascore),1:12])
plot.list<-list()
for(i in 1:length(subpathway)){
  p<-ggplot(data=gsvascore, aes_string(x=as.name("pathway"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=gsvascore, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/subpathway_gene_cor_TCGA.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

plot.list<-list()
for(i in 1:length(subpathway)){
  a<-gsvascore[which(gsvascore$group=="cancer"),]
  p<-ggplot(data=a, aes_string(x=as.name("tumor"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("Tumor score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/tumor_gene_cor_TCGA.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

load("../metabolic_pathway/bulk/GEO_anno.rda")
load("../metabolic_pathway/bulk/GEO_exp.rda")
gseanno<-data.frame(sample=colnames(exp),group=anno)
genemeta<-exp[subpathway,]
genemeta<-as.data.frame(t(genemeta))
genemeta$group<-gseanno$group
genemelt<-reshape2::melt(genemeta,id.var="group")
genemelt$group<-factor(genemelt$group,levels = c("tumor","normal"))
p<-ggplot(data = genemelt, aes(x = variable, y = value))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(2,14))+ggtitle("GSE42568")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.signif",label.y = 13.5)
ggsave("bulk/subpathway_gene_GSE42568.pdf",p,width = 10,height = 4)

p<-ggplot(data = genemelt, aes(x = variable, y = value))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(2,14))+ggtitle("GSE42568")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 13.5)
ggsave("bulk/subpathway_gene_GSE42568_pvalue.pdf",p,width = 12,height = 4)

gsvaobj<-ssgseaParam(as.matrix(exp),geneset,normalize = T)
gsvascore<-gsva(gsvaobj)
gsvascore<-as.data.frame(t(gsvascore))
gsvascore$group<-NA
gsvascore[gseanno$sample,"group"]<-gseanno$group
save(gsvascore,file = "bulk/GSE42568_pathway_score.Rdata")
gsvascore$group<-factor(gsvascore$group,levels = c("tumor","normal"))
p<-ggplot(data = gsvascore, aes(x = group, y = pathway))+
  geom_boxplot(width=0.8,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  geom_jitter(shape=16, position = position_jitter(0.3),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("key sub-pathway score") +theme_classic() + 
  coord_cartesian(ylim = c(0,1))+ggtitle("GSE42568")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 0.9)
ggsave("bulk/subpathway_score_GSE.pdf",p,width = 4,height = 4)

a<-gsvascore[which(gsvascore$group=="tumor"),]
p<-ggplot(data=a, aes(x=pathway, y=tumor))+
  geom_point(color="red")+xlab("key sub-pathway score")+
  ylab("Tumor Score")+ggtitle("GSE42568")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")
ggsave("bulk/subpathway_tumor_cor_GSE42568.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=pathway, y=Tcell0))+
  geom_point(color="red")+xlab("key sub-pathway score")+
  ylab("CD8 T cell subtype 0")+ggtitle("GSE42568")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")
ggsave("bulk/subpathway_Tcell0_cor_GSE42568.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=tumor, y=Tcell0))+
  geom_point(color="red")+xlab("Tumor score")+
  ylab("CD8 T cell subtype 0")+ggtitle("GSE42568")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#去除背景色
  theme(panel.grid=element_blank())+#去除网格线
  stat_cor(data=a, method = "spearman")
ggsave("bulk/Tumor_Tcell0_cor_GSE42568.pdf",p,width = 4,height = 4)

gsvascore<-cbind(gsvascore,genemeta[rownames(gsvascore),1:12])
plot.list<-list()
for(i in 1:length(subpathway)){
  p<-ggplot(data=gsvascore, aes_string(x=as.name("pathway"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=gsvascore, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/subpathway_gene_cor_GSE42568.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

plot.list<-list()
for(i in 1:length(subpathway)){
  a<-gsvascore[which(gsvascore$group=="tumor"),]
  p<-ggplot(data=a, aes_string(x=as.name("tumor"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("Tumor score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#去除背景色
    theme(panel.grid=element_blank())+#去除网格线
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("bulk/tumor_gene_cor_GSE42568.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

#
chen_genelist <- c("X","ADH1A","ADH1B","ADH1C","ADH4","ADH5","ADH6","ADH7","ALDH1B1","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1")
CRISPRGeneEffect_chen <- CRISPRGeneEffect[,chen_genelist]

OncotreeLineage <- Model[which(Model$OncotreeLineage == "Breast"),]
Model_chen <- OncotreeLineage[which(OncotreeLineage$OncotreePrimaryDisease != "Non-Cancerous"),]

CRISPRGeneEffect_chen_breast <- CRISPRGeneEffect_chen[which(CRISPRGeneEffect_chen$X %in% Model_chen$ModelID),]


setwd("../")
plot_list <- list()
for(i in 2:13){
  plot_data <- CRISPRGeneEffect_chen_breast[,c(1,i)] %>%           
    arrange(.[[2]])                      
  soft_pink <- "#E9746A" 
  y_var_name <- colnames(plot_data)[2]
  
  p <- ggplot(plot_data, aes(x = reorder(X, .data[[y_var_name]]), y = .data[[y_var_name]])) +
    
    geom_bar(stat = "identity", fill = soft_pink, color = "white", linewidth = 0.2) +

    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +

    labs(title = paste("Gene Effect of",y_var_name),
         x = "Cell Lines",
         y = "Gene Effect Score") +

    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5), 
      panel.grid.major.x = element_blank() 
    )
    plot_list[[i - 1]] <- p
}

combined_plot <- wrap_plots(plot_list[c(5,8,11)], ncol = 1, nrow = 3)

ggsave("Combined_3_Genes_Waterfall_Breast_cancer.pdf", 
       plot = combined_plot, 
       width = 6, 
       height = 8, 
       limitsize = FALSE) 

percentage_threshold <- 0.75
data_matrix <- as.matrix(CRISPRGeneEffect_chen_breast[, -1])

gene_names <- colnames(data_matrix)

count_positive_per_gene <- colSums(data_matrix > 0, na.rm = TRUE)
non_na_counts_per_gene <- colSums(!is.na(data_matrix))
proportion_positive_per_gene <- count_positive_per_gene / non_na_counts_per_gene
genes_highly_positive <- gene_names[proportion_positive_per_gene >= percentage_threshold]
count_negative_per_gene <- colSums(data_matrix < 0, na.rm = TRUE)
proportion_negative_per_gene <- count_negative_per_gene / non_na_counts_per_gene
genes_highly_negative <- gene_names[proportion_negative_per_gene >= percentage_threshold]
df_highly_positive_genes <- CRISPRGeneEffect_chen_breast[, c("X", genes_highly_positive), drop = FALSE]
df_highly_negative_genes <- CRISPRGeneEffect_chen_breast[, c("X", genes_highly_negative), drop = FALSE]
all_filtered_genes <- unique(c(genes_highly_positive, genes_highly_negative))
df_all_filtered_genes <- CRISPRGeneEffect_chen_breast[, c("X", all_filtered_genes), drop = FALSE]
