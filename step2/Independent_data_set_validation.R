#library(GEOquery)
library(GSVA)
library(survival)
library(survminer)
library(pROC)
library(e1071)
library(genefilter)
library(Biobase)
library(stringr)
#TCGA  and GEO
#TCGA
load("./step1/data/TCGA_exp.rda")
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
ggsave("./subpathway_gene_TCGA.pdf",p,width = 10,height = 4)
p<-ggplot(data = genemelt, aes(x = variable, y = value))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(0,20))+ggtitle("TCGA")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 18)
ggsave("./subpathway_gene_TCGA_pvalue.pdf",p,width = 12,height = 4)

gsvaobj<-ssgseaParam(as.matrix(exp),geneset,normalize = T)
gsvascore<-gsva(gsvaobj)
gsvascore<-as.data.frame(t(gsvascore))
gsvascore$group<-ifelse(as.numeric(str_sub(rownames(gsvascore),14,15))<10,"cancer","normal")
save(gsvascore,file = "./TCGA_pathway_score.Rdata")
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
ggsave("./subpathway_score_TCGA.pdf",p,width = 4,height = 4)

a<-gsvascore[which(gsvascore$group=="cancer"),]
p<-ggplot(data=a, aes(x=pathway, y=tumor))+
  geom_point(color="red")+xlab("key sub-pathway score")+
  ylab("Tumor Score")+ggtitle("TCGA")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#
  theme(panel.grid=element_blank())+#
  stat_cor(data=a, method = "spearman")
ggsave("./subpathway_tumor_cor_TCGA.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=Tcell0, y=pathway))+
  geom_point(color="red")+xlab("CD8 T cell subtype 0")+
  ylab("key sub-pathway score")+ggtitle("TCGA")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#
  theme(panel.grid=element_blank())+#
  stat_cor(data=a, method = "spearman")
ggsave("./subpathway_Tcell0_cor_TCGA.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=tumor, y=Tcell0))+
  geom_point(color="red")+xlab("Tumor score")+
  ylab("CD8 T cell subtype 0")+ggtitle("TCGA")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#
  theme(panel.grid=element_blank())+#
  stat_cor(data=a, method = "spearman")
ggsave("./Tumor_Tcell0_cor_TCGA.pdf",p,width = 4,height = 4)

gsvascore<-cbind(gsvascore,genemeta[rownames(gsvascore),1:12])
plot.list<-list()
for(i in 1:length(subpathway)){
  p<-ggplot(data=gsvascore, aes_string(x=as.name("pathway"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#
    theme(panel.grid=element_blank())+#
    stat_cor(data=gsvascore, method = "spearman")
  plot.list[[i]]<-p
}
pdf("./subpathway_gene_cor_TCGA.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

plot.list<-list()
for(i in 1:length(subpathway)){
  a<-gsvascore[which(gsvascore$group=="cancer"),]
  p<-ggplot(data=a, aes_string(x=as.name("tumor"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("Tumor score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#
    theme(panel.grid=element_blank())+#
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("./tumor_gene_cor_TCGA.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

####GSE42568------------------------------------------------------------------------------------------------
load("./step2/data/GSE42568_anno.rda")
load("./step2/data/GSE42568_exp.rda") 
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
ggsave("./subpathway_gene_GSE42568.pdf",p,width = 10,height = 4)

p<-ggplot(data = genemelt, aes(x = variable, y = value))+
  geom_boxplot(width=0.4,outlier.shape = NA,position = position_dodge(0.5),aes(color=group))+
  scale_color_manual(values = c("red","blue"))+
  xlab("") +ylab("Expression") +theme_classic() + 
  coord_cartesian(ylim = c(2,14))+ggtitle("GSE42568")+
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1, vjust=1)) + 
  theme(legend.key=element_blank())+
  stat_compare_means(aes(group=group),label = "p.format",label.y = 13.5)
ggsave("./subpathway_gene_GSE42568_pvalue.pdf",p,width = 12,height = 4)
geneset<-list()
geneset[["pathway"]]<- c("ADH1A","ADH1B","ADH1C","ADH4","ADH5","ADH6","ADH7",
              "ALDH1B1","ALDH2","ALDH3A2","ALDH7A1","ALDH9A1")
gsvaobj<-ssgseaParam(as.matrix(exp),geneset,normalize = T)
gsvascore<-gsva(gsvaobj)
gsvascore<-as.data.frame(t(gsvascore))
gsvascore$group<-NA
gsvascore[gseanno$sample,"group"]<-gseanno$group
save(gsvascore,file = "./step2/data/GSE42568_pathway_score.Rdata")
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
ggsave("./subpathway_score_GSE.pdf",p,width = 4,height = 4)

a<-gsvascore[which(gsvascore$group=="tumor"),]
p<-ggplot(data=a, aes(x=pathway, y=tumor))+
  geom_point(color="red")+xlab("key sub-pathway score")+
  ylab("Tumor Score")+ggtitle("GSE42568")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#
  theme(panel.grid=element_blank())+#
  stat_cor(data=a, method = "spearman")
ggsave("./subpathway_tumor_cor_GSE42568.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=pathway, y=Tcell0))+
  geom_point(color="red")+xlab("key sub-pathway score")+
  ylab("CD8 T cell subtype 0")+ggtitle("GSE42568")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#
  theme(panel.grid=element_blank())+#
  stat_cor(data=a, method = "spearman")
ggsave("./subpathway_Tcell0_cor_GSE42568.pdf",p,width = 4,height = 4)

p<-ggplot(data=a, aes(x=tumor, y=Tcell0))+
  geom_point(color="red")+xlab("Tumor score")+
  ylab("CD8 T cell subtype 0")+ggtitle("GSE42568")+
  stat_smooth(method="lm",se=FALSE)+
  theme_bw()+#
  theme(panel.grid=element_blank())+#
  stat_cor(data=a, method = "spearman")
ggsave("./Tumor_Tcell0_cor_GSE42568.pdf",p,width = 4,height = 4)

gsvascore<-cbind(gsvascore,genemeta[rownames(gsvascore),1:12])
plot.list<-list()
for(i in 1:length(subpathway)){
  p<-ggplot(data=gsvascore, aes_string(x=as.name("pathway"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("sub-pathway score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#
    theme(panel.grid=element_blank())+#
    stat_cor(data=gsvascore, method = "spearman")
  plot.list[[i]]<-p
}
pdf("./subpathway_gene_cor_GSE42568.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()

plot.list<-list()
for(i in 1:length(subpathway)){
  a<-gsvascore[which(gsvascore$group=="tumor"),]
  p<-ggplot(data=a, aes_string(x=as.name("tumor"), y=as.name(subpathway[i])))+
    geom_point(color="red")+xlab("Tumor score")+
    ylab(subpathway[i])+ggtitle(subpathway[i])+
    stat_smooth(method="lm",se=FALSE)+
    theme_bw()+#
    theme(panel.grid=element_blank())+#
    stat_cor(data=a, method = "spearman")
  plot.list[[i]]<-p
}
pdf("./tumor_gene_cor_GSE42568.pdf",width = 14,height = 10)
wrap_plots(plot.list,ncol = 4)
dev.off()