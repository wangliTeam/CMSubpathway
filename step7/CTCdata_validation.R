library(limma)
library(GSVA)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
####CTC datasets------------------------------------------------------------------------------------------------
#12 genes
core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")

library(data.table)
library(openxlsx)
library(stringr)
library(readxl)
ctc.dectail<-read.xlsx("../step7/CTCdataset.xlsx")
wbc<-read.xlsx(".wbc_GSE51984_gene_FPKM.xlsx")
  
  ctc.expression<-list()
  ctx.deg<-list()
  ctc.group<-list()
  library(limma)
  #1 colorectal cancer
  gse<-read.xlsx("./GSE74369_COMPLETE_RAW_FPKM_TABLE.xlsx") #data download from GEO website
  row.names(gse)<-gse[,1]
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion<-annotion[which(!str_detect(annotion$sample,"tissue")),]
  annotion$group<-ifelse(str_detect(annotion$sample,"normal"),"normal","cancer")
  annotion$group[34:36]<-"normal"
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  gse<-gse[,annotion$sample]
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  coad.deg = topTable(fit, n=Inf,sort.by = "P")
  coad.deg$change<-ifelse(coad.deg$P.Value<0.05&abs(coad.deg$logFC)>1,
                           ifelse(coad.deg$logFC>1,"UP","DOWN"),"NOT")
  table(coad.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[2]]<-gsenor
  names(ctc.expression)[2]<-"GSE74369"
  ctx.deg[[2]]<-coad.deg
  names(ctx.deg)[2]<-"GSE74369"
  ctc.group[[2]]<-annotion
  names(ctc.group)[2]<-"GSE74369"

    #2 pancreatic cancer
  gse<-fread(".GSE40174_human_processed_data.txt.gz")
  gse<-as.data.frame(gse)
  row.names(gse)<-gse$ID
  gse<-gse[,-(1:2)]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-ifelse(str_detect(annotion$sample,"HD"),"normal","cancer")
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  paad.deg = topTable(fit, n=Inf,sort.by = "P")
  paad.deg$change<-ifelse(paad.deg$P.Value<0.05&abs(paad.deg$logFC)>1,
                          ifelse(paad.deg$logFC>1,"UP","DOWN"),"NOT")
  table(paad.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[4]]<-gsenor
  names(ctc.expression)[4]<-"GSE40174"
  ctx.deg[[4]]<-paad.deg
  names(ctx.deg)[4]<-"GSE40174"
  ctc.group[[4]]<-annotion
  names(ctc.group)[4]<-"GSE40174"

  #3 non-small cell lung cancer
  gse<-read.xlsx(".no_samll_lu_GSE74639_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:35],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",10),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  lung.deg = topTable(fit, n=Inf,sort.by = "P")
  lung.deg$change<-ifelse(lung.deg$P.Value<0.05&abs(lung.deg$logFC)>1,
                           ifelse(lung.deg$logFC>1,"UP","DOWN"),"NOT")
  table(lung.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[3]]<-gsenor
  names(ctc.expression)[3]<-"GSE74639"
  ctx.deg[[3]]<-lung.deg
  names(ctx.deg)[3]<-"GSE74639"
  ctc.group[[3]]<-annotion
  names(ctc.group)[3]<-"GSE74639"
  
  #4 prostate cancer GSE67980
  gse<-fread(".GSE67980_readCounts.txt.gz")
  gse<-as.data.frame(gse)
  sample<-fread(".GSE67980_sampleProperties.txt.gz")
  gse<-gse[,c("symbol",sample$title[c(46:122,167:169)])]
  gse<-aggregate(x=gse[2:81],by=list(gse$symbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  gsenor<-normalizeBetweenArrays(gse)
  gsenor<-as.data.frame(gsenor)
  annotion<-data.frame(sample=colnames(gsenor),group=NA)
  annotion$group<-ifelse(str_detect(annotion$sample,"HD"),"normal","cancer")
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  prad.deg = topTable(fit, n=Inf,sort.by = "P")
  prad.deg$change<-ifelse(prad.deg$P.Value<0.05&abs(prad.deg$logFC)>1,
                          ifelse(prad.deg$logFC>1,"UP","DOWN"),"NOT")
  table(prad.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[5]]<-gsenor
  names(ctc.expression)[5]<-"GSE67980"
  ctx.deg[[5]]<-prad.deg
  names(ctx.deg)[5]<-"GSE67980"
  ctc.group[[5]]<-annotion
  names(ctc.group)[5]<-"GSE67980"
  
  ##5 prostate cancer GSE104209  对于只有CTC样本的数据，使用WBC作为对照
  gse<-read.xlsx(".Pro_GSE104209_gene_FPKM.xlsx")
  wbc<-read.xlsx(".wbc_GSE51984_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:37],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",12),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  prad2.deg = topTable(fit, n=Inf,sort.by = "P")
  prad2.deg$change<-ifelse(prad2.deg$P.Value<0.05&abs(prad2.deg$logFC)>1,
                          ifelse(prad2.deg$logFC>1,"UP","DOWN"),"NOT")
  table(prad2.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[6]]<-gsenor
  names(ctc.expression)[6]<-"GSE104209"
  ctx.deg[[6]]<-prad2.deg
  names(ctx.deg)[6]<-"GSE104209"
  ctc.group[[6]]<-annotion
  names(ctc.group)[6]<-"GSE104209"
  
    #6 liver cancer
  gse<-read.xlsx(".LI_GSE117623_gene.xlsx")
  gse<-merge(gse,wbc,by.x = "gene",by.y  = "EnsemblGene_GeneSymbol")
  gse$gene<-str_split(gse$gene,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:70],by=list(gse$gene),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",45),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  liver.deg = topTable(fit, n=Inf,sort.by = "P")
  liver.deg$change<-ifelse(liver.deg$P.Value<0.05&abs(liver.deg$logFC)>1,
                           ifelse(liver.deg$logFC>1,"UP","DOWN"),"NOT")
  table(liver.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[1]]<-gsenor
  names(ctc.expression)[1]<-"GSE117623"
  ctx.deg[[1]]<-liver.deg
  names(ctx.deg)[1]<-"GSE117623"
  ctc.group[[1]]<-annotion
  names(ctc.group)[1]<-"GSE117623"

  #7 pancreatic cancer GSE60407
  gse<-read.xlsx(".panc_GSE60407_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:32],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",7),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  paad2.deg = topTable(fit, n=Inf,sort.by = "P")
  paad2.deg$change<-ifelse(paad2.deg$P.Value<0.05&abs(paad2.deg$logFC)>1,
                           ifelse(paad2.deg$logFC>1,"UP","DOWN"),"NOT")
  table(paad2.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[7]]<-gsenor
  names(ctc.expression)[7]<-"GSE60407"
  ctx.deg[[7]]<-paad2.deg
  names(ctx.deg)[7]<-"GSE60407"
  ctc.group[[7]]<-annotion
  names(ctc.group)[7]<-"GSE60407"
  
  #8 melanoma GSE38495
  gse<-read.xlsx(".mela_GSE38495_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:31],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",6),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  skcm.deg = topTable(fit, n=Inf,sort.by = "P")
  skcm.deg$change<-ifelse(skcm.deg$P.Value<0.05&abs(skcm.deg$logFC)>1,
                          ifelse(skcm.deg$logFC>1,"UP","DOWN"),"NOT")
  table(skcm.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[8]]<-gsenor
  names(ctc.expression)[8]<-"GSE38495"
  ctx.deg[[8]]<-skcm.deg
  names(ctx.deg)[8]<-"GSE38495"
  ctc.group[[8]]<-annotion
  names(ctc.group)[8]<-"GSE38495"
  
  
  #9 breast cancer GSE111842
  gse<-read.xlsx(".BR_GSE111842_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:38],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",13),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca1.deg = topTable(fit, n=Inf,sort.by = "P")
  brca1.deg$change<-ifelse(brca1.deg$P.Value<0.05&abs(brca1.deg$logFC)>1,
                           ifelse(brca1.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca1.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[9]]<-gsenor
  names(ctc.expression)[9]<-"GSE111842"
  ctx.deg[[9]]<-brca1.deg
  names(ctx.deg)[9]<-"GSE111842"
  ctc.group[[9]]<-annotion
  names(ctc.group)[9]<-"GSE111842"
  
  #10 breast cancer GSE109761
  gse<-read.xlsx(".BR_GSE109761_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:84],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",59),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca2.deg = topTable(fit, n=Inf,sort.by = "P")
  brca2.deg$change<-ifelse(brca2.deg$P.Value<0.05&abs(brca2.deg$logFC)>1,
                           ifelse(brca2.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca2.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[10]]<-gsenor
  names(ctc.expression)[10]<-"GSE109761"
  ctx.deg[[10]]<-brca2.deg
  names(ctx.deg)[10]<-"GSE109761"
  ctc.group[[10]]<-annotion
  names(ctc.group)[10]<-"GSE109761"
  
  #11  breast cancer GSE111065
  gse<-read.xlsx(".BR_GSE111065_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:94],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",69),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca3.deg = topTable(fit, n=Inf,sort.by = "P")
  brca3.deg$change<-ifelse(brca3.deg$P.Value<0.05&abs(brca3.deg$logFC)>1,
                           ifelse(brca3.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca3.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[11]]<-gsenor
  names(ctc.expression)[11]<-"GSE111065"
  ctx.deg[[11]]<-brca3.deg
  names(ctx.deg)[11]<-"GSE111065"
  ctc.group[[11]]<-annotion
  names(ctc.group)[11]<-"GSE111065"
  
  #12 breast cancer GSE51827
  gse<-read.xlsx(".BR_GSE51827_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:54],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",29),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca4.deg = topTable(fit, n=Inf,sort.by = "P")
  brca4.deg$change<-ifelse(brca4.deg$P.Value<0.05&abs(brca4.deg$logFC)>1,
                           ifelse(brca4.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca4.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[12]]<-gsenor
  names(ctc.expression)[12]<-"GSE51827"
  ctx.deg[[12]]<-brca4.deg
  names(ctx.deg)[12]<-"GSE51827"
  ctc.group[[12]]<-annotion
  names(ctc.group)[12]<-"GSE51827"
  
  
  #13 breast cancer GSE55807
  gse<-read.xlsx(".BR_GSE55807_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:31],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",6),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca5.deg = topTable(fit, n=Inf,sort.by = "P")
  brca5.deg$change<-ifelse(brca5.deg$P.Value<0.05&abs(brca5.deg$logFC)>1,
                           ifelse(brca5.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca5.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[13]]<-gsenor
  names(ctc.expression)[13]<-"GSE55807"
  ctx.deg[[13]]<-brca5.deg
  names(ctx.deg)[13]<-"GSE55807"
  ctc.group[[13]]<-annotion
  names(ctc.group)[13]<-"GSE55807"

  #14 breast cancer GSE67939
  
  gse<-fread(".GSE67939_readCounts.txt.gz")
  gse<-as.data.frame(gse)
  gene<-fread(".GSE67939_annotation_file.txt.gz")
  gene<-as.data.frame(gene)
  gse<-merge(gse,gene[,c(1,4)],by.x="V1",by.y="ID")
  
  gse<-aggregate(x=gse[,2:18],by=list(gse$symbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  gse<-normalizeBetweenArrays(gse)
  gse<-as.data.frame(gse)
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",15),rep("normal",2))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca6.deg = topTable(fit, n=Inf,sort.by = "P")
  brca6.deg$change<-ifelse(brca6.deg$P.Value<0.05&abs(brca6.deg$logFC)>1,
                           ifelse(brca6.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca6.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[14]]<-gsenor
  names(ctc.expression)[14]<-"GSE67939"
  ctx.deg[[14]]<-brca6.deg
  names(ctx.deg)[14]<-"GSE67939"
  ctc.group[[14]]<-annotion
  names(ctc.group)[14]<-"GSE67939"
  
  #15 breast cancer GSE75367
  gse<-read.xlsx(".BR_GSE75367_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:86],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",61),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca7.deg = topTable(fit, n=Inf,sort.by = "P")
  brca7.deg$change<-ifelse(brca7.deg$P.Value<0.05&abs(brca7.deg$logFC)>1,
                           ifelse(brca7.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca7.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[15]]<-gsenor
  names(ctc.expression)[15]<-"GSE75367"
  ctx.deg[[15]]<-brca7.deg
  names(ctx.deg)[15]<-"GSE75367"
  ctc.group[[15]]<-annotion
  names(ctc.group)[15]<-"GSE75367"
  
  #16  breast cancer GSE86978
  gse<-read.xlsx(".BR_GSE86978_gene_FPKM.xlsx")
  gse<-merge(gse,wbc,by = "EnsemblGene_GeneSymbol")
  gse$EnsemblGene_GeneSymbol<-str_split(gse$EnsemblGene_GeneSymbol,"_",simplify = T)[,2]
  gse<-aggregate(x=gse[,2:102],by=list(gse$EnsemblGene_GeneSymbol),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",77),rep("normal",24))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca8.deg = topTable(fit, n=Inf,sort.by = "P")
  brca8.deg$change<-ifelse(brca8.deg$P.Value<0.05&abs(brca8.deg$logFC)>1,
                           ifelse(brca8.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca8.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[16]]<-gsenor
  names(ctc.expression)[16]<-"GSE86978"
  ctx.deg[[16]]<-brca8.deg
  names(ctx.deg)[16]<-"GSE86978"
  ctc.group[[16]]<-annotion
  names(ctc.group)[16]<-"GSE86978"
  
  #17 breast cancer GSE41245
  gse<-fread(".GSE41245_processed_data.txt.gz")
  gse<-as.data.frame(gse)
  gse<-aggregate(x=gse[,3:32],by=list(gse$ID_REF),FUN = mean)
  row.names(gse)<-gse$Group.1
  gse<-gse[,-1]
  annotion<-data.frame(sample=colnames(gse),group=NA)
  annotion$group<-c(rep("cancer",10),rep("normal",20))
  annotion$group<-factor(annotion$group,levels = c("cancer","normal"))
  design<-model.matrix(~0+annotion$group)
  colnames(design)=levels(annotion$group)
  rownames(design)=annotion$sample
  contrast.matrix<-makeContrasts("cancer-normal",levels=design)
  fit <- lmFit(gse,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit) 
  brca9.deg = topTable(fit, n=Inf,sort.by = "P")
  brca9.deg$change<-ifelse(brca9.deg$P.Value<0.05&abs(brca9.deg$logFC)>1,
                           ifelse(brca9.deg$logFC>1,"UP","DOWN"),"NOT")
  table(brca9.deg$change)
  gsenor<-log2(gse+1)
  ctc.expression[[17]]<-gsenor
  names(ctc.expression)[17]<-"GSE41245"
  ctx.deg[[17]]<-brca9.deg
  names(ctx.deg)[17]<-"GSE41245"
  ctc.group[[17]]<-annotion
  names(ctc.group)[17]<-"GSE41245"
  
  #差异分析用FDR
  for (i in 1:17) {
    ctx.deg[[i]]$change<-ifelse(ctx.deg[[i]]$adj.P.Val<0.05&abs(ctx.deg[[i]]$logFC)>1,
                                ifelse(ctx.deg[[i]]$logFC>1,"UP","DOWN"),"NOT")
  }
  
setwd("../step7/")
save(ctc.expression,file = "CTC17dataset_exorssion.Rdata")
save(ctc.group,file = "CTC17dataset_group.Rdata")

#CTC
load("step7/CTC17dataset_exorssion.Rdata")
load("step7/CTC17dataset_group.Rdata")
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
#ggsave("CTCdata_pathwayscore.pdf",p,width = 7,height = 3)

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
pdf("CTCdata_pathwaytumor_cor.pdf",width = 10,height = 10)
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
pdf("CTCdata_Tcell0tumor_cor.pdf",width = 10,height = 10)
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
pdf("CTCdata_Tcell7tumor_cor.pdf",width = 10,height = 10)
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
pdf("CTCdata_pathwayTcell0_cor.pdf",width = 10,height = 10)
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
pdf("CTCdata_pathwayTcell7_cor.pdf",width = 10,height = 10)
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
pdf("CTCdata_pathwaygene_expression.pdf",width = 12,height = 15)
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
write.csv(pvalue,file = "CTCpvalue.csv")

