##step1:metabolic subpathway ------------------------------
meta=c("00010","00020","00030","00040","00051","00052","00053","00500","00520","00620","00630",
       "00640","00650","00562","00190","00910","00920","00061","00062","00071","00100","00120",
       "00140","00561","00564","00565","00600","00590","00591","00592","01040","00230","00240",
       "00250","00260","00270","00280","00290","00310","00220","00330","00340","00350","00360",
       "00380","00400","00410","00430","00440","00450","00470","00480","00510","00513","00512",
       "00515","00514","00532","00534","00533","00531","00563","00601","00603","00604","00511",
       "00730","00740","00750","00760","00770","00780","00785","00790","00670","00830","00860",
       "00130","00900","00232","00524","00980","00982","00983")
library(igraph)
library(KEGGgraph)

mapkG <- parseKGML2Graph("hsa00010.xml",expandGenes=TRUE, genesOnly = TRUE)
setwd("../data")
mapkNodes <- nodes(mapkG)
mapkEdges <- edges(mapkG)
mapkEdges <- mapkEdges[sapply(mapkEdges, length) > 0]
res <- lapply(1:length(mapkEdges), function(t){
  name <- names(mapkEdges)[t]
  len  <- length(mapkEdges[[t]])
  do.call(rbind, lapply(1:len, function(n){
    c(name, mapkEdges[[t]][n])
  }))
})
result <- data.frame(do.call(rbind, res))
write.table(result,  "edges.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(mapkNodes, "nodes.txt", sep = "\t", row.names = F, col.names = F, quote = F)

library(igraph)
#BiocManager::install("KEGGgraph")
library(KEGGgraph)
#idx<-list.files(path = "../data/kegglist",pattern = "*.xml")
idx <- paste0("hsa",meta,".xml")
setwd("../data/kegglist")
subpathway<-list()
for (i in 1:length(idx)) {
  list<-list()
  pathway <- parseKGML2Graph(idx[i],expandGenes=TRUE, genesOnly = TRUE)
  Nodes <- nodes(pathway)
  #if(length(Nodes)<=30){
  #  list<-list(Nodes)
  #  names(list)<-paste(idx[i],1,sep = "-")
  #  subpathway<-append(subpathway,list)
  #  next
  #}
  Edges <- edges(pathway)
  Edges <- Edges[sapply(Edges, length) > 0]
  if(length(Edges)==0){
    list<-list(Nodes)
    names(list)<-paste(idx[i],1,sep = "-")
    subpathway<-append(subpathway,list)
    next
  }
  res <- lapply(1:length(Edges), function(t){
    name <- names(Edges)[t]
    len  <- length(Edges[[t]])
    do.call(rbind, lapply(1:len, function(n){
      c(name, Edges[[t]][n])
    }))
  })
  result <- data.frame(do.call(rbind, res))
  a<-paste(result$X1,result$X2,sep = '-')
  b<-paste(result$X2,result$X1,sep = '-')
  if(length(intersect(a,b))>0){
    result<-result[-which(a%in%intersect(a,b)),]
  }
  net<-graph_from_data_frame(result,directed = F)
  fc<-fastgreedy.community(net)
  if(modularity(fc)<0.3){
    list<-list(Nodes)
    names(list)<-paste(idx[i],1,sep = "-")
    subpathway<-append(subpathway,list)
    next
  }
  for (j in 1:length(fc)) {
    list<-append(list,fc[j]) 
  }
  for (k in sort(1:length(list),decreasing = T)) {
    if(length(list[[k]])<5){list[k]=NULL}
  }
  if(length(list)==0){
    list<-list(Nodes)
    names(list)<-paste(idx[i],1,sep = "-")
    subpathway<-append(subpathway,list)
    next
  }
  for (l in 1:length(list)) {
    names(list)[l]<-paste(idx[i],l,sep = "-")
    
  }
  subpathway<-append(subpathway,list)
}
##NCBI-GeneID <- KEGGID 
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)
subpathway_d <- data.frame()
subpathway_gene <- data.frame()
for (i in 1:length(subpathway)) {
  name <- names(subpathway)[i]
  geneID <- unlist(strsplit(subpathway[[i]],":"))[2*(1:length(subpathway[[i]]))]
  genesymbol <- bitr(geneID,fromType = 'ENTREZID',
                           toType = c('SYMBOL'),
                           OrgDb='org.Hs.eg.db',)
  temp_data <- cbind(data.frame(name=name),genesymbol)
  subpathway_d <- rbind(subpathway_d,temp_data)
  temp_data2 <- data.frame(name=name,gene=paste(genesymbol$SYMBOL,collapse = ", "))
  subpathway_gene <- rbind(subpathway_gene,temp_data2)
  print(i)
}
write.csv(subpathway_d,"../data_new/subpathway_gene_symbol.csv")
write.csv(subpathway_gene,"../data_new/subpathway_gene.csv")

#
##step2:GSEA(Disease-related metabolic subpathways) ------------------------------
#subpathway_gene <- read.csv("../data/subpathway_gene.csv",header = T,row.names = 1)
subpathway_gene <- read.csv("../data_new/subpathway_gene.csv",header = T,row.names = 1)

library(tidyr)

subpathway_gene <- subpathway_gene %>%
  separate_rows(gene, sep = ", ")#dim(subpathway_gene)[1] 6219    2
dim(subpathway_gene)#2803

load("../data/exp.rda")#BRCA表达谱
##group<-anno
dim(exp) #19962  1151
head(colnames(exp))#"TCGA-B6-A0RH-01A-21R-A115-07" "TCGA-BH-A1FU-11A-23R-A14D-07" "TCGA-BH-A1FU-01A-11R-A14D-07"
group<-c(1:ncol(exp))
tumor<-which((substr(colnames(exp),14,15)<=10)&substr(colnames(exp),1,4)=="TCGA")
group[tumor]<-"tumor"
group[-tumor]<-"normal"


tumor<-exp[,which(group=="tumor")]
normal<-exp[,which(group=="normal")]

nes_conb<-data.frame(ID=meta,row.names = meta)

p_conb<-data.frame(ID=meta,row.names = meta)

adjp_conb <-data.frame(ID=meta,row.names = meta)
library(limma)
library(clusterProfiler)
meta <- unique(subpathway_gene$name)
##group<-c(rep("malig",rep(length(malig))),rep("normal",length(normal)))
for (j in 1:100) {
  tumor_test<-tumor[,sample(ncol(tumor),0.7*ncol(tumor))]
  normal_test<-normal[,sample(ncol(normal),0.7*ncol(normal))]
  
  exp_test<-cbind(tumor_test,normal_test)
  group_test<-c(rep("tumor",ncol(tumor_test)),rep("normal",ncol(normal_test)))

  design <- model.matrix(~0+factor(group_test))
  colnames(design) <- levels(factor(group_test))
  rownames(design) <- colnames(exp_test)
  
  ##contrast.matrix <- makeContrasts(tumor-normal,levels = design)
  
  contrast.matrix <- makeContrasts(tumor-normal,levels = design)
  
  fit <- lmFit(exp_test,design) #非线性最小二乘法
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
  DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
  DEG <- na.omit(DEG)
  #colnames(DEG) "ID"        "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"    
   
  DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                         ifelse(DEG$logFC > 1, "up-regulated",
                                ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
  deg<-DEG
  deg<-deg[!duplicated(deg$ID),]
  rownames(deg)<-deg$ID
  
  deg<-deg[intersect(rownames(deg),unique(subpathway_gene$gene)),]
  
  geneList <- deg$logFC
  names(geneList) <- rownames(deg)
  geneList <- sort(geneList, decreasing = T)
  geneList <- geneList[geneList != 0]
  #TERM2GENE是一个必需的data.frame，第一列为term ID，第二列为对应映射基因；
  egmt <- GSEA(geneList, TERM2GENE=subpathway_gene, verbose=FALSE, pAdjustMethod = "BH",
               nPerm = 100, minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)
  gsea_results <- egmt@result
  #colnames(gsea_results)
  # [1] "ID"              "Description"     "setSize"         "enrichmentScore" "NES"             "pvalue"         
  # [7] "p.adjust"        "qvalue"          "rank"            "leading_edge"    "core_enrichment"
  #write.csv(gsea_results,"../data/gsea_results.csv")
  
  nes<-data.frame(ID=gsea_results$ID,number=gsea_results$NES)
  p<-data.frame(ID=gsea_results$ID,number=gsea_results$pvalue)
  adjp<-data.frame(ID=gsea_results$ID,number=gsea_results$p.adjust)#FDR
  
  rownames(nes)<-gsea_results$ID
  rownames(p)<-gsea_results$ID
  rownames(adjp)<-gsea_results$ID
  colnames(nes) <- c("ID",j)
  colnames(p) <- c("ID",j)
  colnames(adjp) <- c("ID",j)
  
  nes_conb<-merge(nes_conb,nes,by="ID",all = T)
  p_conb<-merge(p_conb,p,by="ID",all = T)
  adjp_conb<-merge(adjp_conb,adjp,by="ID",all = T)
  print(j)
}
write.csv(nes_conb,"../data_new/nes_conb.csv")
write.csv(p_conb,"../data_new/p_conb.csv")
write.csv(adjp_conb,"../data_new/adjp_conb.csv")
#abs (mean (NES))>1 and P < 0.05 in 95 or more repetitions
# nes_mat <- read.csv("../data_new/nes_conb.csv", row.names = 1)
# p_mat <- read.csv("../data_new/p_conb.csv", row.names = 1)
# adjp_mat <- read.csv("../data_new/adjp_conb.csv", row.names = 1)
nes_mat <- nes_conb
p_mat <- p_conb
adjp_mat <- adjp_conb

rownames(nes_mat) <- nes_mat$ID
nes_mat <- nes_mat[,-1]
rownames(p_mat) <- p_mat$ID
p_mat <- p_mat[,-1]
rownames(adjp_mat) <- adjp_mat$ID
adjp_mat <- adjp_mat[,-1]
old_criteria_count <- rowSums(abs(nes_mat) > 1 & p_mat < 0.05, na.rm = TRUE)
new_criteria_count <- rowSums(abs(nes_mat) > 1 & p_mat < 0.05 & adjp_mat < 0.25, na.rm = TRUE)# [PMID: 31611988, 2019; PMID37869642, 2023]

compare_df <- data.frame(
  Subpathway = rownames(nes_mat),
  Old_Pass_Count = old_criteria_count,
  New_Pass_Count = new_criteria_count
)

old_results <- compare_df[compare_df$Old_Pass_Count >= 80, ]
dim(old_results) #19 3
print(old_results)

new_results <- compare_df[compare_df$New_Pass_Count >= 80, ]
print(new_results)

##step3:AUC-------------------------
load("meta_gsva_mat.rda")
exp<-t(gsva_mat)
colnames(exp)<-gsub("-","_",colnames(exp))
exp<-exp[,spw]
group<-c(1:nrow(exp))
tumor<-which((substr(rownames(exp),14,15)<=10))
group[tumor]<-1
group[-tumor]<-0
exp<-as.data.frame(exp)
exp$group<-as.factor(group)
data<-exp
for (j in spw) {
  exp<-data[,c(j,"group")]
  colnames(exp)[1]<-"spw"
  x<-exp
  y<-as.factor(exp$group)
  set.seed(123)
  idx<-sample(nrow(x),nrow(x)*0.75)
  x_train<-x[idx,]
  y_train<-y[idx]
  x_test<-x[-idx,]
  y_test<-y[-idx]
  ###key<-c("hsa00071_2" ,"hsa00010_1" ,"hsa00350_1", "hsa00620_2", "hsa00830_1", "hsa00980_1")
  ####"hsa00350_1" "hsa00983_2" "hsa00982_1" "hsa00240_1"
  obj <- tune(svm, group~spw, data = x_train,
              ranges = list(gamma = 2^(-10:3), cost = 2^(-5:4)),
              tunecontrol = tune.control(sampling = "fix")
  )
  model_best<-svm(group~spw,x_train, scale=FALSE, kernel='radial', gamma=obj$best.parameters$gamma,
                  cost= obj$best.parameters$cost,probability = T)
  y_pred <- predict(model_best, x_test,probability = T)
  res<-confusionMatrix(y_test,y_pred,positive = '1')
  y_pred<-as.ordered(y_pred)
  svm.roc = roc(response = x_test$group,
                predictor = y_pred,
                levels = levels(x_test$group))
  auc<-round(svm.roc$auc,4)
  label <- data.frame(
    x=0.9, 
    y=0.01, 
    label = paste("AUC =",auc)
  )
  color<-"#A2AB28"
  plot<-ggroc(svm.roc,color = color,
              linetype = 1,
              size = 1,
              alpha = 1,
              legacy.axes = T)+
    geom_abline(intercept = 0,
                slope = 1,
                color = "grey",
                size = 1,
                linetype = 2)+
    labs(x = "False Postive Rate",
         y = "True Positive Rate")+
    
    coord_cartesian(xlim = c(0,1),
                    ylim = c(0,1))+
    theme_bw()+
    theme(panel.background = element_rect(fill = "transparent"),
    )+#
    geom_label(data = label, aes(x,y,label = label),colour=color)+
    labs(title = j)+theme(plot.title = element_text(size = 20, face = "bold",  hjust = 0.5))
}
    
####pan-cancer###----------------------
#1B
library(tidyverse)
#https://gitlab.com/cvejic-group/lung/-/blob/master/scRNA-Seq/R_pipeline/Plot_clinical.R?ref_type=heads
pan_sub <- read.csv("../data/pancancer_sub.csv", header = T)
# colnames(pan_sub)
# [1] "cancer"     "Subpathway" "NES"        "P"          "AUC"  
pan_sub$Subpathway <- paste0("path:",pan_sub$Subpathway)
pan_sub$NES <- "y" #需要休息
pan_sub$P <- ifelse(pan_sub$P < 0.05 ,"y","n")
pan_sub$AUC <- ifelse(pan_sub$AUC > 0.75,"y","n")

#哪些
pan_sub$sum <- paste0(pan_sub$NES,pan_sub$P,pan_sub$AUC)

#将 NES P AUC扩展到两列
pan_sub_data <- pan_sub[,c("cancer","Subpathway","NES")]
colnames(pan_sub_data) <- c("cancer","Subpathway","value")
pan_sub_data$level <- "NES"
pan_sub_data <- rbind(pan_sub_data,
                      cbind(pan_sub[,c("cancer","Subpathway")],data.frame(value=pan_sub$P,level="P")),
                      cbind(pan_sub[,c("cancer","Subpathway")],data.frame(value=pan_sub$AUC,level="AUC")))

#pan_sub$specific <- ifelse(pan_sub$P < 0.05 & pan_sub$AUC > 0.8,"y","n")


jitter <- position_jitter(width = 0.5, height = 0.5)
set.seed(123)
fig1b<- ggplot(pan_sub_data, aes(x=Subpathway, y=cancer)) +
  geom_jitter(aes(fill=value, shape=level) , size=1, width = 0.25, height = 0.35)+
  scale_shape_manual(values = c(21, 23,24)) +
  #scale_color_manual( c("#af8dc3")) +
  scale_fill_manual(values = c("grey","#af8dc3")) +
  theme_linedraw() +
  theme(panel.grid.major =  element_blank())+
  geom_vline(xintercept = c(1:(length(unique(pan_sub$Subpathway))-1)+0.5), color="gray", linetype = 2)+
  geom_hline(yintercept = c(1:(length(unique(pan_sub$cancer))-1)+0.5), color="gray", linetype = 2)+
  theme(axis.text.x=element_text(angle=45))
pdf("../f1B(AUC0.75).pdf",width = 13,height = 5)
print(fig1b)
dev.off()




#S1A
load("../data/data.rda")

data$B <- gsub("hsa","path:",data$B)
have <- data[data$psum > 95 & abs(data$value) >1,]
library(ggplot2)
figS1a <- ggplot(data, aes(B, variable)) +
  geom_tile(aes(fill = value), colour = "grey", size = 1)+
  scale_fill_gradient2(low = "#5C5DAF",mid = "white",high = "#EA2E2D") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label=p),col ="black",size = 3)
pdf("../figS1a.pdf",width = 12,height = 5)
print(figS1a)
dev.off()

length(intersect(unique(data$B),unique(have$B)))#53
length(unique(data$B))#156
length(unique(have$B))#53
unique(have$variable)

pan_sub <- read.csv("../data/pancancer_sub.csv", header = T)
# colnames(pan_sub)
# [1] "cancer"     "Subpathway" "NES"        "P"          "AUC"  
pan_sub$Subpathway <- paste0("path:",pan_sub$Subpathway)

S1B_data <- data.frame(table(pan_sub$Subpathway))
S1B_data <- S1B_data[order(S1B_data$Freq),]
S1B_data$prop <- S1B_data$Freq/sum(S1B_data$Freq)

figS1b <- ggplot(S1B_data, aes(x = "", y = prop, fill = Freq)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  #geom_text(aes(y = lab.ypos, label = Freq), color = "white")+
  scale_fill_gradient(low = "#FEE0D2", high = "red") +
  theme_void()
pdf("../figS1b.pdf",width = 7,height = 5)
print(figS1b)
dev.off()


#S1C
library(tidyverse)
#https://gitlab.com/cvejic-group/lung/-/blob/master/scRNA-Seq/R_pipeline/Plot_clinical.R?ref_type=heads
pan_sub <- read.csv("../data/pancancer_sub.csv", header = T)
# colnames(pan_sub)
# [1] "cancer"     "Subpathway" "NES"        "P"          "AUC"  
pan_sub$Subpathway <- paste0("path:",pan_sub$Subpathway)
pan_sub$NES <- "y"
pan_sub$P <- ifelse(pan_sub$P < 0.05 ,"y","n")
pan_sub$AUC <- ifelse(pan_sub$AUC > 0.75,"y","n")

table(unlist(pan_sub[1,c("NES","P","AUC")]))
ynum <- c()
for(i in 1:dim(pan_sub)[1]){
  temp_data <- data.frame(table(unlist(pan_sub[i,c("NES","P","AUC")])))
  temp_ynum <- temp_data[temp_data$Var1=="y","Freq"]
  ynum <- c(ynum,temp_ynum)
}
pan_sub$ynum <- ynum
pan_sub_data <- data.frame(table(pan_sub$cancer,pan_sub$ynum))
colnames(pan_sub_data) <- c("cancer","type","number")
pan_sub_data$type <- as.numeric(pan_sub_data$type)
figS1c <- ggplot(pan_sub_data, aes(x=number,y =cancer,fill = type)) +
  geom_bar(stat="identity",position = "stack")+
  theme_classic() +
  geom_text(aes(label = number),size=3)+
  scale_fill_gradient(low = "#FEE0D2", high = "red")
pdf("../画图文件/figS1c.pdf",width = 5,height = 10)
print(figS1c)
dev.off()
##level2
sum(pan_sub_data[pan_sub_data$type!="3","number"])/sum(pan_sub_data$number)
#0.9493671
ds
####BRCA###----------------------
library(tidyr)
subpathway_gene <- read.csv("../data/subpathway_gene.csv",header = T,row.names = 1)
BRCA_sub <- c("hsa00330-4", "hsa00500-1", "hsa00010-1", "hsa00071-2", "hsa00120-2", 
              "hsa00350-1", "hsa00360-2", "hsa00620-2", "hsa00830-1", "hsa00980-1", 
              "hsa00982-1", "hsa00982-2", "hsa00983-1")
BRCA_sub <- subpathway_gene[subpathway_gene$name %in% BRCA_sub,]
BRCA_sub$description <- unlist(strsplit(BRCA_sub$name,"-"))[2*(1:dim(BRCA_sub)[1])-1]

library(KEGGREST)
for (i in 1:dim(BRCA_sub)[1]) {
  
  pathway_info <- keggGet(BRCA_sub[i,"description"])
  description <- pathway_info[[1]][["PATHWAY_MAP"]]
  BRCA_sub[i,"description"] <- description
}
write.csv(BRCA_sub,"../data/BRCA_subways.csv")


# upset plot
BRCA_sub <- read.csv("../data/BRCA_subways.csv",header = T,row.names = 1)

library(UpSetR)
sub_list <- list()
sub_name <- BRCA_sub$name
for(i in 1:length(sub_name)){
  sub_list[[sub_name[i]]] <- unlist(strsplit(BRCA_sub[i,"gene"],", "))
}

library(UpSetR)

pdf("../f1_13sub_upset.pdf",width = 6,height = 5)
p <- upset(fromList(sub_list),nsets = 27,order.by = "freq",
           point.size = 3.5,
           line.size = 2,
           mb.ratio = c(0.45, 0.45),
           text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),
           matrix.color = "#E09392", main.bar.color = "gray75",
           #sets.bar.color = color_sig[1:length(DE_ceRNA_list)], sets.x.label = "Set Size",
           shade.color = "gray88")
print(p)
dev.off()

#gene
hsa_exam <- c("hsa00010-1","hsa00620-2","hsa00071-2")
library(ggVennDiagram)
library(ggplot2)
exa_list <- list()
hsa_exam <- c("hsa00010-1","hsa00620-2","hsa00071-2")
BRCA_sub_exa <- BRCA_sub[BRCA_sub$name %in% hsa_exam,]
for(i in 1:length(hsa_exam)){
  exa_list[[hsa_exam[i]]] <- unlist(strsplit(BRCA_sub_exa[i,"gene"],", "))
}
pdf("../f1_3sub_venn.pdf",width = 6,height = 5)
p <- ggVennDiagram(exa_list,category.names = names(exa_list),
              label = "count",
              label_color = "black",
              label_alpha = 0,
              edge_lty = "dashed",
              edge_size = 1) +
  scale_fill_gradient(low="white",high = "#F0736D",name = "count")
print(p)
dev.off()
intersect(intersect(exa_list[[1]],exa_list[[2]]),exa_list[[3]])
# [1] "ALDH2"   "ALDH1B1" "ALDH9A1" "ALDH3A2" "ALDH7A1" "ADH1A"   "ADH1B"  
# [8] "ADH1C"   "ADH4"    "ADH5"    "ADH6"    "ADH7" 


# GSEA plot
BRCA_subways <- c("hsa00330-4", "hsa00500-1", "hsa00010-1", "hsa00071-2", "hsa00120-2", 
                  "hsa00350-1", "hsa00360-2", "hsa00620-2", "hsa00830-1", "hsa00980-1", 
                  "hsa00982-1", "hsa00982-2", "hsa00983-1")

nes_conb <- read.csv("../data/nes_conb.csv",header = T,row.names = 1)
p_conb <- read.csv("../data/p_conb.csv",header = T,row.names = 1)

p_conb <- p_conb[p_conb$ID %in%  BRCA_subways,]
  
gsea_results <- read.csv("../data/gsea_results.csv",header = T,row.names = 1)
gsea_results <- gsea_results[gsea_results$ID %in% BRCA_subways, ]
BRCA_sub <- read.csv("../data/BRCA_subways.csv",header = T,row.names = 1)
rownames(BRCA_sub) <- BRCA_sub$name
gsea_results$Description <- BRCA_sub[gsea_results$ID,"description"]
write.csv(gsea_results,"../data/BRCA_gsea_results.csv")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
gseaplot2(egmt,BRCA_subways, pvalue_table = TRUE)

##fig1C
library(ggplot2)
library(RColorBrewer)
fig1c_data<- read.csv("../Supplementary Table 1.csv",header = T)
colnames(fig1c_data) <- c("Subpathway","NES","Per_P","Survive_P","AUC" )
subset(fig1c_data,Survive_P<0.05 & AUC >0.8)

fig1c_data$logp <- -log10(fig1c_data$Survive_P)#-log10(0.05)  1.30103
col <- c(brewer.pal(12,"Set3"),brewer.pal(3,"Dark2")[1])

fig1c_left <- ggplot(fig1c_data,aes(NES,Subpathway),color = Subpathway)+#, fill = color
  geom_bar(aes(NES,Subpathway),data=fig1c_data,position=position_dodge(),width=0.05, stat="identity") +
  geom_point(aes(NES,Subpathway, size=Per_P, color = Subpathway),data = fig1c_data)+
  theme_bw() + scale_color_manual(values  = col) + #,limits=c(-2,2)
   geom_vline(xintercept = 1,linetype=2) + #xlim(-2,2) +
  geom_vline(xintercept = -1,linetype=2) + 
  labs(x = "mean(NES) score", y = "subpathways")+ 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
pdf("../f1C_left.pdf",width = 5,height = 6)
print(fig1c_left)
dev.off()

fig1c_mid <- ggplot(fig1c_data,aes(logp,Subpathway),color = Subpathway)+#, fill = color
  #geom_bar(aes(NES,Subpathway),data=fig1c_data,position=position_dodge(),width=0.05, stat="identity") +
  geom_point(aes(logp,Subpathway, size=logp, color = Subpathway),data = fig1c_data)+
  theme_bw() + scale_color_manual(values  = col) + #,limits=c(-2,2)
  geom_vline(xintercept = -log10(0.05),linetype=2) + #xlim(-2,2) +
  labs(x = "-log10(p-value)", y = "subpathways")+ 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
pdf("../f1C_mid.pdf",width = 5,height = 6)
print(fig1c_mid)
dev.off()

fig1c_right <- ggplot(fig1c_data,aes(AUC,Subpathway),color = Subpathway)+#, fill = color
  #geom_bar(aes(NES,Subpathway),data=fig1c_data,position=position_dodge(),width=0.05, stat="identity") +
  geom_point(aes(AUC,Subpathway, size=AUC, color = Subpathway),data = fig1c_data)+
  theme_bw() + scale_color_manual(values  = col) + #,limits=c(-2,2)
  geom_vline(xintercept = 0.8,linetype=2) + #xlim(-2,2) +
  labs(x = "AUC", y = "subpathways")+ 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
pdf("../f1C_right.pdf",width = 5,height = 6)
print(fig1c_right)
dev.off()

BRCA_sub <- read.csv("../data/BRCA_subways.csv",header = T,row.names = 1)
##
BRCA_sub$geneNum <- lengths(strsplit(BRCA_sub$gene,", "))
fig1c_right2 <- ggplot(BRCA_sub,aes(geneNum,name, fill=name)) + 
  geom_bar(stat="identity")+
  geom_text(aes(label =geneNum))+ theme_bw() + scale_fill_manual(values  = col) + #,limits=c(-2,2)
  labs(x = "gene numbers", y = "subpathways")+ 
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
pdf("../f1C_right_num.pdf",width = 5,height = 6)
print(fig1c_right2)
dev.off()
