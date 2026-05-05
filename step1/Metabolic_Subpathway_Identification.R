library(igraph)
library(KEGGgraph)
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
#KEGG metabolic pathways ---
meta <- c("00010","00020","00030","00040","00051","00052","00053","00500","00520","00620","00630",
          "00640","00650","00562","00190","00910","00920","00061","00062","00071","00100","00120",
          "00140","00561","00564","00565","00600","00590","00591","00592","01040","00230","00240",
          "00250","00260","00270","00280","00290","00310","00220","00330","00340","00350","00360",
          "00380","00400","00410","00430","00440","00450","00470","00480","00510","00513","00512",
          "00515","00514","00532","00534","00533","00531","00563","00601","00603","00604","00511",
          "00730","00740","00750","00760","00770","00780","00785","00790","00670","00830","00860",
          "00130","00900","00232","00524","00980","00982","00983")

kgml_dir <- "./step1/data/kegglist"
dir.create(kgml_dir, showWarnings = FALSE, recursive = TRUE)

for (pathway_id in meta) {
  kgml_file <- file.path(kgml_dir, paste0("hsa", pathway_id, ".xml"))
  if (!file.exists(kgml_file)) {
    tryCatch({
      keggGet(paste0("hsa", pathway_id), "kgml", file = kgml_file)
    }, error = function(e) {
      message(paste("Failed to download:", pathway_id))
    })
  }
}

setwd(kgml_dir)
idx <- paste0("hsa", meta, ".xml")
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

write.csv(subpathway_d, "./step1/data/subpathway_gene_symbol.csv", row.names = FALSE)
write.csv(subpathway_gene, "./step1/data/subpathway_gene.csv", row.names = FALSE)


load("./step1/data/TCGA_exp.rda")
subpathway_gene <- read.csv("./step1/data/subpathway_gene.csv", header = TRUE, row.names = 1)
subpathway_gene <- tidyr::separate_rows(subpathway_gene, gene, sep = ", ")

meta <- unique(subpathway_gene$name)

#
group <- ifelse(substr(colnames(exp), 14, 15) <= "10" & substr(colnames(exp), 1, 4) == "TCGA", "tumor", "normal")
tumor <- exp[, which(group == "tumor")]
normal <- exp[, which(group == "normal")]

#GSEA
nes_conb <- data.frame(ID = meta, row.names = meta)
p_conb <- data.frame(ID = meta, row.names = meta)
adjp_conb <- data.frame(ID = meta, row.names = meta)

set.seed(123)
for (j in 1:100) {
  tumor_test <- tumor[, sample(ncol(tumor), 0.7 * ncol(tumor))]
  normal_test <- normal[, sample(ncol(normal), 0.7 * ncol(normal))]
  
  exp_test <- cbind(tumor_test, normal_test)
  group_test <- c(rep("tumor", ncol(tumor_test)), rep("normal", ncol(normal_test)))
  
  design <- model.matrix(~0 + factor(group_test))
  colnames(design) <- levels(factor(group_test))
  rownames(design) <- colnames(exp_test)
  
  contrast.matrix <- makeContrasts(tumor - normal, levels = design)
  fit <- lmFit(exp_test, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  DEG <- topTable(fit2, coef = 1, n = Inf, sort.by = "logFC")
  DEG <- na.omit(DEG)
  
  DEG$regulate <- ifelse(DEG$adj.P.Val > 0.05, "unchanged",
                         ifelse(DEG$logFC > 1, "up-regulated",
                                ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
  deg <- DEG[!duplicated(DEG$ID), ]
  rownames(deg) <- deg$ID
  deg <- deg[intersect(rownames(deg), unique(subpathway_gene$gene)), ]
  
  geneList <- deg$logFC
  names(geneList) <- rownames(deg)
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- geneList[geneList != 0]
  
  egmt <- GSEA(geneList, TERM2GENE = subpathway_gene, verbose = FALSE,
               pAdjustMethod = "BH", nPerm = 100, minGSSize = 5, maxGSSize = 1000, pvalueCutoff = 1)
  gsea_results <- egmt@result
  
  nes <- data.frame(ID = gsea_results$ID, number = gsea_results$NES)
  p <- data.frame(ID = gsea_results$ID, number = gsea_results$pvalue)
  adjp <- data.frame(ID = gsea_results$ID, number = gsea_results$p.adjust)
  
  rownames(nes) <- gsea_results$ID
  rownames(p) <- gsea_results$ID
  rownames(adjp) <- gsea_results$ID
  colnames(nes) <- c("ID", j)
  colnames(p) <- c("ID", j)
  colnames(adjp) <- c("ID", j)
  
  nes_conb <- merge(nes_conb, nes, by = "ID", all = TRUE)
  p_conb <- merge(p_conb, p, by = "ID", all = TRUE)
  adjp_conb <- merge(adjp_conb, adjp, by = "ID", all = TRUE)
}

# 

rownames(nes_conb) <- nes_conb$ID
rownames(p_conb) <- p_conb$ID
rownames(adjp_conb) <- adjp_conb$ID
nes_mat <- nes_conb[, -1]
p_mat <- p_conb[, -1]
adjp_mat <- adjp_conb[, -1]

old_criteria_count <- rowSums(abs(nes_mat) > 1 & p_mat < 0.05, na.rm = TRUE)
new_criteria_count <- rowSums(abs(nes_mat) > 1 & p_mat < 0.05 & adjp_mat < 0.25, na.rm = TRUE)

compare_df <- data.frame(
  Subpathway = rownames(nes_mat),
  Old_Pass_Count = old_criteria_count,
  New_Pass_Count = new_criteria_count
)

old_results <- compare_df[compare_df$Old_Pass_Count >= 80, ]
new_results <- compare_df[compare_df$New_Pass_Count >= 80, ]


write.csv(nes_conb, "./step1/data/nes_conb.csv", row.names = FALSE)
write.csv(p_conb, "./step1/data/p_conb.csv", row.names = FALSE)
write.csv(adjp_conb, "./step1/data/adjp_conb.csv", row.names = FALSE)
saveRDS(list(old_results = old_results, new_results = new_results), "./step1/data/gsea_screening_results.rds")

load("./step1/data/meta_mat.rda")
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

library(tidyr)
subpathway_gene <- read.csv("../step1/data/subpathway_gene.csv",header = T,row.names = 1)
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
