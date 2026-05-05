library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(GSVA)
library(stringr)
library(Seurat)
library(data.table)
library(ggpubr)
library(ggplot2)
raw_dir <- "./step6/data/GSE195861_RAW"#download the data from GEO 
all_files <- list.files(raw_dir, pattern = "\\.gz$", full.names = TRUE)
sample_names <- unique(gsub("^[^_]+_([^_]+)_.*", "\\1", basename(all_files)))
print(sample_names)
# [1] "NCCBC2"  "NCCBC3"  "NCCBC5"  "NCCBC6"  "NCCBC7"  "NCCBC11" "NCCBC13" "NCCBC14" "P2"      "M2"      "P1"      "M1"     
# [13] "P6"      "M6"      "P5"      "M5"      "P4"      "M4"      "P3"      "M3"    
for (s in sample_names) {
  sample_dir <- file.path(raw_dir, s)
  if (!dir.exists(sample_dir)) dir.create(sample_dir)
  barcode_file <- grep(sprintf(".*%s.*barcodes\\.tsv\\.gz$", s), all_files, value = TRUE)
  feature_file  <- grep(sprintf(".*%s.*features\\.tsv\\.gz$", s), all_files, value = TRUE)
  matrix_file   <- grep(sprintf(".*%s.*matrix\\.mtx\\.gz$", s), all_files, value = TRUE)
  if (length(barcode_file) == 1) {
    file.copy(barcode_file, file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  }
  if (length(feature_file) == 1) {
    file.copy(feature_file, file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
  }
  if (length(matrix_file) == 1) {
    file.copy(matrix_file, file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)
  }
}
sample_dirs <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
Object_list <- list()
for (dir in sample_dirs) {
  sample_name <- basename(dir)
  scrna_data <- Read10X(dir)
  seurat_obj <- CreateSeuratObject(scrna_data, min.cells = 3, min.features = 200)
  seurat_obj[["sample"]] <- sample_name
  Object_list[[sample_name]] <- seurat_obj
}

pbmc<-merge(Object_list[[1]],Object_list[2:length(Object_list)])
library(harmony)
library(dplyr)
library(Seurat)
library(data.table)
library(patchwork)
library(ggplot2)
library(gridExtra)
setwd("./step6/data")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
pbmc <- NormalizeData(pbmc) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
pbmc.harmony<-RunHarmony(pbmc,group.by.vars = "orig.ident" )
pbmc.harmony<- JackStraw(pbmc.harmony, num.replicate = 100)
pbmc.harmony<- ScoreJackStraw(pbmc.harmony, dims = 1:20)
ElbowPlot(pbmc.harmony)
pbmc.harmony <- FindNeighbors(pbmc.harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1.2)
pbmc.harmony <- RunTSNE(pbmc.harmony, reduction = "harmony", dims = 1:15)
pbmc.harmony <- RunUMAP(pbmc.harmony, reduction = "harmony", dims = 1:15)
new.cluster.ids <- c("Endothelial","Tc","B","Myeloid","Fibroblast","Epithelial","Plasmablasts")
markers<-c(Endothelial=c("PECAM1","VWF","STAB1","RAMP2","CLDN5","FLT1"),
           Tc=c("CD3D","CD3E"),
           B=c("MS4A1","CD79A","CD79B"),
           Myeloid=c("CD68","CD14","FCGR3A","LYZ"),
           Fibroblast=c("PDGFRB","COL1A1","COL1A2"),
           Epithelial=c("EPCAM","CDH1","KRT8","KRT18"),
           Plasmablasts=c("JCHAIN"))


marker<-c(
  Tc=c("CD3D","CD3E"),
  B=c("MS4A1","CD79A","CD79B"),
  Myeloid=c("CD68","CD14","FCGR3A","LYZ"),
  Fibroblast=c("PDGFRB","COL1A1","COL1A2"),
  Epithelial=c("EPCAM","CDH1","KRT8","KRT18"))
names(markers)<-new.cluster.ids

VlnPlot(pbmc.harmony, features = markers,group.by = "seurat_clusters",assay = "RNA",
        
        stack=T,pt.size=0,
        flip = T,
        add.noise = T)+#
  theme(axis.text.y = element_blank(), #
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')

Endothelial<-c()
Tc<-c(3,4,5,24,32)
B<-c(2,6,7,16,17)
Myeloid<-c(10,12,25,30)
Fibroblast<-c(31)
Epithelial<-c(0,1,8,9,13,14,15,18:23,26,28,29,33,34)
Plasmablasts<-c()

a<-c(Tc,B,Myeloid,Fibroblast,Epithelial)
b<-pbmc.harmony@meta.data$seurat_clusters
b<-as.vector(b)
b<-as.numeric(b)
list<-list(
  Tc,B,Myeloid,Fibroblast,Epithelial
)
names(list)<-c("Tc","B","Myeloid","Fibroblast","Epithelial")

for (i in 1:length(list)) {
  b[b%in%list[[i]]]<-names(list)[i]
}
b<-as.factor(b)
pbmc.harmony@meta.data$celltype<-b

DimPlot(pbmc.harmony, reduction = "umap",pt.size=0.5, label = F,repel = TRUE)
save(pbmc.harmony,file = "./step6/data/pbmc.harmony.rda")
#####Tc#####
load("./step6/data/pbmc.harmony.rda")
meta<-pbmc.harmony@meta.data
counts<-pbmc.harmony@assays$RNA$counts
cell<-colnames(pbmc.harmony)[which(pbmc.harmony$celltype=="Tc")]
counts<-counts[,cell]
pbmc<-CreateSeuratObject(counts = counts,project = "pbmc3k")
meta<-meta[colnames(counts),]
pbmc$orig.ident<-meta$orig.ident
pbmc<-JoinLayers(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes)
rb.genes <- rownames(pbmc)[grep("^RP[SL]",rownames(pbmc))]
percent.ribo <- colSums(counts[rb.genes,])/Matrix::colSums(counts)*100
pbmc <- AddMetaData(pbmc, percent.ribo, col.name = "percent.ribo")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 3)
pbmc <- subset(pbmc, subset =  percent.mt < 10)
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures =2000, verbose = FALSE)
pbmc<-ScaleData(pbmc)
pbmc<-RunPCA(pbmc,verbose = T)
pbmc.harmony<-RunHarmony(pbmc,group.by.vars = "orig.ident" )
pbmc.harmony<- JackStraw(pbmc.harmony, num.replicate = 100)
pbmc.harmony<- ScoreJackStraw(pbmc.harmony, dims = 1:20)
ElbowPlot(pbmc.harmony)
pbmc.harmony <- FindNeighbors(pbmc.harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1.2)
pbmc.harmony <- RunTSNE(pbmc.harmony, reduction = "harmony", dims = 1:15)
pbmc.harmony <- RunUMAP(pbmc.harmony, reduction = "harmony", dims = 1:15)

DimPlot(pbmc.harmony, reduction = "umap",pt.size=0.5, label = F,repel = TRUE)

deg<-FindAllMarkers(pbmc.harmony,min.pct = 0.5, logfc.threshold = 0.25)

deg %>%
  group_by(cluster) %>% #
  top_n(n = 15, wt = avg_log2FC) %>% #
  dplyr::arrange(desc(avg_log2FC))%>%
  dplyr::arrange(cluster)-> top15 #
markers<-c("SELL","CCR7","TCF7","FOXP3","CXCL13","IL21","PDCD1","IL7R","FGFBP2",
           "LAG3","TIGIT","IFNG","KLRC1","KLRF1",
           "KLRB1","CD3D","CD4","CD8A","CD8B")
markers<-c(   CD=c("CD4","CD8A","CD8B","CD3D","CD3E"),
              memory=c("CD27","CD28","IL7R"),
              tfh=c("TOX2","CXCR5","ICOS"),
              naive=c("SELL","CCR7","LEF1"),
              regulatory=c("FOXP3","TNFRSF4","TIGIT","CTLA4"),
              exhausted=c("PDCD1","LAG3","HAVCR2"),
              cytotoxic=c("FGFBP2","GNLY","GZMH","GZMB","GZMK","GZMM","GZMA"),
              nk=c("NKG7","NCAM1")
)           
markers<-c(   TH17=c("IL17A","RORC"),
              memory=c("IL7R")
) 
VlnPlot(pbmc.harmony, features = markers,assay = "RNA",group.by = "celltype",
        ##group.by = "celltype",
        
        stack=T,pt.size=0,
        flip = T,
        add.noise = T)+#
  theme(axis.text.y = element_blank(), #
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')

gene<-c("CCR7","SELL","LEF1","IL7R","CD52","S100A11","KLRF1","FCER1G","KLRD1",
        "FGFBP2","GNLY","GZMH","GZMB","GZMA","GZMK","NKG7",
        "FCGR3A",
        "CD52","CD160","GPR65","PDCD1","CTLA4","LAG3","HAVCR2",
        "IL2RB","TRGC1",
        "TRDC","TRGC2","FOXP3","CTLA4","TNFRSF4","TIGIT","CD40LG",
        "IL1R2")
CD4TN=c(11)
CD4TREG=c(3,9,12,19)
CD4TEM=c(10)
CD4TCM=c(0,2,17)
CD4TFH=c(4)
CD8TEF=c(1,6,15,16,18)
NK=c(13)
CD8TEX=c(5,8)
c0_CD4TN_GPR183=c(0)
c1_CD8TEF_GZMK=c(1)
c2_γδT_TRGC1=c(2)
c3_NK_B3GNT7=c(3)
c4_CD8TEF_FGFBP2=c(4)
c5_MAIT_SLC7A5=c(5)
c6_CD8T_HSPA1B=c(6)
c7_CD4TREG_CTLA4=c(7)
c8_CD8TEF_FGFBP2=c(8)
c9_CD8Tex_DUSP4=c(9)
c10_NK_CCL3=c(10)
b<-pbmc.harmony@meta.data$seurat_clusters
b<-as.vector(b)
b<-as.numeric(b)

list<-list(
  CD4TN=c(11),
  CD4TREG=c(3,9,12,19),
  CD4TEM=c(10),
  CD4TCM=c(0,2,17),
  CD4TFH=c(4),
  CD8TEF=c(1,6,15,16,18),
  NK=c(13),
  CD8TEX=c(5,8)
  
)

for (i in 1:length(list)) {
  b[b%in%list[[i]]]<-names(list)[i]
}
b<-factor(b,levels = names(list))
pbmc.harmony@meta.data$celltype<-b

b<-pbmc.harmony@meta.data$seurat_clusters
b<-as.vector(b)
b<-as.numeric(b)

list<-list(
  CD4T=c(0,7),
  CD8Tactivate=c(1,2,4,5,6,8,10),
  NK=c(3,10), 
  CD8Texhausted=c(9))

for (i in 1:length(list)) {
  b[b%in%list[[i]]]<-names(list)[i]
}
b<-factor(b,levels = names(list))
pbmc.harmony@meta.data$celltype<-b
####B####
load("./step6/data/pbmc.harmony.rda")
meta<-pbmc.harmony@meta.data
counts<-pbmc.harmony@assays$RNA$counts
cell<-colnames(pbmc.harmony)[which(pbmc.harmony$celltype=="B")]
counts<-counts[,cell]
pbmc<-CreateSeuratObject(counts = counts,project = "pbmc3k")
pbmc<-JoinLayers(pbmc)
meta<-meta[cell,]
pbmc$orig.ident<-meta$orig.ident
Idents(pbmc)<-pbmc$orig.ident
pbmc <- NormalizeData(pbmc) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
pbmc.harmony<-RunHarmony(pbmc,group.by.vars = "orig.ident" )
pbmc.harmony<- JackStraw(pbmc.harmony, num.replicate = 100)
pbmc.harmony<- ScoreJackStraw(pbmc.harmony, dims = 1:20)
ElbowPlot(pbmc.harmony)
pbmc.harmony <- FindNeighbors(pbmc.harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1.2)
pbmc.harmony <- RunTSNE(pbmc.harmony, reduction = "harmony", dims = 1:15)
pbmc.harmony <- RunUMAP(pbmc.harmony, reduction = "harmony", dims = 1:15)

DimPlot(pbmc.harmony, reduction = "umap",pt.size=0.5, label = F,repel = TRUE)

deg<-FindAllMarkers(pbmc.harmony,min.pct = 0.5, logfc.threshold = 0.25)

deg %>%
  group_by(cluster) %>% #
  top_n(n = 15, wt = avg_log2FC) %>% #
  dplyr::arrange(desc(avg_log2FC))%>%
  dplyr::arrange(cluster)-> top15 #
markers<-c(Mature=c("MS4A1","CD19"),
           BFollicular=c("NR4A2","FCER2","CD40","CD86","CD22","CR2"),
           Bm=c("CD83","CD69","AIM2", "TNFRSF13B", "CD27" ),
           Bn=c("TCL1A", "IGHD" , "IL4R","CD22"),
           Plasma=c("MZB1","IGHG1", "CD38" , "JCHAIN","SDC1"),
           regulator=c("GZMB"))
markers<-c("MS4A1","CD19","CD27","MZB1","IGHG1","JCHAIN")


VlnPlot(pbmc.harmony, features = markers,assay = "RNA",group.by = "celltype",
        
        stack=T,pt.size=0,
        flip = T,
        add.noise = T)+#
  theme(axis.text.y = element_blank(), #
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')
Bn=c(0,4,13,14,16,19,20)
Bm=c(2,3,5,6,7,8,9,10,11,12,15,17,18)
Plasma=c(1)
b<-pbmc.harmony@meta.data$seurat_clusters
b<-as.vector(b)
b<-as.numeric(b)
list<-list(
  Bn=c(0,4,13,14,16,19,20),
  Bm=c(2,3,5,6,7,8,9,10,11,12,15,17,18),
  Plasma=c(1)
)

for (i in 1:length(list)) {
  b[b%in%list[[i]]]<-names(list)[i]
}
b<-factor(b,levels = names(list))
pbmc.harmony@meta.data$celltype<-b

####myeloid####
load("./step6/data/pbmc.harmony.rda")
meta<-pbmc.harmony@meta.data
counts<-pbmc.harmony@assays$RNA$counts
cell<-colnames(pbmc.harmony)[which(pbmc.harmony$celltype=="Myeloid")]
counts<-counts[,cell]
pbmc<-CreateSeuratObject(counts = counts,project = "pbmc3k")
pbmc<-JoinLayers(pbmc)
meta<-meta[cell,]
pbmc$orig.ident<-meta$orig.ident
Idents(pbmc)<-pbmc$orig.ident
pbmc <- NormalizeData(pbmc) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
pbmc.harmony<-RunHarmony(pbmc,group.by.vars = "orig.ident" )
pbmc.harmony<- JackStraw(pbmc.harmony, num.replicate = 100)
pbmc.harmony<- ScoreJackStraw(pbmc.harmony, dims = 1:20)
ElbowPlot(pbmc.harmony)
pbmc.harmony <- FindNeighbors(pbmc.harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1.2)
pbmc.harmony <- RunTSNE(pbmc.harmony, reduction = "harmony", dims = 1:15)
pbmc.harmony <- RunUMAP(pbmc.harmony, reduction = "harmony", dims = 1:15)
DimPlot(pbmc.harmony, reduction = "umap",pt.size=0.5, label = F,repel = TRUE)
deg<-FindAllMarkers(pbmc.harmony,min.pct = 0.25, logfc.threshold = 0.25)

deg %>%
  group_by(cluster) %>% #
  top_n(n = 15, wt = avg_log2FC) %>% #
  dplyr::arrange(desc(avg_log2FC))%>%
  dplyr::arrange(cluster)-> top15 #
markers<-c(PDC=c("LILRA4","SLC15A4","PLD4","CCDC50","IL3RA","LY9","SELL","GAS6"),  
           CDC1=c("CADM1","XCR1"),
           CDC2=c("CD1C","CD1E","FCER1A","CLEC10A"),
           CDC1_CCR7_LAMP3=c("CCR7","LAMP3","FSCN1"),
           Mast=c("KIT","CD9"),
           Neutrophils=c("HCAR3","HCAR2","CXCL8","OSM","IFITM2"),
           mono=c("S100A12","CD300E","VCAN","FCN1","EREG"),
           Macrophage=c("C1QA","C1QB","C1QC","FBP1","SPP1","FABP5","APOC1"))
markers<-c(mono=c("S100A12","CD300E","VCAN","FCN1"),
           DC=c("CD1C","CD1E","FCER1A","CLEC10A"),
           Neutrophils=c("FCGR3B"),
           Macro=c("C1QA","C1QB","C1QC","CD68","FCGR3A","CD14","MERTK")
)
VlnPlot(pbmc.harmony, features = gene,assay = "RNA",
        
        stack=T,pt.size=0,
        flip = T,
        add.noise = T)+#
  theme(axis.text.y = element_blank(), #
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')

DC=c(2)
Mono=c(5)
Neutrophils=c(3)
Macrophage=c(0,1,4,6:19)

b<-pbmc.harmony@meta.data$seurat_clusters
b<-as.vector(b)
b<-as.numeric(b)

list<-list(
  DC=c(2),
  Mono=c(5),
  Neutrophils=c(3),
  Macrophage=c(0,1,4,6:19)
  
)
for (i in 1:length(list)) {
  b[b%in%list[[i]]]<-names(list)[i]
}
b<-factor(b,levels = names(list))
pbmc.harmony@meta.data$celltype<-b
meta<-pbmc.harmony@meta.data
counts<-pbmc.harmony@assays$RNA$counts
cell<-colnames(pbmc.harmony)[which(pbmc.harmony$celltype=="Fibroblast")]
counts<-counts[,cell]
meta<-meta[cell,]
assayData <- SingleCellExperiment(assays = list(counts = counts))
scoresAndIndices <- getPathwayScores(counts(assayData), key1) # The third parameter, percentForNormalization, is optional; If not specified, its value is set to 5.
pathwayScoresOfCells <- scoresAndIndices$pathwayScores
pathwayGeneIndices <- scoresAndIndices$index
meta$sipsic<-pathwayScoresOfCells
data<-meta[,c("group","sipsic")]
ggplot(data,aes(group,sipsic,fill=group))+
  geom_violin(scale = "width",alpha=0.8,width=0.5,size=0.8)+
  scale_fill_manual(values = c("#F7903D","4D85BD"))+ #
  stat_compare_means(                     
    method = "t.test",
    paired = F,                             
    label.y = 1,
    size=4.5)+
  xlab("")+                                                   #
  ylab("Fraction")+                                           #
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),                  #
        axis.text.x = element_text(angle=0,size=10,vjust = 1,hjust =1,color = "black"),
        axis.text.y = element_text(size =10),
        legend.position = c(0.9,0.85) )+scale_y_continuous(limits = c(0,1))


ggplot(data,aes(group,sipsic,fill=group))+
  geom_boxplot()+
  scale_fill_manual(values = c("#F7903D","4D85BD"))+ #
  stat_compare_means(                     
    method = "t.test",
    paired = F,                             
    label.y = 1,
    size=4.5)+
  xlab("")+                                                   #
  ylab("Fraction")+                                           #
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),                  #
        axis.text.x = element_text(angle=0,size=10,vjust = 1,hjust =1,color = "black"),
        axis.text.y = element_text(size =10),
        legend.position = c(0.9,0.85) )+scale_y_continuous(limits = c(0,1))      

for (i in 1:length(Metspwlist)) {
  gene<-Metspwlist[[i]]
  id<-mapIds(
    x=org.Hs.eg.db,
    keys = gene,
    keytype = "ENTREZID",
    column = "SYMBOL"
  )
  Metspwlist[[i]]<-id
}

res<-data.frame()
for (i in levels(pbmc.harmony$seurat_clusters)) {
  
  a<-data[which(data$group%in%c(i,"Tc","B","Myeloid")),]
  a$cluster<-ifelse(a$cluster=="Epithelial","Epithelial","immue")
  
  b<- wilcox_test(a,score ~ cluster, detailed = T)
  res<-rbind(res,b)
}

fc<-c()
for (i in levels(pbmc.harmony$seurat_clusters)) {
  a<-data[which(data$group%in%c(i,"Tc","B","Myeloid")),]
  a$cluster<-ifelse(a$cluster=="Epithelial","Epithelial","immue")
  mean1<-mean(a$score[which(a$cluster=="Epithelial")])
  mean2<-mean(a$score[which(a$cluster=="immue")])
  b<-log2(mean1/mean2)
  fc<-c(fc,b)
}
data<-meta[,c("sipsic","group")]
norcol="#EE0000FF"                    #disease group
concol="#3B4992FF"  
ggviolin(data, x="group", y="sipsic", color = "group", 
         
         
         #legend.title=x,
         add.params = list(fill="white"),
         palette = c(concol,norcol),
         width=1, add = "boxplot")+
  stat_compare_means(aes(group=group),
                     method="wilcox.test" #
  )+scale_y_continuous(limits = c(0, 1))

ggviolin(data, x="Gene", y="Expression", color = "Type", 
         ylab="Gene expression",
         xlab="Gene",
         #legend.title=x,
         add.params = list(fill="white"),
         palette = c(concol,norcol),
         width=1, add = "boxplot")