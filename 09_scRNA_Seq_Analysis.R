# ==============================================================================
# 09_scRNA_Seq_Analysis.R
# Description: Single-cell RNA-seq analysis for GSE195861 including QC, clustering,
#              malignant cell identification, SiPSiC scoring, CD8 T cell subtypes,
#              and CytoTRACE differentiation analysis
# ==============================================================================

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GSVA)

# --- Step 1: Load scRNA-seq data ---
# Process from raw 10X data
  sample_dirs <- list.dirs("./data/raw/GSE195861", recursive = FALSE)
  object_list <- list()
  
  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    scrna_data <- Read10X(sample_dir)
    seurat_obj <- CreateSeuratObject(counts = scrna_data, min.cells = 3, min.features = 200)
    seurat_obj[["sample"]] <- sample_name
    object_list[[sample_name]] <- seurat_obj
  }
  
  pbmc.harmony <- merge(object_list[[1]], object_list[2:length(object_list)])
  
  # --- Step 2: Quality Control ---
  pbmc.harmony[["percent.mt"]] <- PercentageFeatureSet(pbmc.harmony, pattern = "^MT-")
  pbmc.harmony <- subset(pbmc.harmony, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  
  # --- Step 3: Normalization and PCA ---
  pbmc.harmony <- NormalizeData(pbmc.harmony) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30, verbose = FALSE)
  
  # --- Step 4: Harmony Integration ---
  pbmc.harmony <- RunHarmony(pbmc.harmony, group.by.vars = "orig.ident")
  
  # --- Step 5: Clustering ---
  pbmc.harmony <- JackStraw(pbmc.harmony, num.replicate = 100)
  pbmc.harmony <- ScoreJackStraw(pbmc.harmony, dims = 1:20)
  pbmc.harmony <- FindNeighbors(pbmc.harmony, reduction = "harmony", dims = 1:15) %>%
    FindClusters(resolution = 1.2)
  
  # --- Step 6: Dimensionality Reduction ---
  pbmc.harmony <- RunTSNE(pbmc.harmony, reduction = "harmony", dims = 1:15)
  pbmc.harmony <- RunUMAP(pbmc.harmony, reduction = "harmony", dims = 1:15)
  
  # Save processed object
  save(pbmc.harmony, file = "./data/processed/pbmc.harmony.rda")


# --- Step 7: Cell Type Annotation ---
markers <- list(
  Endothelial = c("PECAM1", "VWF", "STAB1", "RAMP2", "CLDN5", "FLT1"),
  Tc = c("CD3D", "CD3E"),
  B = c("MS4A1", "CD79A", "CD79B"),
  Myeloid = c("CD68", "CD14", "FCGR3A", "LYZ"),
  Fibroblast = c("PDGFRB", "COL1A1", "COL1A2"),
  Epithelial = c("EPCAM", "CDH1", "KRT8", "KRT18"),
  Plasmablasts = c("JCHAIN")
)

cluster_annotation <- list(
  Tc = c(3, 4, 5, 24, 32),
  B = c(2, 6, 7, 16, 17),
  Myeloid = c(10, 12, 25, 30),
  Fibroblast = c(31),
  Epithelial = c(0, 1, 8, 9, 13, 14, 15, 18:23, 26, 28, 29, 33, 34)
)

cluster_ids <- pbmc.harmony@meta.data$seurat_clusters
cluster_ids <- as.numeric(as.character(cluster_ids))

for (cell_type in names(cluster_annotation)) {
  cluster_ids[cluster_ids %in% cluster_annotation[[cell_type]]] <- cell_type
}

pbmc.harmony@meta.data$celltype <- factor(cluster_ids, levels = c("Tc", "B", "Myeloid", "Fibroblast", "Epithelial"))

# --- Step 8: Malignant Cell Identification ---
Idents(pbmc.harmony) <- "celltype"
deg_allcell <- FindAllMarkers(pbmc.harmony, min.pct = 0.25, logfc.threshold = 0.25)
deg_epi <- FindMarkers(pbmc.harmony, ident.1 = "Epithelial", ident.2 = setdiff(levels(pbmc.harmony@meta.data$celltype), "Epithelial"), min.pct = 0.25, logfc.threshold = 0.25)

tumor_genes <- deg_epi$gene[deg_epi$p_val_adj < 0.05 & deg_epi$avg_log2FC > 0.25]

# --- Step 9: SiPSiC Scoring ---
core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")

sc_data <- GetAssayData(pbmc.harmony, layer = "data")
valid_genes <- intersect(core_module_genes, rownames(sc_data))

gsva_scores <- gsva(sc_data[valid_genes, , drop = FALSE], list(CoreModule = valid_genes),
                    method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
pbmc.harmony@meta.data$sipsic_score <- as.numeric(gsva_scores)

# --- Step 10: CD8 T Cell Subtype Analysis ---
tc_cells <- subset(pbmc.harmony, subset = celltype == "Tc")
tc_cells <- NormalizeData(tc_cells) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunHarmony(group.by.vars = "orig.ident")

tc_cells <- FindNeighbors(tc_cells, reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 1.2)
tc_cells <- RunUMAP(tc_cells, reduction = "harmony", dims = 1:15)

# CD8 T cell subtype annotation
cd8_subtypes <- list(
  CD4TN = c(11),
  CD4TREG = c(3, 9, 12, 19),
  CD8TEF = c(1, 6, 15, 16, 18),
  NK = c(13),
  CD8TEX = c(5, 8)
)

subtype_ids <- tc_cells@meta.data$seurat_clusters
subtype_ids <- as.numeric(as.character(subtype_ids))

for (subtype in names(cd8_subtypes)) {
  subtype_ids[subtype_ids %in% cd8_subtypes[[subtype]]] <- subtype
}

tc_cells@meta.data$cd8_subtype <- factor(subtype_ids, levels = names(cd8_subtypes))

# --- Step 11: CytoTRACE Differentiation Analysis ---
if (requireNamespace("CytoTRACE", quietly = TRUE)) {
  counts <- tc_cells@assays$RNA$counts
  cytotrace_results <- CytoTRACE(counts)
  tc_cells@meta.data$cytotrace_score <- cytotrace_results$CytoTraceScore
}

# --- Step 12: Save Results ---
saveRDS(list(
  pbmc = pbmc.harmony,
  tc_cells = tc_cells,
  deg_allcell = deg_allcell,
  deg_epi = deg_epi,
  tumor_genes = tumor_genes
), "./data/processed/scrna_analysis_results.rds")

cat("scRNA-seq analysis complete.\n")
