# ==============================================================================
# 07_Peripheral_Blood_Validation.R
# Description: Core metabolic module validation in 9 breast cancer peripheral 
#              blood cohorts (CTC datasets) with ssGSEA-tumor cell correlation
# ==============================================================================

library(limma)
library(GSVA)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# --- Step 1: Load core metabolic module genes from script 06 ---
core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")

# --- Step 2: Load WBC normal control ---
wbc <- read.xlsx("./data/raw/CTC/wbc_GSE51984_gene_FPKM.xlsx")

# --- Step 3: Define 9 breast cancer peripheral blood datasets ---
brca_ctc_datasets <- list(
  GSE111842 = list(file = "BR_GSE111842_gene_FPKM.xlsx", n_cancer = 13, format = "xlsx", cols = "2:38"),
  GSE109761 = list(file = "BR_GSE109761_gene_FPKM.xlsx", n_cancer = 59, format = "xlsx", cols = "2:84"),
  GSE111065 = list(file = "BR_GSE111065_gene_FPKM.xlsx", n_cancer = 69, format = "xlsx", cols = "2:94"),
  GSE51827  = list(file = "BR_GSE51827_gene_FPKM.xlsx",  n_cancer = 29, format = "xlsx", cols = "2:54"),
  GSE55807  = list(file = "BR_GSE55807_gene_FPKM.xlsx",  n_cancer = 6,  format = "xlsx", cols = "2:31"),
  GSE67939  = list(file = "GSE67939_readCounts.txt.gz",   n_cancer = 15, format = "txt",   cols = NULL),
  GSE75367  = list(file = "BR_GSE75367_gene_FPKM.xlsx",   n_cancer = 61, format = "xlsx", cols = "2:86"),
  GSE86978  = list(file = "BR_GSE86978_gene_FPKM.xlsx",   n_cancer = 77, format = "xlsx", cols = "2:102"),
  GSE41245  = list(file = "GSE41245_processed_data.txt.gz", n_cancer = 8, format = "txt", cols = NULL)
)

# --- Step 4: Unified processing function ---
process_ctc_dataset <- function(gse_name, dataset_info, wbc) {
  raw_dir <- "./data/raw/CTC"
  
  if (dataset_info$format == "xlsx") {
    gse <- read.xlsx(file.path(raw_dir, dataset_info$file))
    gse <- merge(gse, wbc, by.x = "gene", by.y = "EnsemblGene_GeneSymbol")
    gse$gene <- str_split(gse$gene, "_", simplify = TRUE)[, 2]
    col_range <- eval(parse(text = paste0("gse[,", dataset_info$cols, "]")))
    gse <- aggregate(x = col_range, by = list(gse$gene), FUN = mean)
    rownames(gse) <- gse$Group.1
    gse <- gse[, -1]
  } else {
    gse <- fread(file.path(raw_dir, dataset_info$file))
    gse <- as.data.frame(gse)
    if (gse_name == "GSE67939") {
      gene_annot <- fread(file.path(raw_dir, "GSE67939_annotation_file.txt.gz"))
      gene_annot <- as.data.frame(gene_annot)
      gse <- merge(gse, gene_annot[, c(1, 4)], by.x = "V1", by.y = "ID")
      gse <- aggregate(x = gse[, 2:18], by = list(gse$symbol), FUN = mean)
      rownames(gse) <- gse$Group.1
      gse <- gse[, -1]
      gse <- normalizeBetweenArrays(gse)
      gse <- as.data.frame(gse)
    } else if (gse_name == "GSE41245") {
      gse <- aggregate(x = gse[, 3:32], by = list(gse$ID_REF), FUN = mean)
      rownames(gse) <- gse$Group.1
      gse <- gse[, -1]
    }
  }
  
  n_cancer <- dataset_info$n_cancer
  n_normal <- ncol(gse) - n_cancer
  annotation <- data.frame(
    sample = colnames(gse),
    group = factor(c(rep("cancer", n_cancer), rep("normal", n_normal)),
                   levels = c("cancer", "normal"))
  )
  
  design <- model.matrix(~0 + annotation$group)
  colnames(design) <- levels(annotation$group)
  rownames(design) <- annotation$sample
  
  fit <- lmFit(gse, design) %>%
    contrasts.fit(makeContrasts(cancer - normal, levels = design)) %>%
    eBayes()
  
  deg <- topTable(fit, n = Inf, sort.by = "P")
  deg$change <- ifelse(deg$P.Value < 0.05 & abs(deg$logFC) > 1,
                       ifelse(deg$logFC > 1, "UP", "DOWN"), "NOT")
  
  expr_norm <- log2(gse + 1)
  
  list(expr = expr_norm, deg = deg, annotation = annotation, gse_name = gse_name)
}

# --- Step 5: Process all 9 datasets ---
ctc_results <- lapply(names(brca_ctc_datasets), function(name) {
  process_ctc_dataset(name, brca_ctc_datasets[[name]], wbc)
})
names(ctc_results) <- names(brca_ctc_datasets)

# --- Step 6: Extract core module gene expression across datasets ---
expr_summary <- data.frame()
for (res in ctc_results) {
  core_expr <- res$expr[intersect(core_module_genes, rownames(res$expr)), , drop = FALSE]
  if (nrow(core_expr) == 0) next
  
  core_df <- as.data.frame(t(core_expr))
  core_df$Dataset <- res$gse_name
  core_df$Group <- res$annotation$group
  core_df$Sample <- rownames(core_df)
  expr_summary <- rbind(expr_summary, core_df)
}

# --- Step 7: Differential expression summary for core module genes ---
deg_summary <- data.frame()
for (res in ctc_results) {
  core_deg <- res$deg[rownames(res$deg) %in% core_module_genes, ]
  core_deg$Dataset <- res$gse_name
  core_deg$Gene <- rownames(core_deg)
  deg_summary <- rbind(deg_summary, core_deg)
}

# --- Step 8: Visualization - Heatmap of logFC (Fig.S2D style) ---
deg_plot_data <- deg_summary[, c("Gene", "Dataset", "logFC")]
deg_wide <- dcast(deg_plot_data, Gene ~ Dataset, value.var = "logFC")
rownames(deg_wide) <- deg_wide$Gene
deg_mat <- as.matrix(deg_wide[, -1])

p_heatmap <- pheatmap::pheatmap(deg_mat,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("blue", "white", "red"))(50),
                                 display_numbers = TRUE,
                                 number_format = "%.2f",
                                 fontsize = 10,
                                 main = "Core Module Genes in Peripheral Blood Cohorts")
ggsave("./results/peripheral_blood_logfc_heatmap.pdf", p_heatmap, width = 10, height = 4)

# --- Step 9: ssGSEA scoring and tumor cell correlation ---
ssgsea_results <- data.frame()
for (res in ctc_results) {
  gene_set <- list(CoreModule = intersect(core_module_genes, rownames(res$expr)))
  if (length(gene_set$CoreModule) == 0) next
  
  scores <- gsva(as.matrix(res$expr), gene_set, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
  ssgsea_df <- data.frame(
    Sample = colnames(scores),
    ssGSEA_Score = as.numeric(scores),
    Group = res$annotation$group,
    Dataset = res$gse_name
  )
  ssgsea_results <- rbind(ssgsea_results, ssgsea_df)
}

# Correlation with tumor cell abundance (cancer group proportion as proxy)
cor_results <- data.frame()
for (ds in unique(ssgsea_results$Dataset)) {
  ds_data <- ssgsea_results[ssgsea_results$Dataset == ds, ]
  tumor_fraction <- mean(ds_data$Group == "cancer")
  mean_score <- mean(ds_data$ssGSEA_Score[ds_data$Group == "cancer"])
  normal_score <- mean(ds_data$ssGSEA_Score[ds_data$Group == "normal"])
  
  cor_results <- rbind(cor_results, data.frame(
    Dataset = ds,
    Tumor_Fraction = tumor_fraction,
    Mean_Cancer_Score = mean_score,
    Mean_Normal_Score = normal_score,
    Score_Difference = mean_score - normal_score
  ))
}

p_cor <- ggplot(cor_results, aes(x = Tumor_Fraction, y = Score_Difference)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  xlab("Tumor Cell Fraction") + ylab("ssGSEA Score Difference (Cancer - Normal)") +
  theme_classic() +
  ggtitle("Core Module ssGSEA vs Tumor Cell Abundance")

ggsave("./results/peripheral_blood_ssGSEA_correlation.pdf", p_cor, width = 6, height = 4)

# --- Step 10: Save results ---
saveRDS(list(
  ctc_results = ctc_results,
  expr_summary = expr_summary,
  deg_summary = deg_summary,
  ssgsea_results = ssgsea_results,
  cor_results = cor_results
), "./data/processed/peripheral_blood_validation.rds")

write.csv(deg_summary, "./results/peripheral_blood_deg_summary.csv", row.names = FALSE)
write.csv(cor_results, "./results/peripheral_blood_ssGSEA_correlation.csv", row.names = FALSE)

cat("Peripheral blood validation complete. Datasets processed:", length(ctc_results), "\n")
