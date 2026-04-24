# ==============================================================================
# 02_Differential_Expression_and_GSEA.R
# Description: Differential expression analysis and GSEA-based subpathway 
#              dysregulation assessment with random resampling (100 permutations)
# ==============================================================================

library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Step 1: Load TCGA-BRCA expression data and subpathway definitions ---
load("./data/raw/TCGA/TCGA_exp.rda")
subpathway_gene <- read.csv("./data/processed/subpathway_gene.csv", header = TRUE, row.names = 1)
subpathway_gene <- tidyr::separate_rows(subpathway_gene, gene, sep = ", ")

meta <- unique(subpathway_gene$name)

# --- Step 2: Define tumor/normal groups ---
group <- ifelse(substr(colnames(exp), 14, 15) <= "10" & substr(colnames(exp), 1, 4) == "TCGA", "tumor", "normal")
tumor <- exp[, which(group == "tumor")]
normal <- exp[, which(group == "normal")]

# --- Step 3: 100 random resampling GSEA ---
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

# --- Step 4: Apply screening criteria ---
nes_mat <- nes_conb[, -1]
p_mat <- p_conb[, -1]
adjp_mat <- adjp_conb[, -1]
rownames(nes_mat) <- nes_conb$ID
rownames(p_mat) <- p_conb$ID
rownames(adjp_mat) <- adjp_conb$ID

old_criteria_count <- rowSums(abs(nes_mat) > 1 & p_mat < 0.05, na.rm = TRUE)
new_criteria_count <- rowSums(abs(nes_mat) > 1 & p_mat < 0.05 & adjp_mat < 0.25, na.rm = TRUE)

compare_df <- data.frame(
  Subpathway = rownames(nes_mat),
  Old_Pass_Count = old_criteria_count,
  New_Pass_Count = new_criteria_count
)

old_results <- compare_df[compare_df$Old_Pass_Count >= 80, ]
new_results <- compare_df[compare_df$New_Pass_Count >= 80, ]

# --- Step 5: Save results ---
write.csv(nes_conb, "./data/processed/nes_conb.csv", row.names = FALSE)
write.csv(p_conb, "./data/processed/p_conb.csv", row.names = FALSE)
write.csv(adjp_conb, "./data/processed/adjp_conb.csv", row.names = FALSE)
saveRDS(list(old_results = old_results, new_results = new_results), "./data/processed/gsea_screening_results.rds")

cat("GSEA screening complete. Old criteria:", nrow(old_results), "New criteria:", nrow(new_results), "\n")
