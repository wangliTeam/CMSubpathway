# ==============================================================================
# 11_Supplementary_Validation.R
# Description: Supplementary validation including Expression Atlas datasets 
#              (E-GEOD-45581, E-MTAB-779) and fat interference experiment (GSE246231)
# ==============================================================================

library(limma)
library(ggplot2)
library(ggpubr)

# --- Step 1: Expression Atlas E-GEOD-45581 ---
# 40 untreated breast cancer + 5 normal control samples
egeod45581 <- read.csv("./data/raw/ExpressionAtlas/E-GEOD-45581/expression_matrix.csv", row.names = 1)
group1 <- c(rep("cancer", 40), rep("normal", 5))
design1 <- model.matrix(~0 + factor(group1))
colnames(design1) <- c("cancer", "normal")

fit1 <- lmFit(egeod45581, design1)
fit1 <- contrasts.fit(fit1, makeContrasts(cancer - normal, levels = design1))
fit1 <- eBayes(fit1)
deg1 <- topTable(fit1, n = Inf, sort.by = "P")

core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")
deg1_core <- deg1[rownames(deg1) %in% core_module_genes, ]

# --- Step 2: Expression Atlas E-MTAB-779 ---
# 20 untreated breast carcinoma + 22 normal control samples
emtab779 <- read.csv("./data/raw/ExpressionAtlas/E-MTAB-779/expression_matrix.csv", row.names = 1)
group2 <- c(rep("cancer", 20), rep("normal", 22))
design2 <- model.matrix(~0 + factor(group2))
colnames(design2) <- c("cancer", "normal")

fit2 <- lmFit(emtab779, design2)
fit2 <- contrasts.fit(fit2, makeContrasts(cancer - normal, levels = design2))
fit2 <- eBayes(fit2)
deg2 <- topTable(fit2, n = Inf, sort.by = "P")
deg2_core <- deg2[rownames(deg2) %in% core_module_genes, ]

# --- Step 3: Combined visualization ---
deg1_core$Dataset <- "E-GEOD-45581"
deg2_core$Dataset <- "E-MTAB-779"
combined_deg <- rbind(deg1_core, deg2_core)

p_atlas <- ggplot(combined_deg, aes(x = Gene, y = logFC, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("") + ylab("Log2 Fold Change") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Expression Atlas Validation")

ggsave("./results/expression_atlas_validation.pdf", p_atlas, width = 10, height = 4)

# --- Step 4: Fat Interference Experiment GSE246231 ---
gse246231 <- readRDS("./data/raw/GSE246231/GSE246231_gene_log2.rds")
fat_group <- gse246231[, grepl("fat", colnames(gse246231), ignore.case = TRUE)]
control_group <- gse246231[, grepl("control", colnames(gse246231), ignore.case = TRUE)]

if (ncol(fat_group) > 0 & ncol(control_group) > 0) {
  group3 <- c(rep("fat", ncol(fat_group)), rep("control", ncol(control_group)))
  expr_combined <- cbind(fat_group, control_group)
  design3 <- model.matrix(~0 + factor(group3))
  colnames(design3) <- c("fat", "control")
  
  fit3 <- lmFit(expr_combined, design3)
  fit3 <- contrasts.fit(fit3, makeContrasts(fat - control, levels = design3))
  fit3 <- eBayes(fit3)
  deg3 <- topTable(fit3, n = Inf, sort.by = "P")
  deg3_core <- deg3[rownames(deg3) %in% core_module_genes, ]
} else {
  deg3_core <- data.frame()
}

# --- Step 5: Save ---
saveRDS(list(
  egeod45581 = deg1_core,
  emtab779 = deg2_core,
  gse246231 = deg3_core
), "./data/processed/supplementary_validation.rds")

write.csv(combined_deg, "./results/expression_atlas_results.csv", row.names = FALSE)

cat("Supplementary validation complete.\n")
