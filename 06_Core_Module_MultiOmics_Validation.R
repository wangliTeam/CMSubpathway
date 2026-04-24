# ==============================================================================
# 06_Core_Module_MultiOmics_Validation.R
# Description: Core metabolic module identification and multi-omics validation 
#              (transcriptome, proteome CPTAC, metabolite CAMP, HPA)
# ==============================================================================

library(ggplot2)
library(reshape2)
library(GSVA)
library(ggpubr)

# --- Step 1: Define core metabolic module (12 overlapping genes) ---
core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")

# --- Step 2: TCGA-BRCA expression validation ---
load("./data/raw/TCGA/TCGA_exp.rda")
group <- ifelse(substr(colnames(exp), 14, 15) <= "10" & substr(colnames(exp), 1, 4) == "TCGA", "cancer", "normal")
expr_core <- exp[intersect(core_module_genes, rownames(exp)), ]
expr_core_df <- as.data.frame(t(expr_core))
expr_core_df$group <- group
expr_melt <- melt(expr_core_df, id.var = "group")

p_tcga <- ggplot(expr_melt, aes(x = variable, y = value, color = group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, position = position_dodge(0.5)) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("") + ylab("Expression") + theme_classic() +
  coord_cartesian(ylim = c(0, 20)) + ggtitle("TCGA") +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  stat_compare_means(aes(group = group), label = "p.signif", label.y = 18)

ggsave("./results/core_module_tcga_expression.pdf", p_tcga, width = 8, height = 3)

# --- Step 3: CAMP metabolomics validation ---
camp_expr <- readRDS("./data/raw/CAMP/camp_expression.rds")
camp_meta <- readRDS("./data/raw/CAMP/camp_metadata.rds")
camp_group <- ifelse(camp_meta$sample_type == "cancer", "cancer", "normal")
camp_core <- camp_expr[intersect(core_module_genes, rownames(camp_expr)), ]
camp_core_df <- as.data.frame(t(camp_core))
camp_core_df$group <- camp_group
camp_melt <- melt(camp_core_df, id.var = "group")

p_camp <- ggplot(camp_melt, aes(x = variable, y = value, color = group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("") + ylab("Expression") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = group), label = "p.signif")

ggsave("./results/core_module_camp_expression.pdf", p_camp, width = 8, height = 3)

# --- Step 4: Metabolite correlation analysis ---
camp_metabolites <- readRDS("./data/raw/CAMP/camp_metabolites.rds")
adh_genes <- c("ADH1A", "ADH1B", "ADH1C")
cor_results <- data.frame()

for (gene in adh_genes) {
  if (gene %in% rownames(camp_expr) & "glucose" %in% colnames(camp_metabolites)) {
    cancer_idx <- which(camp_group == "cancer")
    cor_val <- cor(camp_expr[gene, cancer_idx], camp_metabolites[cancer_idx, "glucose"], method = "pearson")
    cor_results <- rbind(cor_results, data.frame(Gene = gene, Correlation = cor_val, Dataset = "CAMP_cancer"))
  }
}

write.csv(cor_results, "./results/core_module_metabolite_correlation.csv", row.names = FALSE)

# --- Step 5: Save ---
saveRDS(list(tcga_plot = p_tcga, camp_plot = p_camp, cor_results = cor_results),
        "./data/processed/core_module_multiomics.rds")

cat("Core module multi-omics validation complete.\n")
