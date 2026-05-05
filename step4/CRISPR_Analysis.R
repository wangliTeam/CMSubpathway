library(ggplot2)
library(ggpubr)
library(data.table)

CRISPRGeneEffect <- fread("./step4/data/DepMap/CRISPRGeneEffect.csv", data.table = FALSE)
Model <- fread("./step4/data/DepMap/Model.csv", data.table = FALSE)

# 12 gens
core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")
chen_genelist <- c("X", core_module_genes)
CRISPRGeneEffect_chen <- CRISPRGeneEffect[, chen_genelist, with = FALSE]

#
OncotreeLineage <- Model[Model$OncotreeLineage == "Breast", ]
Model_chen <- OncotreeLineage[OncotreeLineage$OncotreePrimaryDisease != "Non-Cancerous", ]
CRISPRGeneEffect_chen_breast <- CRISPRGeneEffect_chen[CRISPRGeneEffect_chen$X %in% Model_chen$ModelID, ]

mean_gene_effect <- colMeans(CRISPRGeneEffect_chen_breast[, -1, with = FALSE])
mean_gene_effect_df <- data.frame(Gene = names(mean_gene_effect), Mean_Effect = mean_gene_effect)

p_crispr <- ggplot(mean_gene_effect_df, aes(x = Gene, y = Mean_Effect)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  xlab("") + ylab("Mean Gene Effect Score") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("CRISPR Dependency in Breast Cancer Cell Lines")

ggsave("./results/crispr_dependency.pdf", p_crispr, width = 8, height = 4)

saveRDS(list(mean_effect = mean_gene_effect_df, crispr_data = CRISPRGeneEffect_chen_breast),
        "./data/processed/crispr_results.rds")

