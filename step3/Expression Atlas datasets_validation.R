library(limma)
dataset1_exp <- as.data.frame(data.table::fread('./step3/data/E-GEOD-45581-A-AGIL-28-normalized-expressions.tsv', header = TRUE))
dataset1_group <- as.data.frame(data.table::fread('./step3/data/E-GEOD-45581-experiment-design.tsv', header = TRUE))  
expr_data <- dataset1_exp[, -(1:3)]
expr_data$Gene_Symbol <- dataset1_exp$`Gene Name`
expr_agg <- aggregate(. ~ Gene_Symbol, data = expr_data, FUN = mean)
rownames(expr_agg) <- expr_agg$Gene_Symbol
expr_matrix <- as.matrix(expr_agg[, -1])
if (any(is.na(expr_matrix))) {
  expr_matrix[is.na(expr_matrix)] <- 0
}
group_info <- dataset1_group[, c('Assay', 'Factor Value[disease]')]
colnames(group_info) <- c('Assay', 'Disease')
target_groups <- c("inflammatory breast cancer", "normal")
group_filtered <- group_info[group_info$Disease %in% target_groups, ]
group_filtered$Assay_clean <- gsub("_\\d+$", "", group_filtered$Assay)  
matched_samples <- colnames(expr_matrix) %in% group_filtered$Assay
expr_matched <- expr_matrix[, matched_samples]
sample_disease <- character(ncol(expr_matched))
for (i in 1:ncol(expr_matched)) {
  sample_name <- colnames(expr_matched)[i]
  match_idx <- grep(sample_name, group_filtered$Assay)
  if (length(match_idx) > 0) {
    sample_disease[i] <- group_filtered$Disease[match_idx[1]]
  }
}
valid_idx <- sample_disease != ""
expr_final <- expr_matched[, valid_idx]
sample_disease <- sample_disease[valid_idx]
library(limma)
design <- model.matrix(~ 0 + factor(sample_disease))
colnames(design) <- c("IBC", "Normal")
fit <- lmFit(expr_final, design)
contrast_matrix <- makeContrasts(IBC_vs_Normal = IBC - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results <- topTable(fit2, coef = "IBC_vs_Normal", number = Inf, adjust.method = "BH")
results$Gene_Symbol <- rownames(results)
results$Significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, "Yes", "No")
results$Is_Target_Gene <- results$Gene_Symbol %in% gene_set
results$negLog10P <- -log10(results$adj.P.Val)
            sum(results$adj.P.Val < 0.05 & results$logFC > 1)))
            sum(results$adj.P.Val < 0.05 & results$logFC < -1)))
            sum(results$Is_Target_Gene & results$adj.P.Val < 0.05)))
write.csv(results, "limma_DEA_IBC_vs_Normal.csv", row.names = FALSE)
results <- results %>%
  mutate(
    Direction = case_when(
      logFC >= 1 & Significant == "Yes" ~ "Up",
      logFC <= -1 & Significant == "Yes" ~ "Down",
      TRUE ~ "Stable"
    )
  )
my_colors <- c("Up" = "#BC3C29", "Down" = "#0072B5", "Stable" = "#E6E6E6")
volcano_plot <- ggplot(results, aes(x = logFC, y = negLog10P)) +
  geom_point(data = subset(results, Direction == "Stable"),
             color = my_colors["Stable"], alpha = 0.4, size = 0.8) +
  geom_point(data = subset(results, Direction != "Stable" & Is_Target_Gene == FALSE),
             aes(color = Direction), alpha = 0.3, size = 1.2) +
  geom_point(data = subset(results, Is_Target_Gene == TRUE),
             aes(fill = Direction), shape = 21, color = "black", 
             size = 2.8, stroke = 0.6, alpha = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_text_repel(
    data = subset(results, Is_Target_Gene == TRUE),
    aes(label = Gene_Symbol),
    size = 3.5,
    fontface = "bold.italic",
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey30",
    segment.linewidth = 0.4,
    force = 8
  ) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1, color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "E-GEOD-45581",
    x = expression(bold(log[2]("Fold Change"))),
    y = expression(bold(-log[10]("FDR")))
  ) +
  scale_x_continuous(expand = expansion(mult = 0.1))
ggsave("volcano_plot_IBC_vs_Normal.pdf", plot = volcano_plot, width = 10, height = 8)


dataset2_exp <- as.data.frame(data.table::fread('./step3/data/E-MTAB-779-A-AGIL-28-normalized-expressions.tsv', header = TRUE))
dataset2_group <- as.data.frame(data.table::fread('./step3/data/E-MTAB-779-experiment-design.tsv', header = TRUE))  
expr_data2 <- dataset2_exp[, -(1:3)]
expr_data2$Gene_Symbol <- dataset2_exp$`Gene Name`
expr_agg2 <- aggregate(. ~ Gene_Symbol, data = expr_data2, FUN = mean)
rownames(expr_agg2) <- expr_agg2$Gene_Symbol
expr_matrix2 <- as.matrix(expr_agg2[, -1])
if (any(is.na(expr_matrix2))) {
  expr_matrix2[is.na(expr_matrix2)] <- 0
}
group_info2 <- dataset2_group[, c('Assay', 'Factor Value[disease]')]
colnames(group_info2) <- c('Assay', 'Disease')
target_groups2 <- c("breast carcinoma", "normal")
group_filtered2 <- group_info2[group_info2$Disease %in% target_groups2, ]
matched_samples2 <- colnames(expr_matrix2) %in% group_filtered2$Assay
expr_matched2 <- expr_matrix2[, matched_samples2]
sample_disease2 <- character(ncol(expr_matched2))
for (i in 1:ncol(expr_matched2)) {
  sample_name <- colnames(expr_matched2)[i]
  match_idx <- grep(sample_name, group_filtered2$Assay, fixed = TRUE)
  if (length(match_idx) > 0) {
    sample_disease2[i] <- group_filtered2$Disease[match_idx[1]]
  }
}
valid_idx2 <- sample_disease2 != ""
expr_final2 <- expr_matched2[, valid_idx2]
sample_disease2 <- sample_disease2[valid_idx2]
design2 <- model.matrix(~ 0 + factor(sample_disease2))
colnames(design2) <- c("BC", "Normal")
fit_bc <- lmFit(expr_final2, design2)
contrast_matrix2 <- makeContrasts(BC_vs_Normal = BC - Normal, levels = design2)
fit_bc2 <- contrasts.fit(fit_bc, contrast_matrix2)
fit_bc2 <- eBayes(fit_bc2)
results2 <- topTable(fit_bc2, coef = "BC_vs_Normal", number = Inf, adjust.method = "BH")
results2$Gene_Symbol <- rownames(results2)
results2$Significant <- ifelse(results2$adj.P.Val < 0.05 & abs(results2$logFC) > 1, "Yes", "No")
results2$Is_Target_Gene <- results2$Gene_Symbol %in% gene_set
results2$negLog10P <- -log10(results2$adj.P.Val)
            sum(results2$adj.P.Val < 0.05 & results2$logFC > 1)))
            sum(results2$adj.P.Val < 0.05 & results2$logFC < -1)))
            sum(results2$Is_Target_Gene & results2$adj.P.Val < 0.05)))
write.csv(results2, "limma_DEA_BC_vs_Normal_EMTAB779.csv", row.names = FALSE)
results2 <- results2 %>%
  mutate(
    Direction = case_when(
      logFC >= 1 & Significant == "Yes" ~ "Up",
      logFC <= -1 & Significant == "Yes" ~ "Down",
      TRUE ~ "Stable"
    )
  )
my_colors2 <- c("Up" = "#BC3C29", "Down" = "#0072B5", "Stable" = "#E6E6E6")
volcano_plot2 <- ggplot(results2, aes(x = logFC, y = negLog10P)) +
  geom_point(data = subset(results2, Direction == "Stable"),
             color = my_colors2["Stable"], alpha = 0.4, size = 0.8) +
  geom_point(data = subset(results2, Direction != "Stable" & Is_Target_Gene == FALSE),
             aes(color = Direction), alpha = 0.3, size = 1.2) +
  geom_point(data = subset(results2, Is_Target_Gene == TRUE),
             aes(fill = Direction), shape = 21, color = "black", 
             size = 2.8, stroke = 0.6, alpha = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray60", linewidth = 0.4) +
  geom_text_repel(
    data = subset(results2, Is_Target_Gene == TRUE),
    aes(label = Gene_Symbol),
    size = 3.5,
    fontface = "bold.italic",
    box.padding = 0.6,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey30",
    segment.linewidth = 0.4,
    force = 8
  ) +
  scale_color_manual(values = my_colors2) +
  scale_fill_manual(values = my_colors2) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1, color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "E-MTAB-779",
    x = expression(bold(log[2]("Fold Change"))),
    y = expression(bold(-log[10]("FDR")))
  ) +
  scale_x_continuous(expand = expansion(mult = 0.1))
ggsave("volcano_plot_BC_vs_Normal_EMTAB779.pdf", plot = volcano_plot2, width = 10, height = 8)
library(tidyr)
target_genes_expr1 <- expr_final[rownames(expr_final) %in% gene_set, ]
available_genes1 <- intersect(gene_set, rownames(expr_final))
if (length(available_genes1) < length(gene_set)) {
  missing_genes <- setdiff(gene_set, rownames(expr_final))
}
target_genes_expr1_df <- as.data.frame(t(target_genes_expr1))
target_genes_expr1_df$Sample <- rownames(target_genes_expr1_df)
target_genes_expr1_df$Group <- ifelse(sample_disease[valid_idx] == "inflammatory breast cancer", "IBC", "Normal")
expr_long1 <- target_genes_expr1_df %>%
  pivot_longer(cols = all_of(available_genes1), names_to = "Gene", values_to = "Expression")
boxplot1 <- ggplot(expr_long1, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(aes(color = Group), 
               fill = "white", 
               outlier.shape = NA,      
               linewidth = 0.8,         
               width = 0.6,             
               position = position_dodge(0.8)) +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", 
                     method = "t.test",    
                     hide.ns = TRUE,       
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", " "))) +
  scale_color_manual(values = c("IBC" = "#BC3C29", "Normal" = "#0072B5")) +
  theme_classic(base_size = 14) + 
  theme(
    axis.line = element_line(linewidth = 1.0, color = "black"), 
    axis.ticks = element_line(linewidth = 1.0, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "plain"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(size = 14, face = "plain"),
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    legend.position = "none" 
  )+
  labs(
    title = "E-GEOD-45581",
    x=NULL,
    y = "Expression"
  )
ggsave("boxplot_ADH_ALDH_EGEOD45581.pdf", plot = boxplot1, width = 10, height = 6)
target_genes_expr2 <- expr_final2[rownames(expr_final2) %in% gene_set, ]
available_genes2 <- intersect(gene_set, rownames(expr_final2))
if (length(available_genes2) < length(gene_set)) {
  missing_genes2 <- setdiff(gene_set, rownames(expr_final2))
}
target_genes_expr2_df <- as.data.frame(t(target_genes_expr2))
target_genes_expr2_df$Sample <- rownames(target_genes_expr2_df)
target_genes_expr2_df$Group <- ifelse(sample_disease2[valid_idx2] == "breast carcinoma", "BC", "Normal")
expr_long2 <- target_genes_expr2_df %>%
  pivot_longer(cols = all_of(available_genes2), names_to = "Gene", values_to = "Expression")
boxplot2 <- ggplot(expr_long2, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(aes(color = Group), 
               fill = "white", 
               outlier.shape = NA,      
               linewidth = 0.8,         
               width = 0.6,             
               position = position_dodge(0.8)) +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", 
                     method = "t.test",    
                     hide.ns = TRUE,       
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", " "))) +
  scale_color_manual(values = c("IBC" = "#BC3C29", "Normal" = "#0072B5")) +
  theme_classic(base_size = 14) + 
  theme(
    axis.line = element_line(linewidth = 1.0, color = "black"), 
    axis.ticks = element_line(linewidth = 1.0, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "plain"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(size = 14, face = "plain"),
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    legend.position = "none" 
  ) +
  labs(
    title = "E-MTAB-779",
    x=NULL,
    y = "Expression"
  )
ggsave("boxplot_ADH_ALDH_EMTAB779.pdf", plot = boxplot2, width = 10, height = 6)
library(rstatix)
boxplot_stats1 <- expr_long1 %>%
  group_by(Gene) %>%
  t_test(Expression ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
print(boxplot_stats1[, c("Gene", "p", "p.adj", "p.adj.signif")])
boxplot_stats2 <- expr_long2 %>%
  group_by(Gene) %>%
  t_test(Expression ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
print(boxplot_stats2[, c("Gene", "p", "p.adj", "p.adj.signif")])
gene_order <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7", 
                "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")
stats1_wide <- boxplot_stats1 %>%
  select(Gene, p.adj, p.adj.signif) %>%
  rename(p.adj_EGEOD45581 = p.adj, signif_EGEOD45581 = p.adj.signif)
stats2_wide <- boxplot_stats2 %>%
  select(Gene, p.adj, p.adj.signif) %>%
  rename(p.adj_EMTAB779 = p.adj, signif_EMTAB779 = p.adj.signif)
combined_stats_wide <- full_join(stats1_wide, stats2_wide, by = "Gene")
complete_gene_df <- data.frame(Gene = gene_order)
combined_stats_wide <- left_join(complete_gene_df, combined_stats_wide, by = "Gene")
combined_stats_wide <- combined_stats_wide %>%
  select(Gene, p.adj_EGEOD45581, p.adj_EMTAB779)
colnames(combined_stats_wide) <- c("Gene", "p.adj_EGEOD45581", "p.adj_EMTAB779")
print(combined_stats_wide)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Statistics")
combined_stats_wide_formatted <- combined_stats_wide %>%
  mutate(
    p.adj_EGEOD45581 = ifelse(is.na(p.adj_EGEOD45581), 
                               "NA", 
                               formatC(p.adj_EGEOD45581, format = "e", digits = 2)),
    p.adj_EMTAB779 = ifelse(is.na(p.adj_EMTAB779), 
                             "NA", 
                             formatC(p.adj_EMTAB779, format = "e", digits = 2))
  )
colnames(combined_stats_wide_formatted) <- c("Gene", "E-GEOD-45581", "E-MTAB-779")
writeData(wb, "Statistics", combined_stats_wide_formatted, startRow = 1, startCol = 1)
headerStyle <- createStyle(
  fontSize = 12,
  fontColour = "white",
  fgFill = "#4472C4",
  halign = "center",
  valign = "center",
  textDecoration = "bold"
)
addStyle(wb, "Statistics", headerStyle, rows = 1, cols = 1:3, gridExpand = TRUE)
redStyle <- createStyle(
  fontColour =  "#FF0000",
  fontSize = 11,
  halign = "center"
)
blackStyle <- createStyle(
  fontColour =  "#000000",
  fontSize = 11,
  halign = "center"
)
for (i in 1:nrow(combined_stats_wide)) {
  row_num <- i + 1  
  if (!is.na(combined_stats_wide[i,2]) && combined_stats_wide[i,2] < 0.05) {
    addStyle(wb, "Statistics", redStyle, rows = row_num, cols = 2, gridExpand = TRUE)
  } else {
    addStyle(wb, "Statistics", blackStyle, rows = row_num, cols = 2, gridExpand = TRUE)
  }
  if (!is.na(combined_stats_wide[i,3]) && combined_stats_wide[i,3] < 0.05) {
    addStyle(wb, "Statistics", redStyle, rows = row_num, cols = 3, gridExpand = TRUE)
  } else {
    addStyle(wb, "Statistics", blackStyle, rows = row_num, cols = 3, gridExpand = TRUE)
  }
}
setColWidths(wb, "Statistics", cols = 1, widths = 12)
setColWidths(wb, "Statistics", cols = 2:3, widths = 18)
saveWorkbook(wb, "boxplot_statistics_combined.xlsx", overwrite = TRUE)
combined_boxplot <- boxplot1 + boxplot2 + 
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )
ggsave("boxplot_ADH_ALDH_combined.pdf", plot = combined_boxplot, width = 16, height = 6)
