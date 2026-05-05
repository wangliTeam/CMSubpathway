library(ggplot2)
library(GSVA)
library(survival)
library(survminer)
library(org.Hs.eg.db)
library(dplyr)
library(patchwork)
load("./step1/data/TCGA_exp.rda")
data_timeline_treatment <- read.delim('./step1/data/data_timeline_treatment.txt', stringsAsFactors = FALSE)
data_clinical_patient <- read.delim('./step1/data/data_clinical_patient.txt', stringsAsFactors = FALSE)
gene_set <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7", "ALDH2", "ALDH1B1", "ALDH9A1", "ALDH3A2", "ALDH7A1")
sample_names <- colnames(data_mrna_seq_tpm)[-1]
patient_ids <- substr(sample_names, 1, 12)
is_tumor <- grep(".01A$", sample_names)
data_mrna_seq_tpm <- data_mrna_seq_tpm[, c(1, is_tumor+1)]
patient_ids <- patient_ids[is_tumor]
entrez_ids <- data_mrna_seq_tpm$Entrez_Gene_Id
entrez_ids <- as.character(entrez_ids)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = entrez_ids, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")
valid_idx <- !is.na(gene_symbols)
data_mrna_seq_tpm <- data_mrna_seq_tpm[valid_idx, ]
gene_symbols <- gene_symbols[valid_idx]
dup_symbols <- duplicated(gene_symbols)
if (any(dup_symbols)) {
  data_mrna_seq_tpm$gene_symbol <- gene_symbols
  data_mrna_seq_tpm_agg <- aggregate(. ~ gene_symbol, data = data_mrna_seq_tpm, FUN = mean)
  rownames(data_mrna_seq_tpm_agg) <- data_mrna_seq_tpm_agg$gene_symbol
  data_mrna_seq_tpm <- data_mrna_seq_tpm_agg[, -1]
  gene_symbols <- rownames(data_mrna_seq_tpm)
}
rownames(data_mrna_seq_tpm) <- gene_symbols
data_mrna_seq_tpm <- data_mrna_seq_tpm[, colnames(data_mrna_seq_tpm) != 'Entrez_Gene_Id']
available_genes <- intersect(gene_set, rownames(data_mrna_seq_tpm))

gene_set_list <- list(Subpathway = available_genes)
gsvaPar <- gsvaParam(
  exprData = as.matrix(data_mrna_seq_tpm),
  geneSets = gene_set_list,
  minSize = 1,
  kcdf = "Gaussian"
)
gsva_scores <- gsva(gsvaPar, verbose = FALSE)
subpathway_score <- as.numeric(gsva_scores)
names(subpathway_score) <- colnames(gsva_scores)
clinical_data <- data_clinical_patient %>%
  select(X.Patient.Identifier, Overall.Survival.Status, Overall.Survival..Months.) %>%
  filter(!is.na(Overall.Survival.Status) & !is.na(Overall.Survival..Months.))
colnames(clinical_data) <- c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")
clinical_data$OS_MONTHS <- as.numeric(as.character(clinical_data$OS_MONTHS))
clinical_data$Event <- ifelse(clinical_data$OS_STATUS == "1:DECEASED", 1, 0)
chemo_data <- data_timeline_treatment %>%
  select(PATIENT_ID, THERAPEUTIC_AGENT, TREATMENT_TYPE)
chemo_data$Is_Chemotherapy <- chemo_data$TREATMENT_TYPE == "Chemotherapy"
anthracycline_agents <- c("Doxorubicin", "Epirubicin", "Doxorubicin Hydrochloride", 
                          "Pegylated Liposomal Doxorubicin Hydrochloride")
taxane_agents <- c("Paclitaxel", "Docetaxel", "Nab-paclitaxel", "Tesetaxel")
alkylating_agents <- c("Cyclophosphamide", "Ifosfamide")
platinum_agents <- c("Carboplatin", "Cisplatin")
other_chemo_agents <- c("Capecitabine", "Fluorouracil", "Gemcitabine", "Gemcitabine Hydrochloride",
                        "Methotrexate", "Vinblastine", "Vincristine", "Vinorelbine", "Vinorelbine Tartrate",
                        "Etoposide", "Mitomycin", "Mitoxantrone", "Ixabepilone", "Trabectedin",
                        "Pemetrexed", "Lapatinib", "Everolimus")
chemo_data$Is_Chemotherapy <- chemo_data$TREATMENT_TYPE == "Chemotherapy"
chemo_data$Anthracycline <- chemo_data$THERAPEUTIC_AGENT %in% anthracycline_agents & chemo_data$Is_Chemotherapy
chemo_data$Taxane <- chemo_data$THERAPEUTIC_AGENT %in% taxane_agents & chemo_data$Is_Chemotherapy
chemo_data$Alkylating <- chemo_data$THERAPEUTIC_AGENT %in% alkylating_agents & chemo_data$Is_Chemotherapy
chemo_data$Platinum <- chemo_data$THERAPEUTIC_AGENT %in% platinum_agents & chemo_data$Is_Chemotherapy
chemo_summary <- chemo_data %>%
  group_by(PATIENT_ID) %>%
  summarise(
    Received_Any_Chemo = max(Is_Chemotherapy),  
    Anthracycline = max(Anthracycline),
    Taxane = max(Taxane),
    Alkylating = max(Alkylating),
    Platinum = max(Platinum),
    .groups = 'drop'
  )
chemo_summary$Chemo_Group <- with(chemo_summary, case_when(
  Received_Any_Chemo == FALSE ~ "Chemotherapy-naive",
  Anthracycline == TRUE & Taxane == TRUE ~ "AC-T (Anthracycline+Taxane)",
  Anthracycline == TRUE & Alkylating == TRUE ~ "AC (Anthracycline+Cyclophosphamide)",
  Anthracycline == TRUE ~ "Anthracycline-based",
  Taxane == TRUE & Platinum == TRUE ~ "TP (Taxane+Platinum)",
  Taxane == TRUE ~ "Taxane-based",
  Platinum == TRUE ~ "Platinum-based",
  Alkylating == TRUE ~ "Cyclophosphamide-based",
  TRUE ~ "Other Chemotherapy"
))
analysis_data <- data.frame(
  Patient_ID = patient_ids,
  Sample_ID = names(subpathway_score),
  Subpathway_Score = subpathway_score,
  stringsAsFactors = FALSE
)
analysis_data <- merge(analysis_data, clinical_data, 
                       by.x = "Patient_ID", by.y = "PATIENT_ID", 
                       all.x = TRUE)
chemo_join <- chemo_summary[, c("PATIENT_ID", "Chemo_Group")]
analysis_data <- merge(analysis_data, chemo_join, 
                       by.x = "Patient_ID", by.y = "PATIENT_ID", 
                       all.x = TRUE)
analysis_data <- analysis_data[!is.na(analysis_data$OS_STATUS) & !is.na(analysis_data$Subpathway_Score), ]
analysis_data$Chemo_Group[is.na(analysis_data$Chemo_Group)] <- "Chemotherapy-naive"
chemo_table <- table(analysis_data$Chemo_Group)
print(chemo_table)
for (group in names(chemo_table)) {
              group, 
              chemo_table[group], 
              chemo_table[group] / sum(chemo_table) * 100))
}
naive_patients <- analysis_data[analysis_data$Chemo_Group == "Chemotherapy-naive", ]
if (nrow(naive_patients) > 0) {
  naive_patient_ids <- naive_patients$Patient_ID
  naive_treatment_data <- data_timeline_treatment[data_timeline_treatment$PATIENT_ID %in% naive_patient_ids, ]
  if (nrow(naive_treatment_data) > 0) {
    treatment_type_table <- table(naive_treatment_data$TREATMENT_TYPE)
    print(treatment_type_table)
    agent_table <- sort(table(naive_treatment_data$THERAPEUTIC_AGENT), decreasing = TRUE)
    top_agents <- head(agent_table, 20)
    for (i in 1:length(top_agents)) {
      agent_name <- names(top_agents)[i]
      count <- as.numeric(top_agents[i])
      percentage <- count / nrow(naive_treatment_data) * 100
                  substr(agent_name, 1, 50), 
                  count, 
                  percentage))
    }
    if (length(agent_table) > 20) {
    }
    naive_treatment_stats <- data.frame(
      Treatment_Agent = names(agent_table),
      Count = as.numeric(agent_table),
      Percentage = as.numeric(agent_table) / sum(agent_table) * 100
    )
    write.csv(naive_treatment_stats, "chemotherapy_naive_treatment_details.csv", row.names = FALSE)
  } else {
  }
} else {
}
plot_survival_curve <- function(data, chemo_group, output_prefix) {
  subgroup_data <- data[data$Chemo_Group == chemo_group, ]
  if (nrow(subgroup_data) < 10) {
    return(NULL)
  }
  subgroup_data$OS_MONTHS <- as.numeric(as.character(subgroup_data$OS_MONTHS))
  subgroup_data$Event <- as.numeric(subgroup_data$Event)
  subgroup_data <- subgroup_data[!is.na(subgroup_data$OS_MONTHS) & !is.na(subgroup_data$Event) & !is.na(subgroup_data$Subpathway_Score), ]
  if (nrow(subgroup_data) < 10) {
    return(NULL)
  }
  cutpoint_result <- surv_cutpoint(
    data = subgroup_data,
    time = "OS_MONTHS",
    event = "Event",
    variables = "Subpathway_Score"
  )
  best_cutpoint <- cutpoint_result$cutpoint$cutpoint
  subgroup_data$Risk_Group <- ifelse(subgroup_data$Subpathway_Score >= best_cutpoint, "High", "Low")
  subgroup_data$Risk_Group <- factor(subgroup_data$Risk_Group, levels = c("Low", "High"))
  subgroup_data$surv_obj <- Surv(subgroup_data$OS_MONTHS, subgroup_data$Event)
  fit <- survfit(surv_obj ~ Risk_Group, data = subgroup_data)
  cox_fit <- coxph(surv_obj ~ Risk_Group, data = subgroup_data)
  cox_summary <- summary(cox_fit)
  hr <- cox_summary$conf.int[1]
  hr_lower <- cox_summary$conf.int[3]
  hr_upper <- cox_summary$conf.int[4]
  p_value <- cox_summary$logtest["pvalue"]
  p <- ggsurvplot(fit, 
                  data = subgroup_data,
                  pval = TRUE,
                  pval.coord = c(0, 0.2),  
                  conf.int = TRUE,  
                  risk.table = FALSE,
                  legend.title = "Group",
                  legend.labs = c("low-score", "high-score"),
                  palette = c("#194280", "#c22f23"),  # blue and red
                  title = chemo_group,
                  xlab = "Time",
                  ylab = "Survival Probability",
                  ggtheme = theme_bw() +
                    theme(
                      panel.grid.major = element_blank(),  
                      panel.grid.minor = element_blank(),  
                      panel.border = element_blank(),      
                      axis.line = element_line(color = "black", size = 0.5),  
                      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                      legend.position = "top"
                    ),
                  table.height = 0.25)
  return(list(
    stats = data.frame(
      Chemo_Group = chemo_group,
      N = nrow(subgroup_data),
      High_N = sum(subgroup_data$Risk_Group == "High"),
      Low_N = sum(subgroup_data$Risk_Group == "Low"),
      Events = sum(subgroup_data$Event),
      HR = hr,
      HR_Lower = hr_lower,
      HR_Upper = hr_upper,
      P_Value = p_value,
      Best_Cutpoint = best_cutpoint
    ),
    plot = p$plot
  ))
}
results_list <- list()
plot_list <- list()
chemo_groups <- unique(analysis_data$Chemo_Group)
analysis_data$OS_MONTHS <- as.numeric(analysis_data$OS_MONTHS)
for (group in chemo_groups) {
  result <- plot_survival_curve(analysis_data, group, "KM_curve")
  if (!is.null(result)) {
    results_list[[group]] <- result$stats
    plot_list[[group]] <- result$plot
  }
}
all_chemo_data <- analysis_data[analysis_data$Chemo_Group != "Chemotherapy-naive", ]
if (nrow(all_chemo_data) >= 10) {
  cutpoint_all_chemo <- surv_cutpoint(
    data = all_chemo_data,
    time = "OS_MONTHS",
    event = "Event",
    variables = "Subpathway_Score"
  )
  best_cutpoint_all_chemo <- cutpoint_all_chemo$cutpoint$cutpoint
  all_chemo_data$Risk_Group <- ifelse(all_chemo_data$Subpathway_Score >= best_cutpoint_all_chemo, "High", "Low")
  all_chemo_data$Risk_Group <- factor(all_chemo_data$Risk_Group, levels = c("Low", "High"))
  surv_obj_chemo <- Surv(all_chemo_data$OS_MONTHS, all_chemo_data$Event)
  fit_chemo <- survfit(surv_obj_chemo ~ Risk_Group, data = all_chemo_data)
  cox_chemo <- coxph(surv_obj_chemo ~ Risk_Group, data = all_chemo_data)
  cox_summary_chemo <- summary(cox_chemo)
  hr_chemo <- cox_summary_chemo$conf.int[1]
  hr_lower_chemo <- cox_summary_chemo$conf.int[3]
  hr_upper_chemo <- cox_summary_chemo$conf.int[4]
  p_value_chemo <- cox_summary_chemo$logtest["pvalue"]
  p_chemo <- ggsurvplot(fit_chemo,
                        data = all_chemo_data,
                        pval = TRUE,
                        pval.coord = c(0, 0),  
                        conf.int = FALSE,  
                        risk.table = TRUE,
                        legend.title = "Group",
                        legend.labs = c("low-score", "high-score"),
                        palette = c("#194280", "#c22f23"),
                        title = "All Chemotherapy Patients",
                        xlab = "Time (Months)",
                        ylab = "Overall Survival Probability",
                        ggtheme = theme_bw() +
                          theme(
                            panel.grid.major = element_blank(),  
                            panel.grid.minor = element_blank(),  
                            panel.border = element_blank(),      
                            axis.line = element_line(color = "black", size = 0.5),  
                            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                            legend.position = "top"
                          ),
                        table.height = 0.25)
  results_list[["All_Chemotherapy"]] <- data.frame(
    Chemo_Group = "All Chemotherapy Patients",
    N = nrow(all_chemo_data),
    High_N = sum(all_chemo_data$Risk_Group == "High"),
    Low_N = sum(all_chemo_data$Risk_Group == "Low"),
    Events = sum(all_chemo_data$Event),
    HR = hr_chemo,
    HR_Lower = hr_lower_chemo,
    HR_Upper = hr_upper_chemo,
    P_Value = p_value_chemo,
    Best_Cutpoint = best_cutpoint_all_chemo
  )
  plot_list[["All_Chemotherapy"]] <- p_chemo$plot
} else {
}
if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 2) + 
    plot_annotation(
      title = 'Survival Analysis by Chemotherapy Subgroups',
      subtitle = paste('Analyzed on:', Sys.Date()),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  ggsave("KM_curve_All_Chemo_Groups.pdf", 
         plot = combined_plot, 
         width = 12, 
         height = 6 * ceiling(length(plot_list) / 2),
         limitsize = FALSE)
} else {
}
results_df <- do.call(rbind, results_list)
print(results_df)
chemo_summary_stats <- analysis_data %>%
  group_by(Chemo_Group) %>%
  summarise(
    Total_Patients = n(),
    Events = sum(Event, na.rm = TRUE),
    Event_Rate = sum(Event, na.rm = TRUE) / n() * 100,
    Mean_Score = mean(Subpathway_Score, na.rm = TRUE),
    Median_Score = median(Subpathway_Score, na.rm = TRUE)
  )
print(as.data.frame(chemo_summary_stats))
significant_groups <- results_df[results_df$P_Value < 0.05, ]
if (nrow(significant_groups) > 0) {
  for (i in 1:nrow(significant_groups)) {
                i, 
                significant_groups$Chemo_Group[i],
                significant_groups$P_Value[i],
                significant_groups$HR[i]))
  }
  significant_plot_list <- list()
  for (group_name in significant_groups$Chemo_Group) {
    for (key in names(plot_list)) {
      if (grepl(group_name, key, fixed = TRUE) || group_name == key) {
        significant_plot_list[[key]] <- plot_list[[key]]
        break
      }
    }
  }
  if (length(significant_plot_list) > 0) {
    combined_plot <- wrap_plots(significant_plot_list, ncol = 2) + 
      plot_annotation(
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
      )
    ggsave("KM_curve_Significant_Groups.pdf", 
           plot = combined_plot, 
           width = 8, 
           height = 4 * ceiling(length(significant_plot_list) / 2),
           limitsize = FALSE)
  }
} else {
}
write.csv(results_df, "survival_analysis_all_results.csv", row.names = FALSE)
write.csv(as.data.frame(chemo_summary_stats), "chemo_group_statistics.csv", row.names = FALSE)
write.csv(analysis_data, "analysis_data_with_scores.csv", row.names = FALSE)