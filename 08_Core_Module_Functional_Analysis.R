# ==============================================================================
# 08_Core_Module_Functional_Analysis.R
# Description: Core module functional characteristics analysis including immune 
#              infiltration, chemotherapy subgroup survival, and pathway enrichment
# ==============================================================================

library(GSVA)
library(survival)
library(survminer)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# --- Step 1: Load TCGA-BRCA data ---
load("./data/raw/TCGA/TCGA_exp.rda")
load("./data/raw/TCGA/data_clinical_patient.txt")
load("./data/raw/TCGA/data_timeline_treatment.txt")

core_module_genes <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH5", "ADH6", "ADH7",
                       "ALDH1B1", "ALDH2", "ALDH3A2", "ALDH7A1", "ALDH9A1")

# --- Step 2: GSVA scoring ---
gene_set_list <- list(CoreModule = core_module_genes)
gsva_scores <- gsva(as.matrix(exp), gene_set_list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
score <- as.numeric(gsva_scores)
names(score) <- colnames(gsva_scores)

# --- Step 3: Clinical data preparation ---
clinical <- data_clinical_patient[, c("X.Patient.Identifier", "Overall.Survival.Status", "Overall.Survival..Months.")]
colnames(clinical) <- c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")
clinical$OS_MONTHS <- as.numeric(clinical$OS_MONTHS)
clinical$Event <- ifelse(clinical$OS_STATUS == "1:DECEASED", 1, 0)

chemo_data <- data_timeline_treatment %>%
  select(PATIENT_ID, THERAPEUTIC_AGENT, TREATMENT_TYPE)
chemo_data$Is_Chemo <- chemo_data$TREATMENT_TYPE == "Chemotherapy"

anthracycline <- c("Doxorubicin", "Epirubicin")
taxane <- c("Paclitaxel", "Docetaxel")
chemo_data$Anthracycline <- chemo_data$THERAPEUTIC_AGENT %in% anthracycline & chemo_data$Is_Chemo
chemo_data$Taxane <- chemo_data$THERAPEUTIC_AGENT %in% taxane & chemo_data$Is_Chemo

chemo_summary <- chemo_data %>%
  dplyr::group_by(PATIENT_ID) %>%
  dplyr::summarise(Received_Chemo = max(Is_Chemo), Anthracycline = max(Anthracycline), Taxane = max(Taxane), .groups = "drop")

chemo_summary$Chemo_Group <- dplyr::case_when(
  Received_Chemo == FALSE ~ "Chemotherapy-naive",
  Anthracycline == TRUE & Taxane == TRUE ~ "AC-T",
  Anthracycline == TRUE ~ "Anthracycline",
  Taxane == TRUE ~ "Taxane",
  TRUE ~ "Other"
)

# --- Step 4: Survival analysis by chemotherapy subgroups ---
analysis_data <- data.frame(Patient_ID = substr(names(score), 1, 12), Score = score)
analysis_data <- merge(analysis_data, clinical, by.x = "Patient_ID", by.y = "PATIENT_ID", all.x = TRUE)
analysis_data <- merge(analysis_data, chemo_summary[, c("PATIENT_ID", "Chemo_Group")],
                       by.x = "Patient_ID", by.y = "PATIENT_ID", all.x = TRUE)
analysis_data$Chemo_Group[is.na(analysis_data$Chemo_Group)] <- "Chemotherapy-naive"

chemo_groups <- unique(analysis_data$Chemo_Group)
surv_results <- list()

for (grp in chemo_groups) {
  sub <- analysis_data[analysis_data$Chemo_Group == grp, ]
  sub <- sub[!is.na(sub$OS_MONTHS) & !is.na(sub$Event) & !is.na(sub$Score), ]
  if (nrow(sub) < 10) next
  
  cutpoint <- surv_cutpoint(sub, time = "OS_MONTHS", event = "Event", variables = "Score")
  best_cut <- cutpoint$cutpoint$cutpoint
  sub$Risk <- factor(ifelse(sub$Score >= best_cut, "High", "Low"), levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS_MONTHS, Event) ~ Risk, data = sub)
  cox_fit <- coxph(Surv(OS_MONTHS, Event) ~ Risk, data = sub)
  p_val <- summary(cox_fit)$logtest["pvalue"]
  hr <- summary(cox_fit)$conf.int[1]
  
  surv_results[[grp]] <- data.frame(Group = grp, N = nrow(sub), HR = hr, P_Value = p_val, Cutpoint = best_cut)
}

surv_df <- do.call(rbind, surv_results)

# --- Step 5: Immune infiltration analysis ---
high_group <- analysis_data[analysis_data$Score >= median(analysis_data$Score, na.rm = TRUE), ]
low_group <- analysis_data[analysis_data$Score < median(analysis_data$Score, na.rm = TRUE), ]

high_samples <- high_group$Patient_ID
low_samples <- low_group$Patient_ID

# Differential expression between high/low metabolic activity groups
high_expr <- exp[, substr(colnames(exp), 1, 12) %in% high_samples]
low_expr <- exp[, substr(colnames(exp), 1, 12) %in% low_samples]

design <- model.matrix(~0 + factor(c(rep("High", ncol(high_expr)), rep("Low", ncol(low_expr)))))
colnames(design) <- c("High", "Low")
fit <- lmFit(exp[, c(colnames(high_expr), colnames(low_expr))], design)
fit2 <- contrasts.fit(fit, makeContrasts(High - Low, levels = design))
fit2 <- eBayes(fit2)
deg <- topTable(fit2, n = Inf, sort.by = "P")

# Pathway enrichment
deg_sig <- deg[deg$adj.P.Val < 0.05, ]
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(deg_sig), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
ego <- enrichGO(gene = na.omit(entrez_ids), OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.25)

# --- Step 6: Save ---
saveRDS(list(surv_results = surv_df, deg = deg, go_enrichment = ego),
        "./data/processed/core_module_functional.rds")
write.csv(surv_df, "./results/chemo_subgroup_survival.csv", row.names = FALSE)

cat("Core module functional analysis complete.\n")
