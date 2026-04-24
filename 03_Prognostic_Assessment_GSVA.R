# ==============================================================================
# 03_Prognostic_Assessment_GSVA.R
# Description: GSVA scoring for prognostic evaluation (log-rank test)
# ==============================================================================

library(GSVA)
library(survival)
library(survminer)
library(ggplot2)
library(reshape2)

# --- Step 1: Load data ---
load("./data/raw/TCGA/TCGA_exp.rda")
gsea_results <- readRDS("./data/processed/gsea_screening_results.rds")
sig_subpathways <- gsea_results$old_results$Subpathway

subpathway_gene <- read.csv("./data/processed/subpathway_gene.csv", header = TRUE, row.names = 1)
subpathway_list <- split(subpathway_gene$gene, subpathway_gene$name)
subpathway_list <- subpathway_list[sig_subpathways]

# --- Step 2: GSVA scoring ---
group <- ifelse(substr(colnames(exp), 14, 15) <= "10" & substr(colnames(exp), 1, 4) == "TCGA", "tumor", "normal")
gsva_scores <- gsva(as.matrix(exp), subpathway_list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
gsva_mat <- as.data.frame(t(gsva_scores))
gsva_mat$group <- group

# --- Step 3: Survival analysis for each subpathway ---
load("./data/raw/TCGA/data_clinical_patient.txt")
clinical <- data_clinical_patient[, c("X.Patient.Identifier", "Overall.Survival.Status", "Overall.Survival..Months.")]
colnames(clinical) <- c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")
clinical$OS_MONTHS <- as.numeric(clinical$OS_MONTHS)
clinical$Event <- ifelse(clinical$OS_STATUS == "1:DECEASED", 1, 0)

surv_results <- data.frame()
for (spw in sig_subpathways) {
  score <- gsva_scores[spw, ]
  names(score) <- colnames(gsva_scores)
  
  surv_data <- data.frame(
    Sample = names(score),
    Score = score,
    Patient_ID = substr(names(score), 1, 12)
  )
  surv_data <- merge(surv_data, clinical, by.x = "Patient_ID", by.y = "PATIENT_ID", all.x = TRUE)
  surv_data <- surv_data[!is.na(surv_data$OS_MONTHS) & !is.na(surv_data$Event), ]
  
  if (nrow(surv_data) < 10) next
  
  cutpoint <- surv_cutpoint(surv_data, time = "OS_MONTHS", event = "Event", variables = "Score")
  best_cut <- cutpoint$cutpoint$cutpoint
  surv_data$Risk <- ifelse(surv_data$Score >= best_cut, "High", "Low")
  surv_data$Risk <- factor(surv_data$Risk, levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS_MONTHS, Event) ~ Risk, data = surv_data)
  cox_fit <- coxph(Surv(OS_MONTHS, Event) ~ Risk, data = surv_data)
  cox_sum <- summary(cox_fit)
  p_val <- cox_sum$logtest["pvalue"]
  
  surv_results <- rbind(surv_results, data.frame(
    Subpathway = spw,
    Cutpoint = best_cut,
    HR = cox_sum$conf.int[1],
    P_Value = p_val
  ))
}

# --- Step 4: Save ---
saveRDS(list(gsva_mat = gsva_mat, surv_results = surv_results), "./data/processed/gsva_survival_results.rds")
write.csv(surv_results, "./results/prognostic_assessment_results.csv", row.names = FALSE)

cat("Prognostic assessment complete.\n")
