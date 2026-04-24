# ==============================================================================
# 05_BRCA_Validation_GSE42568.R
# Description: Independent BRCA dataset validation (GSE42568) and Subpathway-CorSP 
#              benchmark comparison
# ==============================================================================

library(GEOquery)
library(GSVA)
library(survival)
library(survminer)
library(pROC)
library(e1071)

# --- Step 1: Load GSE42568 ---
gse <- getGEO("GSE42568", GSEMatrix = TRUE)[[1]]
expr <- exprs(gse)
pdata <- pData(gse)

# Define tumor/normal groups based on sample characteristics
group <- ifelse(tolower(pdata$characteristics_ch1.1) == "tumor", "tumor", "normal")

# --- Step 2: Load subpathway gene sets ---
subpathway_gene <- read.csv("./data/processed/subpathway_gene.csv", header = TRUE, row.names = 1)
subpathway_list <- split(subpathway_gene$gene, subpathway_gene$name)

# --- Step 3: GSVA scoring ---
gsva_scores <- gsva(expr, subpathway_list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

# --- Step 4: Survival analysis for BRCA risk subpathways ---
brca_risk_subpathways <- c("hsa00071-2", "hsa00620-2", "hsa00830-1")
surv_results <- data.frame()

for (spw in brca_risk_subpathways) {
  if (!spw %in% rownames(gsva_scores)) next
  
  score <- gsva_scores[spw, ]
  surv_data <- data.frame(Score = score, Sample = names(score))
  
  # Match with clinical data
  surv_data$OS_MONTHS <- as.numeric(pdata[match(names(score), rownames(pdata)), "geo_accession"])
  surv_data$Event <- as.numeric(pdata[match(names(score), rownames(pdata)), "status"])
  surv_data <- surv_data[!is.na(surv_data$OS_MONTHS) & !is.na(surv_data$Event), ]
  
  if (nrow(surv_data) < 10) next
  
  cutpoint <- surv_cutpoint(surv_data, time = "OS_MONTHS", event = "Event", variables = "Score")
  best_cut <- cutpoint$cutpoint$cutpoint
  surv_data$Risk <- factor(ifelse(surv_data$Score >= best_cut, "High", "Low"), levels = c("Low", "High"))
  
  fit <- survfit(Surv(OS_MONTHS, Event) ~ Risk, data = surv_data)
  cox_fit <- coxph(Surv(OS_MONTHS, Event) ~ Risk, data = surv_data)
  p_val <- summary(cox_fit)$logtest["pvalue"]
  
  surv_results <- rbind(surv_results, data.frame(Subpathway = spw, Cutpoint = best_cut, P_Value = p_val))
}

# --- Step 5: SVM classification ---
auc_results <- data.frame()
for (spw in brca_risk_subpathways) {
  if (!spw %in% rownames(gsva_scores)) next
  
  x <- data.frame(spw = gsva_scores[spw, ], group = factor(group))
  set.seed(123)
  idx <- sample(nrow(x), nrow(x) * 0.75)
  
  model <- svm(group ~ spw, x[idx, ], scale = FALSE, kernel = "radial", probability = TRUE)
  pred <- predict(model, x[-idx, ], probability = TRUE)
  
  roc_obj <- roc(response = x$group[-idx], predictor = as.numeric(as.ordered(pred)))
  auc_results <- rbind(auc_results, data.frame(Subpathway = spw, AUC = round(as.numeric(auc(roc_obj)), 4)))
}

# --- Step 6: Save ---
saveRDS(list(gsva_scores = gsva_scores, surv_results = surv_results, auc_results = auc_results),
        "./data/processed/gse42568_validation.rds")
write.csv(surv_results, "./results/gse42568_survival.csv", row.names = FALSE)
write.csv(auc_results, "./results/gse42568_auc.csv", row.names = FALSE)

cat("GSE42568 validation complete.\n")
