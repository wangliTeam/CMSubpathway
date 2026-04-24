# ==============================================================================
# 04_Classification_Assessment_SVM.R
# Description: SVM classification performance evaluation (AUC) for subpathways
# ==============================================================================

library(e1071)
library(pROC)
library(ggplot2)

# --- Step 1: Load GSVA scores ---
gsva_data <- readRDS("./data/processed/gsva_survival_results.rds")
gsva_mat <- gsva_data$gsva_mat
exp_data <- t(gsva_mat[, -ncol(gsva_mat)])

# --- Step 2: Define groups ---
group_labels <- ifelse(substr(rownames(exp_data), 14, 15) <= "10", 1, 0)
exp_df <- as.data.frame(exp_data)
exp_df$group <- as.factor(group_labels)

# --- Step 3: SVM + AUC for each subpathway ---
sig_subpathways <- colnames(exp_data)
auc_results <- data.frame()

for (spw in sig_subpathways) {
  x <- exp_df[, c(spw, "group")]
  colnames(x)[1] <- "spw"
  y <- as.factor(x$group)
  
  set.seed(123)
  idx <- sample(nrow(x), nrow(x) * 0.75)
  x_train <- x[idx, ]
  y_train <- y[idx]
  x_test <- x[-idx, ]
  y_test <- y[-idx]
  
  obj <- tune(svm, group ~ spw, data = x_train,
              ranges = list(gamma = 2^(-10:3), cost = 2^(-5:3)),
              tunecontrol = tune.control(sampling = "fix"))
  
  model_best <- svm(group ~ spw, x_train, scale = FALSE, kernel = "radial",
                    gamma = obj$best.parameters$gamma,
                    cost = obj$best.parameters$cost, probability = TRUE)
  
  y_pred <- predict(model_best, x_test, probability = TRUE)
  y_pred <- as.ordered(y_pred)
  
  svm_roc <- roc(response = x_test$group, predictor = y_pred, levels = levels(x_test$group))
  auc_val <- round(as.numeric(auc(svm_roc)), 4)
  
  auc_results <- rbind(auc_results, data.frame(Subpathway = spw, AUC = auc_val))
}

# --- Step 4: Save ---
saveRDS(auc_results, "./data/processed/svm_auc_results.rds")
write.csv(auc_results, "./results/classification_assessment_results.csv", row.names = FALSE)

cat("Classification assessment complete.\n")
