library(tibble)
library(org.Hs.eg.db)
library(switchBox)
library(readr)
library(readxl)
library(ggpubr)
library(caret)
library(pROC)         
library(mltools)      
library(dplyr)
library(stringr)
library(rstatix)

geneexp <- read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names = FALSE)
rownames(geneexp) <- geneexp[, 1]
geneexp <- geneexp[, -1]
fimm <- read.csv("GeneExpMatrices/fimm_logged_tpm.csv")
fimm <- column_to_rownames(fimm, var = "X")

aml_dat <- fimm[, sapply(colnames(fimm), startsWith, "AML")]
inc_genes <- intersect(rownames(geneexp), rownames(aml_dat))
aml_dat <- aml_dat[inc_genes, ]

ext_drugs <- as.data.frame(read_excel("../DataSets/File_3.2_Drug_response_DSS_sDSS_164S_17Healthy.xlsx", skip = 2))
ext_drugs <- ext_drugs[!is.na(ext_drugs$sDSS), ]

ext_drugs <- ext_drugs %>%
  mutate(Sensitivity_Call = if_else(sDSS < 8.7, "Resistant", "Sensitive"))

GetACC = function(folder, pth_suffix) {
  available_files <- list.files(folder, pattern = pth_suffix)
  if(length(available_files) == 0) return(NULL)
    current_drugs <- sapply(available_files, function(x) str_split(x, "_")[[1]][1])
  
  res <- list()
  for (drug in current_drugs) {
    model_path <- paste0(folder, drug, pth_suffix)
    model <- readRDS(model_path)
    rel_ext <- subset(ext_drugs, Chemical_compound == drug)
    dd <- intersect(rel_ext$Sample_ID, colnames(aml_dat))
    test_dat <- t(aml_dat[, dd])
    test_labels <- rel_ext$Sensitivity_Call[match(dd, rel_ext$Sample_ID)]
    lvls <- model$levels 
    actuals <- factor(test_labels, levels = lvls)
    
    pred_class <- predict(model, test_dat)
    pred_prob  <- predict(model, test_dat, type = "prob")[, lvls[2]]
    
    conf <- confusionMatrix(data = pred_class, 
                            reference = actuals, 
                            positive = lvls[2])
    
    others <- conf$byClass[c("Sensitivity", "Specificity", 
                             "Balanced Accuracy")]
    
    roc_obj <- roc(response = actuals, predictor = pred_prob, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    res[[drug]] <- data.frame(Drug = drug, t(others), AUC = auc_val, 
                              stringsAsFactors = FALSE)
  }
  
  if(length(res) == 0) return(NULL)
  return(do.call(rbind, res))
}

path <- "AllModelsExternal/SingleDrugModels/"
methods <- c("elasticnet", "linearsvm", "rbfsvm", "randomforest")
balancing_options <- c("no", "yes")
all_acc_list <- list()
for (m in methods) {
  for (b in balancing_options) {
    suffix <- paste0("_", m, "_", b, ".RDS")
    
    # Extract metrics for this model type
    df <- GetACC(path, suffix)
    
    if (!is.null(df)) {
      # Label metadata
      df$Method <- switch(m, 
                          "elasticnet" = "ElasticNet",
                          "linearsvm" = "linearSVM",
                          "rbfsvm" = "rbfSVM",
                          "randomforest" = "Random Forest")
      
      df$Balancing <- ifelse(b == "yes", "Balanced", "Not Balanced")
      
      all_acc_list[[paste(m, b, sep = "_")]] <- df
    }
  }
}

single_acc_other <- do.call(rbind, all_acc_list)

# --- 5. Save Results ---
dir.create("External_Accuracies", showWarnings = FALSE)
saveRDS(single_acc_other, file = "External_Accuracies/SingleDrug_otheracc.RDS")
