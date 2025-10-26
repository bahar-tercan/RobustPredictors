library(caret)
library(readr)
library(tibble)
library(randomForest) # REQUIRED for method = "rf"
library(glmnet)       # REQUIRED for method = "glmnet"
library(e1071)        # REQUIRED for method = "svmLinear" and "svmRadial"
library(themis)       # REQUIRED for sampling = "smote"

Param_tuning=function(trainData, balancing, model, drug){ 
 if (balancing=="no"){ 
   ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    search="random"
  )
 }
  if (balancing=="yes"){ 
    ctrl <- trainControl(
      method = "cv",
      number = 5,
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      search="random",
      sampling="smote"
    )
  }
  if (model=="elasticnet"){ 
    trained_model <- train(
      Class ~ .,
      data = trainData,
      method = "glmnet",
      trControl = ctrl,
      tuneLength = 10,
      preProcess = c("center", "scale"),
      metric = "ROC"
    )
  }
  if (model=="linearsvm"){ 
       trained_model <- train(Class ~ ., 
         data = trainData,
         method = "svmLinear",    
         trControl = ctrl,
         tuneLength = 10,
         metric = "ROC",          
         preProcess = c("center", "scale"))
  }
  if (model=="rbfsvm"){ 
    trained_model <- train(
      Class ~ ., 
      data = trainData,
      method = "svmRadial", 
      trControl = ctrl,
      preProcess = c("center", "scale"),
      metric = "ROC" ,
      tuneLength = 10
    )
  }
  if (model=="randomforest"){ 
    trained_model <- train(
      Class ~ .,
      data = trainData,
      method = "rf",
      trControl = ctrl,
      metric = "ROC",
      tuneLength = 10,
    )
  } 
  filename=paste0("SingleDrugModels/", drug, "_", model, "_", balancing, ".RDS")
  saveRDS(trained_model, file=filename)
}

drug_res=readRDS("Outputs/RelevantDrugRes.RDS")
drugs=unique(drug_res$inhibitor)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
waves1_2=unlist(cli_data[cli_data$cohort=="Waves1+2","dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(cli_data[cli_data$cohort=="Waves3+4","dbgap_rnaseq_sample"]);names(waves3_4)=NULL
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
dir.create("SingleDrugModels")

drugs=unique(drug_res$inhibitor)
for (drug1 in drugs){  
  drug_status=subset(drug_res, inhibitor==drug1)
  inc_samples_wave12=intersect(waves1_2, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
  inc_samples_wave34=intersect(waves3_4, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
  w12_labels=unlist(drug_status[match(inc_samples_wave12, drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
  w34_labels=unlist(drug_status[match(inc_samples_wave34, drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
  if (length(table(w12_labels))==2 & all(table(w12_labels)>=20)){ 
    if (length(table(w34_labels))==2 & all(table(w34_labels)>=10)){ 
      wave12_data=geneexp[,inc_samples_wave12]
      train_data=cbind.data.frame(t(wave12_data), Class=w12_labels)
      train_data$Class <- as.factor(train_data$Class)
      train_data$Class <- relevel(train_data$Class, ref = "Resistant")
      Param_tuning(train_data, balancing="no",  model="randomforest", drug1)
      Param_tuning(train_data, balancing="yes", model="randomforest", drug1)
      Param_tuning(train_data, balancing="no",  model="linearsvm", drug1)
      Param_tuning(train_data, balancing="yes", model="linearsvm", drug1)
      Param_tuning(train_data, balancing="no",  model="rbfsvm", drug1)
      Param_tuning(train_data, balancing="yes", model="rbfsvm", drug1)
      Param_tuning(train_data, balancing="no",  model="elasticnet", drug1)
      Param_tuning(train_data, balancing="yes", model="elasticnet", drug1)
      }
  }
}
