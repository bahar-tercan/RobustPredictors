library(caret)
library(readr)
library(tibble)
library(readxl)    # REQUIRED for read_excel()
library(dplyr)     # REQUIRED for mutate()
library(randomForest) # REQUIRED for method = "rf"
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
  filename=paste0("ComboDrugModels/", drug, "_", model, "_", balancing, ".RDS")
  saveRDS(trained_model, file=filename)
}
dir.create("ComboDrugModels")
cli=read_excel("DataSets/beataml_wv1to4_clinical.xlsx")
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
vene_comb=read_excel("DataSets/bcd-23-0014_table_s9_suppst9.xlsx",col_types='text')
vene_comb=subset(vene_comb, dxAtSpecimenAcquisition=="ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS")
vene_comb$probit_auc=as.numeric(vene_comb$probit_auc)
vene_comb=mutate(vene_comb, SensitivityCall=ifelse(probit_auc<100, "Sensitive", "Resistant"))
combo_drugs=subset(vene_comb, drug_type=="combo")
combo_drugs=subset(combo_drugs, !is.na(beataml_dbgap_rnaseq_sample))
combo_drugs$Cohort=unlist(cli[match(combo_drugs$beataml_dbgap_rnaseq_sample, cli$dbgap_rnaseq_sample), "cohort"])
combo_drugs=subset(combo_drugs, Cohort%in%c("Waves1+2", "Waves3+4"))

waves1_2=unlist(combo_drugs[combo_drugs$Cohort=="Waves1+2","beataml_dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(combo_drugs[combo_drugs$Cohort=="Waves3+4","beataml_dbgap_rnaseq_sample"]);names(waves3_4)=NULL
waves1_2=unique(waves1_2[!is.na(waves1_2)])
waves3_4=unique(waves3_4[!is.na(waves3_4)])
drugs=unique(combo_drugs$drug)
train_sens_call_num=list()
test_sens_call_num=list()
drug_data=combo_drugs
j=0
z=list()
results=list()
models=list()
train_labels=list()
test_labels=list()

for (drug1 in drugs){  
  j=j+1
  drug_status=subset(combo_drugs, drug==drug1)
  inc_samples_wave12=intersect(waves1_2, intersect(drug_status$beataml_dbgap_rnaseq_sample, colnames(geneexp)))
  inc_samples_wave34=intersect(waves3_4, intersect(drug_status$beataml_dbgap_rnaseq_sample, colnames(geneexp)))
  w12_labels=unlist(drug_status[match(inc_samples_wave12, drug_status$beataml_dbgap_rnaseq_sample), "SensitivityCall"])
  w34_labels=unlist(drug_status[match(inc_samples_wave34, drug_status$beataml_dbgap_rnaseq_sample), "SensitivityCall"])
  if (length(table(w12_labels))==2 & all(table(w12_labels)>=20)){ 
    if (length(table(w34_labels))==2 & all(table(w34_labels)>=10)){ 
      wave12_data=geneexp[,inc_samples_wave12]
      wave34_data=geneexp[,inc_samples_wave34]
      train_labels[[drug1]]=w12_labels
      test_labels[[drug1]]=w34_labels
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
