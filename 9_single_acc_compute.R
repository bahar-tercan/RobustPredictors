library(caret)
library(readr)
library(tibble)
library(randomForest)
library(readxl)
library(reshape2)
library(dplyr)
library(caret)
library(readr)
library(tibble)
library(glmnet)
library(pROC)
library(mltools)

drug_res=readRDS("Outputs/RelevantDrugRes.RDS")
drugs=unique(drug_res$inhibitor)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
waves1_2=unlist(cli_data[cli_data$cohort=="Waves1+2","dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(cli_data[cli_data$cohort=="Waves3+4","dbgap_rnaseq_sample"]);names(waves3_4)=NULL
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]

GetACC=function(folder, pth){
  i=0;
  res=list()
  for (drug1 in drugs){  
    drug_status=subset(drug_res, inhibitor==drug1)
    inc_samples_wave12=intersect(waves1_2, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
    inc_samples_wave34=intersect(waves3_4, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
    w12_labels=unlist(drug_status[match(inc_samples_wave12, drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
    w34_labels=unlist(drug_status[match(inc_samples_wave34, drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
    if (length(table(w12_labels))==2 & all(table(w12_labels)>=20)){ 
      if (length(table(w34_labels))==2 & all(table(w34_labels)>=10)){    
        i=i+1
        wave12_data=geneexp[,inc_samples_wave12]
        wave34_data=geneexp[,inc_samples_wave34]
        #train_data=cbind.data.frame(t(wave12_data), Class=w12_labels)
        #train_data$Class <- as.factor(train_data$Class)
        #train_data$Class <- relevel(train_data$Class, ref = "Resistant")
       if(file.exists(paste0(folder, drug1, pth))){ 
          model=readRDS(paste0(folder, drug1, pth))
          pred_class=predict(model, t(wave34_data))  
          pred_class=relevel(pred_class, ref = "Resistant")
          pred_prob  <- predict(model, t(wave34_data), type = "prob")[,"Sensitive"]
          conf <- confusionMatrix(data=as.factor(pred_class), 
                                reference=relevel(as.factor(w34_labels), ref = "Resistant"),
                                positive = "Sensitive")
         # print(conf)
          others=conf$byClass[c("Sensitivity", "Specificity",
                      "Balanced Accuracy", "F1",
                      "Pos Pred Value","Neg Pred Value")]
        
          roc_obj <- roc(response = as.factor(w34_labels),
                       predictor = pred_prob)
        auc=auc(roc_obj)
        mcc_value <- mcc(preds = as.factor(pred_class), 
                         actuals =as.factor(w34_labels))
        
        res[[i]]=c(Drug=drug1, others, AUC=auc, MCC=mcc_value)
        }
      }
    }
  }
  res=do.call(rbind.data.frame,res)
  colnames(res)=c("Drug","Sensitivity", "Specificity",
                  "Balanced Accuracy", "F1",
                  "Pos Pred Value","Neg Pred Value", "AUC", "MCC" )
  return(res)
}  


elastic_nobal=GetACC("SingleDrugModels/","_elasticnet_no.RDS")
elastic_bal=GetACC("SingleDrugModels/","_elasticnet_yes.RDS")

linearsvm_nobal=GetACC("SingleDrugModels/","_linearsvm_no.RDS")
linearsvm_bal=GetACC("SingleDrugModels/","_linearsvm_yes.RDS")

rbfsvm_nobal=GetACC("SingleDrugModels/","_rbfsvm_no.RDS")
rbfsvm_bal=GetACC("SingleDrugModels/","_rbfsvm_yes.RDS")

randomforest_nobal=GetACC("SingleDrugModels/","_randomforest_no.RDS")
randomforest_bal=GetACC("SingleDrugModels/","_randomforest_yes.RDS")

elastic_nobal$Method="ElasticNet"
elastic_bal$Method="ElasticNet"
elastic_nobal$Balancing="Not Balanced"
elastic_bal$Balancing="Balanced"

linearsvm_nobal$Method="linearSVM"
linearsvm_bal$Method="linearSVM"
linearsvm_nobal$Balancing="Not Balanced"
linearsvm_bal$Balancing="Balanced"

rbfsvm_nobal$Method="rbfSVM"
rbfsvm_bal$Method="rbfSVM"
rbfsvm_nobal$Balancing="Not Balanced"
rbfsvm_bal$Balancing="Balanced"

randomforest_nobal$Method="Random Forest"
randomforest_bal$Method="Random Forest"
randomforest_nobal$Balancing="Not Balanced"
randomforest_bal$Balancing="Balanced"

single_acc_other=rbind(elastic_nobal,elastic_bal,linearsvm_nobal,linearsvm_bal, 
      rbfsvm_nobal, rbfsvm_bal, randomforest_nobal,randomforest_bal)

saveRDS(single_acc_other , file="Accuracies/SingleDrug_otheracc.RDS")
