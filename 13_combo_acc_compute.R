library(caret)
library(readr)
library(tibble)
library(randomForest)
library(readxl)
library(reshape2)
library(dplyr)
library(glmnet)
library(pROC)
library(mltools)
library(e1071)

cli=read_excel("../DataSets/beataml_wv1to4_clinical.xlsx")

geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
ava_s=colnames(geneexp)

vene_comb=read_excel("../DataSets/bcd-23-0014_table_s9_suppst9.xlsx",col_types='text')
vene_comb=subset(vene_comb, dxAtSpecimenAcquisition=="ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS")
vene_comb$probit_auc=as.numeric(vene_comb$probit_auc)
vene_comb=mutate(vene_comb, SensitivityCall=ifelse(probit_auc>100, "Resistant", "Sensitive"))
combo_drugs=subset(vene_comb, drug_type=="combo")
f=unique(combo_drugs[!(combo_drugs$beataml_dbgap_rnaseq_sample%in%ava_s),c("beataml_dbgap_subject_id",
                                                                           "beataml_dbgap_rnaseq_sample")])
zz=cli[cli$dbgap_subject_id%in%unlist(f[complete.cases(f),"beataml_dbgap_subject_id"]),"dbgap_rnaseq_sample"]
xxx=unique(unlist(combo_drugs[which(is.na(combo_drugs$beataml_dbgap_rnaseq_sample)),
                              "beataml_dbgap_subject_id"]))
xx=subset(cli, dbgap_subject_id%in%xxx)
yy=xx[!is.na(xx$dbgap_rnaseq_sample),]
ttm=yy[yy$manuscript_rnaseq=="yes",c("dbgap_subject_id", "dbgap_rnaseq_sample")]
id=setdiff(yy$dbgap_subject_id, ttm$dbgap_subject_id)
m=yy[yy$dbgap_subject_id%in%id,c("dbgap_subject_id","dbgap_rnaseq_sample")]
m=m[m$dbgap_rnaseq_sample%in%ava_s, c("dbgap_subject_id", "dbgap_rnaseq_sample")]
to_imp=rbind(ttm, m)

for (i in 1:nrow(to_imp)){ 
  combo_drugs[combo_drugs$beataml_dbgap_subject_id==to_imp[i,1], "beataml_dbgap_rnaseq_sample"]=to_imp[i,2]
}
combo_drugs=subset(combo_drugs, !is.na(beataml_dbgap_rnaseq_sample))
combo_drugs$Cohort=unlist(cli[match(combo_drugs$beataml_dbgap_rnaseq_sample, cli$dbgap_rnaseq_sample), "cohort"])
combo_drugs=subset(combo_drugs, Cohort%in%c("Waves1+2", "Waves3+4"))

waves1_2=unlist(combo_drugs[combo_drugs$Cohort=="Waves1+2","beataml_dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(combo_drugs[combo_drugs$Cohort=="Waves3+4","beataml_dbgap_rnaseq_sample"]);names(waves3_4)=NULL
waves1_2=waves1_2[!is.na(waves1_2)]
waves3_4=waves3_4[!is.na(waves3_4)]

drugs=unique(combo_drugs$drug)
j=0

GetACC=function(folder, pth){
  i=0;
  res=list()
  for (drug1 in drugs){  
    drug_status=subset(combo_drugs, drug==drug1)
    inc_samples_wave12=intersect(waves1_2, intersect(drug_status$beataml_dbgap_rnaseq_sample, colnames(geneexp)))
    inc_samples_wave34=intersect(waves3_4, intersect(drug_status$beataml_dbgap_rnaseq_sample, colnames(geneexp)))
    w12_labels=unlist(drug_status[match(inc_samples_wave12, drug_status$beataml_dbgap_rnaseq_sample), "SensitivityCall"])
    w34_labels=unlist(drug_status[match(inc_samples_wave34, drug_status$beataml_dbgap_rnaseq_sample), "SensitivityCall"])
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


elastic_nobal=GetACC("ComboDrugModels/","_elasticnet_no.RDS")
elastic_bal=GetACC("ComboDrugModels/","_elasticnet_yes.RDS")

linearsvm_nobal=GetACC("ComboDrugModels/","_linearsvm_no.RDS")
linearsvm_bal=GetACC("ComboDrugModels/","_linearsvm_yes.RDS")

rbfsvm_nobal=GetACC("ComboDrugModels/","_rbfsvm_no.RDS")
rbfsvm_bal=GetACC("ComboDrugModels/","_rbfsvm_yes.RDS")

randomforest_nobal=GetACC("ComboDrugModels/","_randomforest_no.RDS")
randomforest_bal=GetACC("ComboDrugModels/","_randomforest_yes.RDS")

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

combo_acc_other=rbind(elastic_nobal,elastic_bal,linearsvm_nobal,linearsvm_bal, 
                       rbfsvm_nobal, rbfsvm_bal, randomforest_nobal,randomforest_bal)

dir.create("Accuracies", showWarnings = FALSE)
saveRDS(combo_acc_other, file="Accuracies/combo_other_classifiers.RDS")
