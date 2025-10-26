library(readxl)
library(readr)
library(reshape2)
library(dplyr)
library(splitTools)
library(switchBox)
library(caret)
library(mltools)
library(themis)
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
      dd=cbind.data.frame(t(wave12_data),  labels=as.factor(w12_labels))
      
      mm <- smote(dd, var = "labels", over_ratio = 1) 
      wave12_data<- t(mm[, !(names(mm) %in% "labels")])
      w12_labels=unlist(mm[,"labels"])
      x<- SWAP.Train.KTSP(as.matrix(wave12_data),
                          as.factor(w12_labels), krange=1:15)
      models[[drug1]]=x
      res=SWAP.GetKTSP.Result(x, as.matrix(wave34_data), w34_labels,
                              c("Sensitive", "Resistant"), 
                              predictions=TRUE, decision_values=TRUE)
      results[[drug1]]=res
    }
  }
}


alls=list()
for (i in 1:length(results)){
  accs=results[[i]]$stats
  conf=caret::confusionMatrix(data=as.factor(results[[i]]$predictions), 
                              reference=as.factor(test_labels[[i]]), positive="Sensitive")
  f1_measure=conf$byClass[c("F1","Pos Pred Value","Neg Pred Value")]    
  mcc_value <- mcc(preds = as.factor(results[[i]]$predictions), 
                   actuals =as.factor(test_labels[[i]]))
  alls[[i]]=c(accs, f1_measure, MCC=mcc_value)  
}
accuracies=do.call(rbind, alls)
rownames(accuracies)=names(results)
#dir.create("Accuracies")
saveRDS(accuracies, file="Accuracies/kTSP_balanced_combo.RDS")




