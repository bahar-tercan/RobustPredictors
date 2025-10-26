library(splitTools)
library(switchBox)
library(themis)
library(mltools)
library(caret)

drug_res=readRDS("Outputs/RelevantDrugRes.RDS")
drugs=unique(drug_res$inhibitor)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
waves1_2=unlist(cli_data[cli_data$cohort=="Waves1+2","dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(cli_data[cli_data$cohort=="Waves3+4","dbgap_rnaseq_sample"]);names(waves3_4)=NULL
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
#accuracies=list()
#rel_inds=matrix(NA, nrow=length(drugs) , ncol=5)

results=list()
models=list()
train_labels=list()
test_labels=list()
for (drug1 in drugs){
  drug_status=subset(drug_res, inhibitor==drug1)
  inc_samples_wave12=intersect(waves1_2, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
  inc_samples_wave34=intersect(waves3_4, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
  w12_labels=unlist(drug_status[match(inc_samples_wave12, drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
  w34_labels=unlist(drug_status[match(inc_samples_wave34, drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
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
                              reference=as.factor(test_labels[[i]]),
                              positive="Sensitive")
  f1_measure=conf$byClass[c("F1","Pos Pred Value","Neg Pred Value")]    
  mcc_value <- mcc(preds = as.factor(results[[i]]$predictions), 
                   actuals =as.factor(test_labels[[i]]))
  alls[[i]]=c(accs, f1_measure, MCC=mcc_value)  
}
accuracies=do.call(rbind, alls)
rownames(accuracies)=names(results)
#dir.create("Accuracies")
saveRDS(accuracies, file="Accuracies/kTSP_single_balanced.RDS")



