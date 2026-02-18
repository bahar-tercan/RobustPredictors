library(readxl)
library(readr)
library(reshape2)
library(dplyr)
library(splitTools)
library(switchBox)
library(caret)
library(mltools)


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



combo_model=list()
for (drug1 in drugs){ 
  rel_res=subset(combo_drugs, drug==drug1)
  inc_samples=intersect(rel_res$beataml_dbgap_rnaseq_sample, colnames(geneexp))
  alls=geneexp[,inc_samples]
  all_labels=unlist(rel_res[match(colnames(alls), rel_res$beataml_dbgap_rnaseq_sample), "SensitivityCall"])
  if (length(table(all_labels))==2 & all(table(all_labels)>=20)){ 
    x<-  SWAP.Train.KTSP(as.matrix(alls), as.factor(all_labels), krange=1:15)
    combo_model[[drug1]]=x
  }
}
saveRDS(combo_model, file="Outputs/combodrug_models.RDS")

