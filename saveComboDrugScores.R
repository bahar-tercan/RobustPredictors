library(readxl)
library(tibble)
library(org.Hs.eg.db)
library(readr)
library(ggpubr)
library(caret)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(splitTools)
library(switchBox)
library(mltools)
library(rstatix)
library(gridExtra)
library(reshape2)

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


drugs=unique(combo_drugs$drug)

scores_combo=list()
for (drug1 in drugs){ 
  rel_res=subset(combo_drugs, drug==drug1)
  inc_samples=intersect(rel_res$beataml_dbgap_rnaseq_sample, colnames(geneexp))
  alls=geneexp[,inc_samples]
  all_labels=unlist(rel_res[match(colnames(alls), rel_res$beataml_dbgap_rnaseq_sample), 
                            "SensitivityCall"])
  if (length(table(all_labels))==2 & all(table(all_labels)>=20)){ 
    x<-  SWAP.Train.KTSP(as.matrix(alls), as.factor(all_labels), krange=1:15)
    ktspStatDefault <- SWAP.KTSP.Statistics(inputMat = as.matrix(alls), classifier = x)
    scores_combo[[drug1]]=cbind.data.frame(Drug= drug1, 
                                           Sample=names(ktspStatDefault$statistics),
                                           Stat=ktspStatDefault$statistics,
                                           SensitivityCall= all_labels)
  }
}
dir.create("PlotData", showWarnings = FALSE) 

saveRDS(scores_combo, file="PlotData/combo_drug_scores.RDS")
