library(readxl)
library(tibble)
library(org.Hs.eg.db)
library(switchBox)
library(readr)
library(ggpubr)
library(caret)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(splitTools)
library(mltools)
library(rstatix)
library(gridExtra)
drug_res=readRDS("Outputs/RelevantDrugRes.RDS")
drugs=unique(drug_res$inhibitor)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
single_model=list()
for (drug in drugs){ 
  rel_res=subset(drug_res, inhibitor==drug)
  inc_samples=intersect(rel_res$dbgap_rnaseq_sample, colnames(geneexp))
  alls=geneexp[,inc_samples]
  all_labels=unlist(rel_res[match(colnames(alls), rel_res$dbgap_rnaseq_sample), "SensitivityCall"])
  if (length(table(all_labels))==2 & all(table(all_labels)>=20)){ 
    x<-  SWAP.Train.KTSP(as.matrix(alls), as.factor(all_labels), krange=1:15)
    single_model[[drug]]=x
  }
}
saveRDS(single_model, file="Outputs/singledrug_models.RDS")
