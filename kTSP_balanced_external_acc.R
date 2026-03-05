library(tibble)
library(org.Hs.eg.db)
library(switchBox)
library(readr)
library(readxl)
library(ggpubr)
library(caret)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(splitTools)
library(mltools)
library(rstatix)
library(stringr)
library(themis)


geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
drug_res=readRDS("Outputs/RelevantDrugRes.RDS")

fimm=read.csv("GeneExpMatrices/fimm_logged_tpm.csv")
fimm=column_to_rownames(fimm, var="X")
aml_dat= fimm[,sapply(colnames(fimm), startsWith, "AML")]
inc_genes=intersect(rownames(geneexp), rownames(aml_dat))
aml_dat=aml_dat[inc_genes, ]
geneexp=geneexp[inc_genes, ]
path="SingleDrugModels/"
pths=list.files(path)
drugs=unique(sapply(pths, function(x){str_split(x, "_")}[[1]][1]))
ext_drugs=as.data.frame(read_excel("../DataSets/File_3.2_Drug_response_DSS_sDSS_164S_17Healthy.xlsx", skip=2))
ext_drugs=ext_drugs[!is.na(ext_drugs$sDSS),]
ext_drugs <- ext_drugs %>%
  mutate(Sensitivity_Call = if_else(sDSS < 8.7, "Resistant", "Sensitive"))

all_models=list.files(path)
results=list()
i=0
all_models=all_models[!sapply(all_models, startsWith, "Cytarabine")]
drugs=unique(sapply(all_models, function(x){str_split(x,"_")[[1]][1]}))

models=list()
test_labels_list=list()
i=0
for (drug in drugs){
  i=i+1
  drug_status=subset(drug_res, inhibitor==drug)
  inc_samples= intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp))
  geneexp2=geneexp[,inc_samples]
  lvls <- c("Resistant", "Sensitive")
  all_labels=unlist(drug_status[match(colnames(geneexp2), drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
  all_labels=factor(all_labels, levels = lvls)
  rel_ext=subset(ext_drugs, Chemical_compound==drug)
  dd=intersect(rel_ext$Sample_ID, colnames(aml_dat))
  test_dat=aml_dat[,dd]
  test_labels=unlist(rel_ext[match(dd,rel_ext$Sample_ID),"Sensitivity_Call"])
  actuals <- factor(test_labels, levels = lvls)
  dd=intersect(rel_ext$Sample_ID, colnames(aml_dat))
  test_labels_list[[i]]=test_labels
  
  dd2=cbind.data.frame(t(geneexp2),  labels=as.factor(all_labels))
  mm <- smote(dd2, var = "labels", over_ratio = 1) 
  geneexp2<- t(mm[, !(names(mm) %in% "labels")])
  all_labels=unlist(mm[,"labels"])
  
  x<- SWAP.Train.KTSP(as.matrix(geneexp2),
                      as.factor(all_labels), krange=1:15)
  models[[drug]]=x
  res=SWAP.GetKTSP.Result(x, as.matrix(test_dat), test_labels,
                          c("Sensitive", "Resistant"), predictions=TRUE, decision_values=TRUE)
  results[[i]]=res
}
names(results)=drugs

alls=list()
for (i in 1:length(results)){
  accs=results[[i]]$stats
  conf=caret::confusionMatrix(data=as.factor(results[[i]]$predictions), 
                              reference=as.factor(test_labels_list[[i]]),
                              positive="Sensitive")
  alls[[i]]= accs 
}
accuracies=do.call(rbind, alls)
rownames(accuracies)=names(results)

saveRDS(accuracies, file="External_Accuracies/kTSP_single_balanced.RDS")
saveRDS(lapply(test_labels_list2,table), file="External_Accuracies/ext_response_dist.RDS")
