library(switchBox)
drug_res=readRDS("Outputs/RelevantDrugRes.RDS")
drugs=unique(drug_res$inhibitor)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]
scores_single=list()

for (drug in drugs){ 
  rel_res=subset(drug_res, inhibitor==drug)
  inc_samples=intersect(rel_res$dbgap_rnaseq_sample, colnames(geneexp))
  alls=geneexp[,inc_samples]
  all_labels=unlist(rel_res[match(colnames(alls), rel_res$dbgap_rnaseq_sample), "SensitivityCall"])
  all_labels <- factor(all_labels, levels = c("Resistant", "Sensitive"))
  if (length(table(all_labels))==2 & all(table(all_labels)>=20)){ 
    x<-  SWAP.Train.KTSP(as.matrix(alls), as.factor(all_labels), krange=1:15)
    ktspStatDefault <- SWAP.KTSP.Statistics(inputMat = as.matrix(alls), classifier = x)
    scores_single[[drug]]=cbind.data.frame(Drug= drug, 
                             Sample=names(ktspStatDefault$statistics),Stat=ktspStatDefault$statistics,
                             SensitivityCall= all_labels)
  }
}
dir.create("PlotData", showWarnings = FALSE)

saveRDS(scores_single, file="PlotData/single_drug_scores.RDS")
