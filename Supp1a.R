library(ggplot2)
library(readxl)
library(readr)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(tibble)
rel_drugs=c()

cli_data=readRDS("Outputs/RelevantClinical.RDS")
drug_data=readRDS("Outputs/RelevantDrugRes.RDS")
dat=read.csv("GeneExpMatrices/beataml_logged_tpm.csv")
dat=column_to_rownames(dat, var="X")
geneexp=dat
waves1_2=unlist(cli_data[cli_data$cohort=="Waves1+2","dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(cli_data[cli_data$cohort=="Waves3+4","dbgap_rnaseq_sample"]);names(waves3_4)=NULL
waves1_2=waves1_2[!is.na(waves1_2)]
waves3_4=waves3_4[!is.na(waves3_4)]
drugs=unique(drug_data$inhibitor)


i=0
for (drug in drugs){  
  drug_status=subset(drug_data, inhibitor==drug)
  inc_samples_wave12=intersect(waves1_2, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
  inc_samples_wave34=intersect(waves3_4, intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp)))
  w12_labels=drug_status[match(inc_samples_wave12,drug_status$dbgap_rnaseq_sample), "SensitivityCall"]
  w34_labels=drug_status[match(inc_samples_wave34,drug_status$dbgap_rnaseq_sample), "SensitivityCall"]
  if (length(table(w12_labels))==2 & all(table(w12_labels)>=20)){
    if (length(table(w34_labels))==2 & all(table(w34_labels)>=10)){
      i=i+1
      rel_drugs[i]=drug
    }
  }
}
drug_data=subset(drug_data,inhibitor%in%rel_drugs)
rel_drug=merge(drug_data,cli_data, by="dbgap_rnaseq_sample")
rel_drug$auc=as.numeric(rel_drug$auc) 
colnames(rel_drug)[colnames(rel_drug)=="cohort"] ="Cohort"

Figure1a=ggplot(rel_drug, aes(x=inhibitor ,y= auc, color=Cohort)) +geom_boxplot() +
  scale_color_manual(values=c("#66A61EFF", "#E6AB02FF")) +
  theme_bw() + scale_x_discrete(guide = guide_axis(angle = 90))  + xlab("Inhibitor") +
  theme(text = element_text(family = "Times New Roman", size=12)) +
  ylab("Area Under the Curve (AUC)") + 
  xlab("Drug")+
  theme(legend.position="bottom") 
saveRDS(file="PlotData/supp_figure1a.RDS",  Figure1a)
  


