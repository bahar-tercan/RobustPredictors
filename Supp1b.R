library(ggplot2)
library(readxl)
library(readr)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(tibble)
rel_drugs=c()
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

drug_data=combo_drugs
i=0
for (drug1 in drugs){  
  drug_status=subset(drug_data, drug==drug1)
  inc_samples_wave12=intersect(waves1_2, intersect(drug_status$beataml_dbgap_rnaseq_sample, colnames(geneexp)))
  inc_samples_wave34=intersect(waves3_4, intersect(drug_status$beataml_dbgap_rnaseq_sample, colnames(geneexp)))
  w12_labels=unlist(drug_status[match(inc_samples_wave12,drug_status$beataml_dbgap_rnaseq_sample),
                                "SensitivityCall"])
  w34_labels=unlist(drug_status[match(inc_samples_wave34,drug_status$beataml_dbgap_rnaseq_sample), 
                                "SensitivityCall"])
  w12_labels=w12_labels[!is.na(w12_labels)]  
  w34_labels=w34_labels[!is.na(w34_labels)]
  if (length(table(w12_labels))==2 & all(table(w12_labels)>=20)){
    if (length(table(w34_labels))==2 & all(table(w34_labels)>=10)){
      i=i+1
      rel_drugs[i]=drug1
    }
  }
}
# --- 4. Final Visualization ---
rel_drug_plot <- subset(combo_drugs, drug %in% rel_drugs)
rel_drug_plot$auc <- as.numeric(rel_drug_plot$probit_auc)

Figure1b <- ggplot(rel_drug_plot, aes(x = drug, y = auc, color = Cohort)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_manual(values = c("#66A61EFF", "#E6AB02FF")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 12),
        legend.position = "bottom") +
  labs(x = "Drug", 
       y = "Area Under the Curve (AUC)")

# View plot
saveRDS(file="PlotData/supp_figure1b.RDS",  Figure1b)
