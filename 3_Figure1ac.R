library(ggplot2)
library(readxl)
library(readr)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(tibble)
rel_drugs=c()
i=0

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
train_sens_call_num=list()
test_sens_call_num=list()

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
      train_sens_call_num[[i]]=table(w12_labels)
      test_sens_call_num[[i]]=table(w34_labels)
    }
  }
}
train_calls= do.call(rbind, train_sens_call_num)
rownames(train_calls)=rel_drugs
test_calls = do.call(rbind, test_sens_call_num)
rownames(test_calls)=rel_drugs

train_calls=cbind(Cohort="Waves 1+2", melt(train_calls))
test_calls=cbind(Cohort="Waves 3+4", melt(test_calls))

all_sensitivity_calls=rbind(train_calls, test_calls)

colnames(all_sensitivity_calls)=c("Cohort", "Drug", "Sensitivity Call", "NumberSamples")
drug_data=subset(drug_data,inhibitor%in%rel_drugs)
rel_drug=merge(drug_data,cli_data, by.x="dbgap_subject_id", by.y= "dbgap_subject_id")
rel_drug$auc=as.numeric(rel_drug$auc) 
colnames(rel_drug)[colnames(rel_drug)=="cohort"] ="Cohort"

Figure1a=ggplot(rel_drug, aes(x=inhibitor ,y= auc, color=Cohort)) +geom_boxplot() +
  scale_color_manual(values=c("#66A61EFF", "#E6AB02FF")) +
  theme_bw() + scale_x_discrete(guide = guide_axis(angle = 90))  + xlab("Inhibitor") +
  theme(text = element_text(family = "Times New Roman", size=12)) +
    ylab("Area Under the Curve (AUC)") + 
  theme(legend.position="bottom") 

Figure1c <- ggplot(all_sensitivity_calls, aes(x = Drug, y = NumberSamples, fill = `Sensitivity Call`)) + 
  geom_col(position = position_dodge(width = 0.8)) +
  facet_grid(Cohort ~ Drug, scales = 'free_x', space = 'free_x') +
  labs(x = "", y = "Number of Samples", fill = "Sensitivity Call", title = "") +
  scale_fill_manual(values=c("#8DA0CB", "#E78AC3"))+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_blank(),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    #legend.box.background = element_rect(),
    legend.position = c(0.4, 0.23), legend.justification = c(0, 0)
  )
dir.create("PlotData")
saveRDS(Figure1a, file="PlotData/figure1a.RDS")
saveRDS(Figure1c, file="PlotData/figure1c.RDS")


