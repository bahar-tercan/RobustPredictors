library(ggplot2)
library(readxl)
library(readr)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ggpubr)

rel_drugs=c()
i=0
cli=read_excel("../DataSets/beataml_wv1to4_clinical.xlsx")
dat=read.csv("GeneExpMatrices/beataml_logged_tpm.csv")
dat=column_to_rownames(dat, var="X")
ava_s=colnames(dat)
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

geneexp=dat
waves1_2=unlist(combo_drugs[combo_drugs$Cohort=="Waves1+2","beataml_dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(combo_drugs[combo_drugs$Cohort=="Waves3+4","beataml_dbgap_rnaseq_sample"]);names(waves3_4)=NULL
waves1_2=waves1_2[!is.na(waves1_2)]
waves3_4=waves3_4[!is.na(waves3_4)]
drugs=unique(combo_drugs$drug)
train_sens_call_num=list()
test_sens_call_num=list()
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
drug_data=subset(drug_data,drug%in%rel_drugs)
#rel_drug=merge(drug_data,cli_data, by.x="dbgap_subject_id", by.y= "dbgap_subject_id")
rel_drug=drug_data
rel_drug$auc=as.numeric(rel_drug$probit_auc) 
#colnames(rel_drug)[colnames(rel_drug)=="cohort"] ="Cohort"

Figure1b=ggplot(rel_drug, aes(x=drug ,y= auc, color=Cohort)) +geom_boxplot() +
  scale_color_manual(values=c("#66A61EFF", "#E6AB02FF")) +
  theme_bw() + scale_x_discrete(guide = guide_axis(angle = 90))  + xlab("Inhibitor") +
  theme(text = element_text(family = "Times New Roman", size=12)) +
  ylab("Area Under the Curve (AUC)") + 
  theme(legend.position="bottom") 


Figure1d <- ggplot(all_sensitivity_calls, aes(x = Drug, y = NumberSamples, fill = `Sensitivity Call`)) + 
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
   # legend.position = c(0.5, 0.23), legend.justification = c(0, 0)
   legend.position="bottom"
  )

saveRDS(rel_drug, file="Outputs/combo_drugs.RDS")
saveRDS(Figure1b, file="PlotData/Figure1b.RDS")
saveRDS(Figure1d, file="PlotData/Figure1d.RDS")


