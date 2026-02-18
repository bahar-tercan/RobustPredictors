library(readxl)
library(dplyr)

`%notin%`=Negate(`%in%`)
cli_data=read_excel("../DataSets/beataml_wv1to4_clinical.xlsx")
rel_cli=subset(cli_data, used_manuscript_analyses=="yes"& 
   dxAtSpecimenAcquisition%in%c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"))

for_normal_samples=read_excel("../DataSets/beataml_waves1to4_sample_mapping.xlsx")
healthy_samples=for_normal_samples[sapply(for_normal_samples$rna_control, startsWith, "Healthy"),]
healthy_samples=unlist(healthy_samples[!is.na(healthy_samples$rna_control),"dbgap_rnaseq_sample"])

rel_cli2=subset(rel_cli, rel_cli$dbgap_rnaseq_sample%notin%healthy_samples)


drug_response=read.csv2("../DataSets/beataml_probit_curve_fits_v4_dbgap.txt", sep="\t")
rel_drug_res=subset(drug_response, paper_inclusion=="TRUE")
rel_drug_res=rel_drug_res[,c("dbgap_subject_id", "dbgap_dnaseq_sample", 
                             "dbgap_rnaseq_sample", "inhibitor", "ic50", "auc")]
rel_drug_res$auc=as.numeric(rel_drug_res$auc)
rel_drug_res$ic50=as.numeric(rel_drug_res$ic50)
rel_drug_res1=mutate(rel_drug_res, SensitivityCall=ifelse(auc<100, "Sensitive", "Resistant"))
dir.create("Outputs")
saveRDS(rel_cli2,  file="Outputs/RelevantClinical.RDS")
saveRDS(rel_drug_res1,  file="Outputs/RelevantDrugRes.RDS")
