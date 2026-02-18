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

drug_res=readRDS("Outputs/RelevantDrugRes.RDS")
drugs=unique(drug_res$inhibitor)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
geneexp=read.csv("GeneExpMatrices/beataml_logged_tpm.csv", check.names=FALSE)
rownames(geneexp)=geneexp[,1]
geneexp=geneexp[,-1]

fimm=read.csv("GeneExpMatrices/fimm_logged_tpm.csv")
fimm=column_to_rownames(fimm, var="X")
aml_dat= fimm[,sapply(colnames(fimm), startsWith, "AML")]
inc_genes=intersect(rownames(geneexp), rownames(aml_dat))
aml_dat=aml_dat[inc_genes, ]
geneexp=geneexp[inc_genes, ]
results=list()
for (drug1 in drugs){
  drug_status=subset(drug_res, inhibitor==drug1)
  inc_samples= intersect(drug_status$dbgap_rnaseq_sample, colnames(geneexp))
  geneexp2=geneexp[,inc_samples]
  all_labels=unlist(drug_status[match(colnames(geneexp2), drug_status$dbgap_rnaseq_sample), "SensitivityCall"])
  if (length(table(all_labels))==2 & all(table(all_labels)>=20)){ 
    x<- SWAP.Train.KTSP(as.matrix(geneexp2),
                        as.factor(all_labels), krange=1:15)
    res=SWAP.KTSP.Classify(as.matrix(aml_dat),  x)
    results[[drug1]]=as.character(res)
  }
}


datt <- list()
for (i in seq_along(names(results))){
  drug_name <- names(results)[i]
  datt[[i]] <- data.frame(
    Sample = colnames(aml_dat),
    Drug = drug_name,
    PredictedLabel = results[[drug_name]],
    stringsAsFactors = FALSE
  )
}
ddd <- do.call(rbind, datt)

ext_drugs=as.data.frame(read_excel("DataSets/File_3.2_Drug_response_DSS_sDSS_164S_17Healthy.xlsx", skip=2))
merged=merge(ext_drugs, ddd, by.x=c("Sample_ID", "Chemical_compound") , by.y=c("Sample", "Drug") )
na_counts <- merged %>%
  group_by(Chemical_compound) %>%
  summarise(
    na_count = sum(is.na(sDSS))
  )

merged=subset(merged, Chemical_compound!="Cytarabine")
#Cytarabine's all values are NA 

merged$PredictedLabel <- factor(merged$PredictedLabel, levels = c("Resistant", "Sensitive"))

stat.test <- merged %>%
  group_by(Chemical_compound) %>%
  wilcox_test(sDSS ~ PredictedLabel, alternative = "less") %>%
  adjust_pvalue() %>%
  add_significance("p.adj")%>%
  add_xy_position(x = "PredictedLabel")

# Fix: Correct save path
dir.create("Accuracies")
saveRDS(stat.test, file = "Accuracies/fimm_validation_stats.RDS")
custom_colors <- c("Resistant" = "#8DA0CB", "Sensitive" = "#E78AC3")

bxp <- ggboxplot(merged, x = "PredictedLabel", y = "sDSS", fill = "PredictedLabel",
                 facet.by = "Chemical_compound",  nrow = 2 )
plt=bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj") +
  scale_fill_manual(values=custom_colors) +ylim(c(-20, 50))+
  theme(legend.position="bottom") +
  theme(strip.text = element_text(size = 8),
        # Rotate x-axis labels 45 degrees
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.ticks.x = element_line() # Restore the axis ticks
  ) + xlab("") + labs(fill = "Predicted Sensitivity Call")
ggsave(plt, filename = "../Plots/Figure4.pdf", dpi=600,height=5)

