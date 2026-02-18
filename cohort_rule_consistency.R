library(tidyr)
library(dplyr)
library(ggplot2)

sd=readRDS("PlotData/single_drug_scores.RDS")
sd1=do.call(rbind, sd)
cd=readRDS("PlotData/combo_drug_scores.RDS")
cd1=do.call(rbind, cd)
alld=rbind(sd1, cd1)

cli_data=readRDS("Outputs/RelevantClinical.RDS")

merged=merge(alld, cli_data, by.x="Sample", by.y="dbgap_rnaseq_sample")
ddd=unique(merged$Drug)
drugs=list()
i=0
for (drug1 in ddd){ 
  xx=subset(merged, Drug==drug1)
  inc_samples_wave12=unique(unlist(subset(xx, cohort=="Waves1+2", select="Sample")))
  inc_samples_wave34=unique(unlist(subset(xx, cohort=="Waves3+4", select="Sample")))
  w12_labels=unlist(xx[match(inc_samples_wave12, xx$Sample), "SensitivityCall"])
  w34_labels=unlist(xx[match(inc_samples_wave34, xx$Sample), "SensitivityCall"])
  if (length(table(w12_labels))==2 & all(table(w12_labels)>=20)){ 
    if (length(table(w34_labels))==2 & all(table(w34_labels)>=20)){ 
      i=i+1
      drugs[[i]]=drug1 
    }
  }
}
merged_drug=subset(merged, Drug%in%unlist(drugs))
merged_drug$SenC=ifelse(merged_drug$SensitivityCall=="Sensitive", "S", "R")


merged_drug <- merged_drug %>%
  mutate(SensitivityCall = factor(SensitivityCall, levels = c("Resistant", "Sensitive")))

p_cohort <- ggplot(merged_drug, 
                   aes(x = SensitivityCall, # Compare S vs R on the x-axis
                       y = Stat, 
                       fill = cohort)) +
  geom_boxplot(outlier.size = 0.5) +
  # Facet by Drug so each inhibitor gets its own panel
  facet_wrap(~ Drug, scales = "free_y") + 
  theme_bw() +
  labs(x = "In-vitro Response", 
       y = "k-TSP Score (Sum of Votes)", 
       fill = "Dataset Cohort") +
  scale_fill_manual(values = c("#66A61EFF", "#E6AB02FF")) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(
      size = 6,             # Very small text for drug names
      margin = ggplot2::margin(1, 0, 1, 0), # Shrink the grey box around the name
      face = "bold"
    ),    axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(p_cohort, file="PlotData/cohort_consistency.RDS")
