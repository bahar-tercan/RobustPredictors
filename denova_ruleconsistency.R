library(tidyr)
library(dplyr)
library(ggplot2) 
library(ggpubr)  

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
merged=subset(merged, isDenovo%in%c("TRUE", "FALSE"))
for (drug1 in ddd){ 
  xx=subset(merged, Drug==drug1)
  inc_samples_wave12=unique(unlist(subset(xx, isDenovo=="TRUE", select="Sample")))
  inc_samples_wave34=unique(unlist(subset(xx, isDenovo=="FALSE",select="Sample")))
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
merged_combined <- merged_drug %>%
  unite(col = "Drug_Sensitivity_Group", 
        c("Drug", "SenC"), 
        sep = "-", 
        remove = FALSE) # Keep the original columns

p2=ggplot(merged_combined, 
       aes(x = Drug_Sensitivity_Group, 
           y = Stat,
           fill = isDenovo)) + 
  
  geom_boxplot() +theme_bw()+ xlab("") +
  theme(legend.position="left") +
  scale_fill_manual(values=palette.colors(palette = "Set2")[1:2])+
  ylab("Sum of Votes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 , size = 8)) 
  # 1. Remove major and minor lines for both axes


p3=readRDS("PlotData/cohort_consistency.RDS")
p3=p3 +theme(legend.position="left") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1 , size = 8)) 
p_all <- ggarrange(
  p3, p2, 
  labels = c("A", "B"), # <-- This is the correct argument in ggarrange
  nrow = 2              # Arranges p3 on top (A) and p2 on bottom (B)
)
#p_all=grid.arrange(p3,p2,nrow=2, labels=c("A", "B"))
#ggsave(p_all, file="rule_consistency.pdf", width=10)
dir.create("ManuscriptFigures", showWarnings = FALSE)

ggsave(p_all, file="ManuscriptFigures/rule_consistency.TIFF", dpi=500, width=10)
