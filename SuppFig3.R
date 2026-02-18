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
  inc_samples_denovo=unique(unlist(subset(xx, isDenovo=="TRUE", select="Sample")))
  inc_samples_notdenovo=unique(unlist(subset(xx, isDenovo=="FALSE",select="Sample")))
  denovo_labels=unlist(xx[match(inc_samples_denovo, xx$Sample), "SensitivityCall"])
  notdenovo_labels=unlist(xx[match(inc_samples_notdenovo, xx$Sample), "SensitivityCall"])
  if (length(table(denovo_labels))==2 & all(table(denovo_labels)>=20)){ 
    if (length(table(notdenovo_labels))==2 & all(table(notdenovo_labels)>=20)){ 
      i=i+1
      drugs[[i]]=drug1 
    }
  }
}
merged_drug=subset(merged, Drug%in%unlist(drugs))

merged_drug$SenC=ifelse(merged_drug$SensitivityCall=="Sensitive", "S", "R")


# --- Clean the Drug Grouping ---
merged_drug <- merged_drug %>%
  # Ensure SensitivityCall is a factor for consistent plotting order
  mutate(SensitivityCall = factor(SensitivityCall, 
                                  levels = c("Resistant", "Sensitive")))

# --- Improved Plot ---
p_denovo <- ggplot(merged_drug, 
                   aes(x = SensitivityCall, 
                       y = Stat, 
                       fill = isDenovo)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Drug, scales = "free_y") + 
  theme_bw() +
  labs(x = "In-vitro Response", 
       y = "k-TSP Score (Sum of Votes)", 
       fill = "Is De Novo?") +
  scale_fill_manual(values=palette.colors(palette = "Set2")[1:2], 
                    
                    labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(
      size = 6,             
      margin = ggplot2::margin(1, 0, 1, 0), 
      face = "bold"
    ),    axis.text.x = element_text(angle = 45, hjust = 1)
  )+theme(legend.position="left") 

saveRDS(p_denovo, file="PlotData/denovo_consistency.RDS")

p3=readRDS("PlotData/cohort_consistency.RDS")
p3=p3 +theme(legend.position="left") +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1 , size = 8))  

p_denovo <- p_denovo + theme(
  axis.text.y = element_text(size = 6),      
  axis.title.y = element_text(size= 10),     
  
  strip.text.x = element_text(size = 5,      # Even smaller drug names
                              margin = ggplot2::margin(t = 2, b = 2)), 
  
  axis.text.x = element_text(size = 6),     
  
  panel.spacing = unit(0.05, "lines"),
  plot.margin =ggplot2::margin(5, 5, 5, 5)          
)

p3 <- p3 + theme(
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size = 10),
  strip.text.x = element_text(size = 5, margin = ggplot2::margin(t = 2, b = 2)),
  axis.text.x = element_text(size = 6),
  panel.spacing = unit(0.05, "lines")
)

p_all = ggarrange(
  p3, 
  NULL, 
  p_denovo, 
  nrow = 3, 
  heights = c(1, 0.1, 1), 
  labels = c("A", "", "B")   
) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(p_all, file="../Plots/Fig_S3.pdf", dpi=600, width=10)
