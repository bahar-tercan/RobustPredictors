library(ggplot2)
library(dplyr)
library(tidyverse)
library(paletteer)
library(ggpubr)
library(tibble)
single_ktsp=readRDS("Accuracies/kTSP_single.RDS")
single_ktsp_bal=readRDS("Accuracies2/kTSP_single_balanced.RDS")
single_ktsp=rownames_to_column(as.data.frame(single_ktsp), var="Drug")
single_ktsp_bal=rownames_to_column(as.data.frame(single_ktsp_bal), var="Drug")

combo_ktsp=readRDS("Accuracies/kTSP_not_balanced_combo.RDS")
combo_ktsp_bal=readRDS("Accuracies/kTSP_balanced_combo.RDS")
combo_ktsp=rownames_to_column(as.data.frame(combo_ktsp), var="Drug")
combo_ktsp_bal=rownames_to_column(as.data.frame(combo_ktsp_bal), var="Drug")

single_ktsp$Method="kTSP"
single_ktsp_bal$Method="kTSP"
combo_ktsp$Method="kTSP"
combo_ktsp_bal$Method="kTSP"

single_ktsp$Balancing="Not Balanced"
single_ktsp_bal$Balancing="Balanced"
combo_ktsp$Balancing="Not Balanced"
combo_ktsp_bal$Balancing="Balanced"

single_ktsp$Type="Single"
single_ktsp_bal$Type="Single"
combo_ktsp$Type="Combination"
combo_ktsp_bal$Type="Combination"

`%notin%`=Negate(`%in%`)
combo_acc_other=readRDS("Accuracies/combo_other_classifiers.RDS")
single_acc_other=readRDS("Accuracies/SingleDrug_otheracc.RDS")
combo_acc_other$Type="Combination"
single_acc_other$Type="Single"

ktsp=rbind(single_ktsp, single_ktsp_bal,combo_ktsp, combo_ktsp_bal)
ktsp=ktsp[,colnames(ktsp)!="accuracy"]
colnames(ktsp)[colnames(ktsp)=="auc"]="AUROC"
colnames(ktsp)[colnames(ktsp)=="sensitivity"]="Sensitivity" 
colnames(ktsp)[colnames(ktsp)=="specificity"]="Specificity" 
colnames(ktsp)[colnames(ktsp)=="balanced_accuracy"]="Balanced Accuracy"


ktsp=ktsp[, c("Drug","Sensitivity", "Specificity", "AUROC", "Balanced Accuracy", "Method", "Balancing")]
others=rbind(single_acc_other, combo_acc_other)
colnames(others)[colnames(others)=="AUC"]="AUROC"
others=others[,c("Drug","Sensitivity", "Specificity", "AUROC", "Balanced Accuracy", "Method", "Balancing")]

alls=rbind.data.frame(others, ktsp)
alls$Sensitivity=as.numeric(alls$Sensitivity)
alls$Specificity=as.numeric(alls$Specificity)
alls$AUROC=as.numeric(alls$AUROC)
alls$`Balanced Accuracy`=as.numeric(alls$`Balanced Accuracy`)

df_long <- alls %>%
  pivot_longer(
    cols = c("Sensitivity","Specificity", "Balanced Accuracy","AUROC"),   # columns to make long
    names_to = "Metric",                  # new column for metric names
    values_to = "Accuracy"                # new column for values
  )


plt1=ggplot(df_long, aes(x=Method, y=Accuracy, fill=Balancing)) +
  geom_boxplot()+
  facet_grid(~Metric) +   theme_bw()  +
  theme(axis.text.x = element_text(angle = 90, 
                  vjust = 0.5, hjust = 1)) + 
  labs(fill = "Balancing") + 
  theme(legend.position = "left") +
  scale_fill_manual(values = c("#56B4E9", "#0072B2")) + xlab("Classifier")

df_long$B=NA
df_long$B[df_long$Balancing=="Not Balanced"]="NB"
df_long$B[df_long$Balancing=="Balanced"]="B"
df_long$Merged=paste(df_long$Method, df_long$B)

p <- ggplot(df_long, aes(x = Drug, y = Merged)) +
  geom_tile(aes(fill = Accuracy), color = "black") +  
  facet_wrap(~Metric) + 
  theme_bw() +
  scale_fill_viridis_c(option = "B", direction = 1) + 
  labs(y = "Classifier and Balancing Status") + 
  theme(
    legend.position = "left",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 10) 
    
  )

pp = ggarrange(
  p, 
  NULL, 
  plt1, 
  nrow = 3, 
  heights = c(1.6, 0.1, 1), 
  labels = c("A", "", "B")   
) +
theme(plot.background = element_rect(fill = "white", color = NA), base_size = 10, base_family = "Arial") 
ggsave(
  filename = "Plots/Fig2.pdf",
  plot = pp,
 # device = cairo_ps,  
  width = 7.5,         
  height = 7,
  units = "in",
  dpi = 300            
)

   
