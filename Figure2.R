library(ggplot2)
library(dplyr)
library(tidyverse)
library(paletteer)
library(ggpubr)
library(tibble)
single_ktsp=readRDS("Accuracies/kTSP_single.RDS")
single_ktsp_bal=readRDS("Accuracies/kTSP_single_balanced.RDS")
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
colnames(ktsp)[colnames(ktsp)=="auc"]="AUC"
colnames(ktsp)[colnames(ktsp)=="sensitivity"]="Sensitivity" 
colnames(ktsp)[colnames(ktsp)=="specificity"]="Specificity" 
colnames(ktsp)[colnames(ktsp)=="balanced_accuracy"]="Balanced Accuracy"


ktsp=ktsp[, c("Drug","Sensitivity", "Specificity", "AUC", "Balanced Accuracy", "Method", "Balancing")]
others=rbind(single_acc_other, combo_acc_other)
others=others[,c("Drug","Sensitivity", "Specificity", "AUC", "Balanced Accuracy", "Method", "Balancing")]

alls=rbind.data.frame(others, ktsp)
alls$Sensitivity=as.numeric(alls$Sensitivity)
alls$Specificity=as.numeric(alls$Specificity)
#alls$F1=as.numeric(alls$F1)
#alls$MCC=as.numeric(alls$MCC)
alls$AUC=as.numeric(alls$AUC)
alls$`Balanced Accuracy`=as.numeric(alls$`Balanced Accuracy`)
#alls$`Pos Pred Value`=as.numeric(alls$`Pos Pred Value`)
#alls$`Neg Pred Value`=as.numeric(alls$`Neg Pred Value`)

#alls2=alls[,colnames(alls)%notin%c("F1", "MCC")]
df_long <- alls %>%
  pivot_longer(
    cols = c("Sensitivity","Specificity", "Balanced Accuracy","AUC"),   # columns to make long
    names_to = "Metric",                  # new column for metric names
    values_to = "Accuracy"                # new column for values
  )

#plt_dat$Full[plt_dat$Full=="Area Under the ROC Curve (AUC)"]="Area Under the ROC Curve"


plt1=ggplot(df_long, aes(x=Method, y=Accuracy, fill=Balancing)) +
  geom_boxplot()+
  facet_grid(~Metric) +   theme_bw()  +
  theme(axis.text.x = element_text(angle = 90, 
                  vjust = 0.5, hjust = 1)) + 
  labs(fill = "Balancing") + 
  theme(legend.position = "left") +
  scale_fill_manual(values = c("#56B4E9", "#0072B2")) +xlab("")

df_long$B=NA
df_long$B[df_long$Balancing=="Not Balanced"]="NB"
df_long$B[df_long$Balancing=="Balanced"]="B"
df_long$Merged=paste(df_long$Method, df_long$B)
p <- ggplot(df_long, aes(x= Drug, y=Merged)) +
  geom_tile(aes(fill = Accuracy)) +  
  facet_wrap(~Metric) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  
  scale_fill_viridis_c(option = "B", direction = 1) +
  labs(title = "", y = "Method") +theme(legend.position="left")+ylab("")+ 
  theme(
    axis.text.x = element_text(size = 7), # Set size to 10 points
    axis.text.y = element_text(size = 7)
  )


pp=ggarrange(p, plt1,nrow=2, heights = c(1.6,1), labels=c("A", "B"))
#ggsave(pp, file="acc.pdf", width=10, height=8)

ggsave(pp, file="ManuscriptFigures/Figure2.TIFF",  dpi=500, width=10, height=8)


   
