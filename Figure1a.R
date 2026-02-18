library(dplyr)
library(tidyr)
library(ggplot2)
cli_data=readRDS("Outputs/RelevantClinical.RDS")
drug_data=readRDS("Outputs/RelevantDrugRes.RDS")
dat <- read.csv("GeneExpMatrices/beataml_logged_tpm.csv", row.names = 1)
drug_data$SensitivityCall <- factor(drug_data$SensitivityCall)
results_list <- list()
drugs=unique(drug_data$inhibitor)
waves1_2=unlist(cli_data[cli_data$cohort=="Waves1+2","dbgap_rnaseq_sample"]);names(waves1_2)=NULL
waves3_4=unlist(cli_data[cli_data$cohort=="Waves3+4","dbgap_rnaseq_sample"]);names(waves3_4)=NULL
waves1_2=waves1_2[!is.na(waves1_2)]
waves3_4=waves3_4[!is.na(waves3_4)]
for (drug in drugs){  
  drug_status <- subset(drug_data, inhibitor == drug)
  
  inc_w12 <- intersect(waves1_2, intersect(drug_status$dbgap_rnaseq_sample, colnames(dat)))
  inc_w34 <- intersect(waves3_4, intersect(drug_status$dbgap_rnaseq_sample, colnames(dat)))
  
  w12_labels <- drug_status$SensitivityCall[match(inc_w12, drug_status$dbgap_rnaseq_sample)]
  w34_labels <- drug_status$SensitivityCall[match(inc_w34, drug_status$dbgap_rnaseq_sample)]
  
  if (length(unique(w12_labels)) == 2 && all(table(w12_labels) >= 20) &&
      length(unique(w34_labels)) == 2 && all(table(w34_labels) >= 10)) {
    
     results_list[[drug]] <- list(
      train = table(w12_labels),
      test = table(w34_labels)
    )
  }
}

all_calls <- bind_rows(lapply(names(results_list), function(d) {
  data.frame(Drug = d, 
             Cohort = "Waves 1+2", 
             as.data.frame(results_list[[d]]$train)) %>%
    rename(Sensitivity = w12_labels, Count = Freq) %>%
    bind_rows(
      data.frame(Drug = d, 
                 Cohort = "Waves 3+4", 
                 as.data.frame(results_list[[d]]$test)) %>%
        rename(Sensitivity = w34_labels, Count = Freq)
    )
}))

colnames(all_calls)=c( "Drug","Cohort", "Sensitivity Call", "NumberSamples")
Figure1a <- ggplot(all_calls, aes(x = Drug, y = NumberSamples, fill = `Sensitivity Call`)) + 
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



