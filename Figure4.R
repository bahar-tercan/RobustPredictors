#PCA plot of the three cohorts
library(readxl)
library(readr)
library(tibble)
library(ggplot2)

`%notin%`=Negate(`%in%`)

beataml_exp=as.data.frame(read.csv("GeneExpMatrices/beataml_logged_tpm.csv"))
beataml_exp=column_to_rownames(beataml_exp, var="X")

cli_data=as.data.frame(read_excel("../DataSets/beataml_wv1to4_clinical.xlsx"))
rel_cli=subset(cli_data, manuscript_rnaseq=="yes"& 
       dxAtSpecimenAcquisition%in%c("ACUTE MYELOID LEUKAEMIA (AML) AND RELATED PRECURSOR NEOPLASMS"))

for_normal_samples=read_excel("../DataSets/beataml_waves1to4_sample_mapping.xlsx")
healthy_samples=for_normal_samples[sapply(for_normal_samples$rna_control, startsWith, "Healthy"),]
healthy_samples=unlist(healthy_samples[!is.na(healthy_samples$rna_control),"dbgap_rnaseq_sample"])

rel_cli2=subset(rel_cli, rel_cli$dbgap_rnaseq_sample%notin%healthy_samples)

waves1_2_samples= unlist(subset(rel_cli2, cohort=="Waves1+2", select="dbgap_rnaseq_sample"))
waves1_2=waves1_2_samples[!is.na(waves1_2_samples)]
waves3_4_samples= unlist(subset(rel_cli2, cohort=="Waves3+4", select="dbgap_rnaseq_sample"))
waves3_4=waves3_4_samples[!is.na(waves3_4_samples)]

dat=read.csv("GeneExpMatrices/fimm_logged_tpm.csv")
dat2=column_to_rownames(dat, var="X")
aml_dat= dat2[,sapply(colnames(dat2), startsWith, "AML")]
common_genes=intersect(rownames(aml_dat), rownames(beataml_exp))
aml_fimm=aml_dat[common_genes,]
rel_beat=beataml_exp[common_genes,]

#pca plot

dt1=cbind(rel_beat[,waves1_2], rel_beat[,waves3_4], aml_fimm)
p3 <- prcomp(t(dt1), retx=TRUE, center=TRUE, scale=TRUE)
summary(p3)
p3scores <- p3$x

plt_dat=cbind.data.frame( PC1=p3scores[,1], PC2=p3scores[,2], 
                          Cohort=rep.int(c("Beat AML Waves 1+2", "Beat AML Waves 3+4", "FPMTB"), 
                                         c(length(waves1_2) , length(waves3_4), ncol(aml_fimm))))


plt=ggplot(plt_dat, aes(x= PC1 , y=PC2, color=Cohort)) + geom_point() + theme_bw() +
  theme(legend.position = c(0.3, 0.08),
        legend.justification = c(0, 0),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),   
        legend.box.background = element_rect()) +
    scale_color_manual(values=c("#FFD92F", "#E5C494", "#B3B3B3"))

dir.create("ManuscriptFigures")
ggsave(plt, file="ManuscriptFigures/Figure4.TIFF", dpi=500, width=5, height=4)



