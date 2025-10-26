library(readxl)
library(ggplot2)
library(ggpubr)
library(forcats)
library(gridExtra)
library(ggplotify)
library(grid)
vene_comb=read_excel("../DataSets/bcd-23-0014_table_s9_suppst9.xlsx",col_types='text')
rels=subset(vene_comb, ((drug_type=="single agent")& (drug!="VEN")), select="drug") 
plts=list()

rels2=subset(vene_comb, drug_type=="combo", select="drug") 


for (drug_ in unique(unlist(rels))){
 a=subset(vene_comb, drug%in%c("VEN", drug_, paste  ("VEN + ", drug_, sep="")), select=c(patient_id, drug, probit_auc))
  a$probit_auc=as.numeric(a$probit_auc)
 plts[[drug_]]=ggplot(a, aes(drug, probit_auc)) + 
   geom_boxplot()+ 
    theme_bw() +xlab("") +ylab("")  + theme(text = element_text(size=16)) 
}
num_plots <- length(plts)
ncol_value <- ceiling(sqrt(num_plots))  # Adjust based on the total plot count
plt_all=ggarrange(plotlist=plts[1:num_plots],ncol = 4, nrow = ceiling(num_plots / 4)) 
plt_all <- annotate_figure(plt_all, left = text_grob("Area Under the Curve (AUC)", rot = 90, size = 18))

#pdf("Figures/combo_altogether_aucs.pdf", width=26, height=26) 
#  grid.newpage()  
#  pushViewport(viewport(angle = 90)) 
#  print(plt_all, newpage = FALSE)  
#dev.off()
              
dir.create("ManuscriptFigures", showWarnings = FALSE)
tiff("ManuscriptFigures/Figure3.TIFF",width=12000, res = 600, height=12000) 
grid.newpage()  
#pushViewport(viewport(angle = 90)) 
print(plt_all, newpage = FALSE)  
dev.off()


#jpeg("Figures_paper/JPEG/Figure3.jpeg", width=1200, height=1200) 
#grid.newpage()  
#pushViewport(viewport(angle = 90)) 
#print(plt_all, newpage = FALSE)  
#dev.off()


