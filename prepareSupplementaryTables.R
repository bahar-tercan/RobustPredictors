library(switchBox)
library(writexl)
library(tibble)
library(ggplot2)
library(ggpubr)

single_drug_models=readRDS("Outputs/singledrug_models.RDS")
combo_drug_models=readRDS("Outputs/combodrug_models.RDS")

GetRules=function(models){
  scores=list()
  drug_names=names(models)
  for (i in 1:length(models)){
    model=models[[i]]
    tm=model$TSPs
    scores[[i]]=paste(tm[,1], tm[,2], sep=">")
  }
  max_len <- max(unlist(lapply(scores, length)))
  score_list_padded <- lapply(scores, function(x){c(x, rep(NA, max_len - length(x)))})
  result <- do.call(cbind, score_list_padded)
  colnames(result)=drug_names
  return(result)
}
single_drug_predictors = GetRules(single_drug_models)
combo_drug_predictors =  GetRules(combo_drug_models)
dir.create("SupplementaryFiles")

write_xlsx(
  x = as.data.frame(single_drug_predictors),
  path = "SupplementaryFiles/singledrug_predictors.xlsx"
)
write_xlsx(
  x = as.data.frame(combo_drug_predictors),
  path = "SupplementaryFiles/combodrug_predictors.xlsx"
)

