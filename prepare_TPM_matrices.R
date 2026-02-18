# convert gene expression matrices to TPM
library(readr)
library(rtracklayer)
library(tibble)
library(DGEobj.utils)

beataml_gexp=read_tsv("../DataSets/beataml_waves1to4_norm_exp_dbgap.txt")
beataml_gexp=subset(beataml_gexp, biotype=="protein_coding")
beataml_gexp2=beataml_gexp[,5:ncol(beataml_gexp)]
beataml_gexp2=2^beataml_gexp2
beataml_tpm=apply(beataml_gexp2, 2, function(x){(10^6)*x/sum(x)})
beataml_tpm=log2(beataml_tpm)
rownames(beataml_tpm)=beataml_gexp$stable_id
beataml_anno=beataml_gexp[,1:4]
dir.create("GeneExpMatrices")
write.csv(beataml_tpm, "GeneExpMatrices/beataml_logged_tpm.csv")
write.csv(beataml_anno, "GeneExpMatrices/beataml_gene_anno.csv")

fimm_count=read.csv("../DataSets/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv")

gtf_file="GeneExpMatrices/gencode.v36.annotation.gff3.gz"
gc_obj <- rtracklayer::import(con = gtf_file)
gc_obj$length <- rtracklayer::width(gc_obj)
geneAnnotation <- as.data.frame(gc_obj)
geneAnnotation <- geneAnnotation[geneAnnotation$type=="gene",]
geneAnnotation$gene_id <- geneAnnotation$ID
geneAnnotation <- geneAnnotation[,c("gene_id", "gene_name", "hgnc_id", "gene_type", "length")]
gene_ann=subset(geneAnnotation, gene_type=="protein_coding")
gene_ann$Ensembl=sapply(gene_ann$gene_id, function(x){strsplit(x, "\\.")[[1]][1]})
inc_genes=intersect(gene_ann$Ensembl,fimm_count$X)
fimm_count=column_to_rownames(fimm_count, var="X")
fimm_count=fimm_count[inc_genes, ]
rel_ann=gene_ann[match(inc_genes, gene_ann$Ensembl),]
tpm_fimm=convertCounts(as.matrix(fimm_count),"TPM", rel_ann$length)
write.csv(log2(tpm_fimm+1), "GeneExpMatrices/fimm_logged_tpm.csv")

