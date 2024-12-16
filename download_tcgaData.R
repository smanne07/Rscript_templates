library(TCGA2STAT)
library(dplyr)
dt_key=Sys.Date()
download_tcga=function(tcga_key){
  tcga_rna_clinical <- getTCGA(disease=tcga_key, data.type="RNASeq2", type="RPKM", clinical=TRUE)
  
  saveRDS(tcga_rna_clinical,paste0('tcga_rna_clinical_',tcga_key,"_",dt_key,'.rds'))
  
}  
#download_tcga(tcga_key="")

lapply(list("COAD","SKCM","LUAD","LUSC"),download_tcga)
