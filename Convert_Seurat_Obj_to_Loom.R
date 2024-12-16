#module load hdf5/1.8.21

library(loomR)
library(Seurat)
scRNA_merged <- readRDS(file = "scRNAmerged_SeuratObj_res1.2_renamed.rds")
meta.data=scRNA_merged@meta.data

##Seurat object with NA values in metadata causes errors in loom file creation
str(meta.data)

##Check possible NA values
#table(meta.data$infection,exclude=F)

 #None   Arm  Cl13 
 #6709 11454 15729 
#table(meta.data$timepoint,exclude=F)

 #  d8   d15   d30  <NA> 
 #8713  7980 10490  6709 
#table(meta.data$treatment,exclude=F)

#aPDL1  <NA> 
# 4319 29573 

#Replace the NA values

meta.data$timepoint=as.character(meta.data$timepoint)
meta.data[is.na(meta.data$timepoint),]$timepoint="d0"

meta.data$treatment=as.character(meta.data$treatment)
meta.data[is.na(meta.data$treatment),]$treatment="None"

scRNA_merged@meta.data=meta.data
#Got the error when reading loom file in scanpy
#https://github.com/theislab/scanpy/issues/598
#ValueError: column index exceeds matrix dimensions
scRNA_merged@graphs <- list()
loom_obj_sct <- as.loom(x = scRNA_merged, assay = "SCT", filename = "Armonly_scRNAmerged_SeuratObj_res1.2_renamed.loom", verbose = TRUE)