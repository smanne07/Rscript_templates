setwd("~/Desktop/sasi_work/MA_RNA_Seq/Genesets_various")
source('~/Desktop/sasi_work/test_R_webapps/Multi_comparison/functions_SM_wlab.R')

gmtdata<-read.delim("20181008_ThelperRUNXGeneSets_FromRamin.gmx",header = F)


Th1_sig=as.character(gmtdata[3:202,3])
Th2_sig=as.character(gmtdata[3:202,5])
Th17_sig=as.character(gmtdata[3:202,1])

geneOverlap_circos3(Gene_list = list(Th1_sig,Th2_sig,Th17_sig))

library(Vennerable)

venn_list<-list(Th1_sig=Th1_sig,Th2_sig=Th2_sig,Th17_sig=Th17_sig)

Venn_obj <- Venn(venn_list)


plovenn=compute.Venn(Venn_obj,doWeights = FALSE)

gpList=VennThemes(plovenn)

a=plot(Venn_obj,show=list(Faces=FALSE),gpList=gpList)


Th1_exc_sig=Venn_obj@IntersectionSets$`100`
Th2_exc_sig=Venn_obj@IntersectionSets$`010`
Th17_exc_sig=Venn_obj@IntersectionSets$`001`

#Create gmt

filename<-"GSE14308_Th1_Th2_Th17_excl_sigs.gmt"

l1 = c("Th1_exc_sig","\t",filename,"\t",paste(Th1_exc_sig,collapse="\t"),"\n")

l2 = c("Th2_exc_sig","\t",filename,"\t",paste(Th2_exc_sig,collapse="\t"),"\n")

l3 = c("Th17_exc_sig","\t",filename,"\t",paste(Th17_exc_sig,collapse="\t"))

cat(l1, file=filename, append=F, sep="")
cat(l2, file=filename, append=T, sep="")
cat(l3, file=filename, append=T, sep="")
