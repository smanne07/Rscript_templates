##Create input files for GSEA,ranked gene lists or genesets using top200 genes by significance

##Create ranked gene list for a comparison from a genestats file
library(dplyr)
rdata<-read.csv("TFHTEFF_genes_stats.csv",header = TRUE)

#Change log2FoldChange and padj depending upon columns correspondingly from the stats file
rdata = rdata %>% mutate(log2FCsign = sign(log2FoldChange))
rdata = rdata %>% mutate(logP = -log10(padj))


rdata = rdata %>% mutate(rankMetric = logP/log2FCsign)
#assign the gene symbol if it is the rownames of gene stats table or any other column
rdata$geneSymbol=rownames(rdata)

#Ordering the ranked list based on metric used

rdata = rdata %>% select(geneSymbol, rankMetric) %>% arrange(desc(rankMetric))

#Write the ranked file used as input in GSEA preranked
write.table(rdata, file = "GSE58596_TFHTEFF.rnk", sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

##Create gmt geneset file for a comparison from a genestats file

#Choose file name to write the genes to a geneset
filename<-"TFHTEFF_top_bott_200_by_fdr.gmt"
#The results are already sorted if not
rdata=rdata[order(rdata$padj),]

rdata_up<-rdata[rdata$logFC>0,]

rdata_dn<-rdata[rdata$logFC<0,]

#Change the number of top significant genes accordingly
i=200

rdata_up_topsig<-rdata_up[1:i,]

rdata_dn_topsig<-rdata_dn[1:i,]

genes_up<-as.character(rdata_up_topsig$geneSymbol)

genes_dn<-as.character(rdata_dn_topsig$geneSymbol)

l1 = c("Sig_genes_up","\t",filename,"\t",paste(genes_up,collapse="\t"),"\n")

l2 = c("Sig_genes_dn","\t",filename,"\t",paste(genes_dn,collapse="\t"))

cat(l1, file=filename, append=F, sep="")
cat(l2, file=filename, append=T, sep="")
