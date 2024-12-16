# Last edit by Josephine Gile 20191211

library(tidyverse)
library(Gviz)
library(GenomicRanges)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)

base_dir1 <- "/Users/josephine/Dropbox\ (Personal)/SCIENCE/WHERRY\ LAB/shared_folders/SM_JG/KP_OK_jogiles_processed/bws/merged_bws/"
#base_dir2 <- "/Volumes/EJW_Giles_13/ATAC_L007_L008_L009/processedFiles/merged_deblBams/merged_bws/"
wd <- "/Users/josephine/Dropbox\ (Personal)/SCIENCE/WHERRY\ LAB/shared_folders/SM_JG/KP_OK_jogiles_processed/bws/gviz/"

fileAdd <- ""

#################

# Coordinates for genes of interest
oi1 <- c("chr4",6659387, 7005459, 0,7.5,"Tox")

# Generate list of all genes -- have to add element for each
chr_oi <- c(oi1[1])
start_oi <- c(oi1[2])
end_oi <- c(oi1[3])
ymin_oi <- c(oi1[4])
ymax_oi <- c(oi1[5])
geneName_oi <- c(oi1[6])


# Read in data for choosing peaks to highlight
#results <- read_csv(file = paste0(wd, "CXCR3_poi.csv"), col_names = TRUE)


#################

# Set Schemes

getOption("Gviz.scheme")

scheme <- getScheme()
scheme$GeneRegionTrack$thinBoxFeature<-TRUE
scheme$GeneRegionTrack$transcriptAnnotation <- "symbol"
addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")

gen<-"mm10"

#################
## sample color scheme
# Naive #33a02c, Arm D8 #a6cee3, Arm D15 #5289C7, Arm D30 #1740B6, Cl13 D8 #F7EE55, Cl13 D15 #F1932D,Cl13 D30 #DC050C, Cl13 D30+aPDL1, #901DB3
#sample_col <- c("#33a02c", "#a6cee3", "#5289C7", "#1740B6", "#FFD700", "#F1932D", "#DC050C")
sample_col <- c("#33a02c", "#a6cee3", "#1740B6", "#DC050C", "#901DB3")
#################

# BW files to be used
bwFile1 <- paste0(base_dir1, "naive_merged_sort.bw")
bwFile2 <- paste0(base_dir1, "effector_merged_sort.bw")
bwFile3 <- paste0(base_dir1, "memory_merged_sort.bw")
bwFile4 <- paste0(base_dir1, "exhausted_merged_sort.bw")
bwFile5 <- paste0(base_dir1, "aPDL1_merged_sort.bw")



# Specify colors for each track
col1 <- sample_col[1]
col2 <- sample_col[2]
col3 <- sample_col[3]
col4 <- sample_col[4]
col5 <- sample_col[5]

# coord list -- (can make vectors to combine for multiple genes)
coord_list <- list(chrs = chr_oi,strts = start_oi,ends = end_oi,ymins = ymin_oi, ymaxs = ymax_oi,fnames = geneName_oi)



#######################################

i <- 1


# Create data tracks for each of our bigwig files
cexax <- 0.7 # axis text size
cexti <- 0.6 # title text size


chr <- coord_list[[1]][[i]]
start_pos <- as.numeric(coord_list[[2]][[i]])
end_pos <- as.numeric(coord_list[[3]][[i]])
y_min <- as.numeric(coord_list[[4]][[i]])
y_lim <- as.numeric(coord_list[[5]][[i]])
fname <- coord_list[[6]][[i]]



dTrack1 <- DataTrack(range = bwFile1, genome = gen,type = "histogram", chromosome = chr, name = "Naive",fill.histogram = col1,col.histogram = col1,ylim=c(y_min,y_lim),fontsize=20)
dTrack2 <- DataTrack(range = bwFile2, genome = gen,type = "histogram", chromosome = chr, name = "Effector",fill.histogram = col2,col.histogram = col2,ylim=c(y_min,y_lim),fontsize=20)
dTrack3 <- DataTrack(range = bwFile3, genome = gen,type = "histogram", chromosome = chr, name = "Memory",fill.histogram = col3,col.histogram = col3,ylim=c(y_min,y_lim),fontsize=20)
dTrack4 <- DataTrack(range = bwFile4, genome = gen,type = "histogram", chromosome = chr, name = "Exhausted",fill.histogram = col4,col.histogram = col4,ylim=c(y_min,y_lim),fontsize=20)
dTrack5 <- DataTrack(range = bwFile5, genome = gen,type = "histogram", chromosome = chr, name = "Exhausted+aPDL1",fill.histogram = col5,col.histogram = col5,ylim=c(y_min,y_lim),fontsize=20)

# set display parameters for biomTrack for dTracks (histograms)
# Plot the tracks changing the from and to for a particular location and include the type of tracks needed.
displayPars(dTrack1) <- list(baseline=0,lty.baseline=1,col.baseline=col1,col.grid="black",col.axis="black",cex.axis=cexax,cex.title=cexti,fontcolor.title="black",col.border.title="black",col.frame="white",background.title="white")
displayPars(dTrack2) <- list(baseline=0,lty.baseline=1,col.baseline=col2,col.grid="black",col.axis="black",cex.axis=cexax,cex.title=cexti,fontcolor.title="black",col.border.title="black",col.frame="white",background.title="white")
displayPars(dTrack3) <- list(baseline=0,lty.baseline=1,col.baseline=col3,col.grid="black",col.axis="black",cex.axis=cexax,cex.title=cexti,fontcolor.title="black",col.border.title="black",col.frame="white",background.title="white")
displayPars(dTrack4) <- list(baseline=0,lty.baseline=1,col.baseline=col4,col.grid="black",col.axis="black",cex.axis=cexax,cex.title=cexti,fontcolor.title="black",col.border.title="black",col.frame="white",background.title="white")
displayPars(dTrack5) <- list(baseline=0,lty.baseline=1,col.baseline=col5,col.grid="black",col.axis="black",cex.axis=cexax,cex.title=cexti,fontcolor.title="black",col.border.title="black",col.frame="white",background.title="white")


# Load genomic axis track and the sequence track if needed
axisTrack <- Gviz::GenomeAxisTrack(scale=TRUE)
sTrack <- SequenceTrack(Mmusculus)

#ensembl <- useMart("ensembl")
#datasets <- listDatasets(ensembl)

bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

biomTrack <- BiomartGeneRegionTrack(genome = gen,chromosome = chr, start = as.numeric(start_pos), end = as.numeric(end_pos) ,name = "ENSEMBL",cex=0.9,fontsize=12,cex.group=1,stacking="squish", biomart = bm)





#biomTrack <- BiomartGeneRegionTrack(genome = gen,chromosome = chr, start = as.numeric(start_pos), end = as.numeric(end_pos) ,name = "ENSEMBL",filter = list(with_ox_refseq_mrna = TRUE),cex=0.9,fontsize=12,cex.group=1,stacking="squish")
#biomTrack <- BiomartGeneRegionTrack(genome = gen,chromosome = chr, start = as.numeric(start_pos), end = as.numeric(end_pos) ,name = "ENSEMBL", biomart=useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", version = "GRCm38.p6", host = "http://useast.ensembl.org"), filter = list(with_ox_refseq_mrna = TRUE),cex=0.9,fontsize=12,cex.group=1,stacking="squish")




# set display parameters for biomTrack
displayPars(biomTrack) <-list(collapseTranscripts=FALSE,fill="black",col.grid="black",C_segment="black",D_segment="black",J_segment="black",Mt_rRNA="black",Mt_tRNA="black",Mt_tRNA_pseudogene="black",V_segment="black",miRNA="black",miRNA_pseudogene="black",misc_RNA="black",misc_RNA_pseudogene="black",protein_coding="black",pseudogene="black",rRNA="black",rRNA_pseudogene="black",retrotransposed="black",scRNA="black",scRNA_pseudogene="black",snRNA="black",snRNA_pseudogene="black",snoRNA="black",snoRNA_pseudogene="black",tRNA_pseudogene="black",utr3="black",utr5="black",geneSymbols=TRUE,col.line="black",col="black",cexfontcolor.title="black",col.border.title="black",col.frame="white",background.title="white")

# set display parameters for axisTrack
displayPars(axisTrack)<-list(col.frame="black",col.grid="black",col="black",fill="black",col.line="black",fontcolor="black",col.axis="black",fontcolor.title="black")






# Select significant peaks from results to highlight
#oi_peaks <- results # %>% filter(geneName == fname) %>% filter(start > start_pos) %>% filter(end < end_pos)

#dTracks_plusBoxes <- HighlightTrack(trackList = list(dTrack1,dTrack2,dTrack3,dTrack4,dTrack5,dTrack6,dTrack7,dTrack8,dTrack9,dTrack10,dTrack11), start = oi_peaks$start, end = oi_peaks$end, chromosome = chr)

#displayPars(dTracks_plusBoxes) <- list(col = "#DCDCDC",fill = "#DCDCDC")




png(filename = paste0(wd,fileAdd,fname,"_noBoxes.png"),width=6000,height=3000,units="px",res = 300)
plotTracks(list(axisTrack,sTrack,dTrack1,dTrack2,dTrack3,dTrack4,dTrack5, biomTrack),from = as.numeric(start_pos), to = as.numeric(end_pos) ,chromosome = chr, scale = 0.5)
dev.off()

#png(filename = paste0(wd,fileAdd,fname,"_wBoxes.png"),width=6000,height=3000,units="px",res = 300)
#plotTracks(list(axisTrack,sTrack,dTracks_plusBoxes,biomTrack),from = as.numeric(start_pos), to = as.numeric(end_pos) ,chromosome = chr, scale = 0.5)
#dev.off()


setEPS()
postscript(paste0(wd,fileAdd,fname,"_noBoxes.eps"),family="sans",colormodel = "rgb",width = 20,height=18,pointsize = 24,paper = "special")
plotTracks(list(axisTrack,sTrack,dTrack1,dTrack2,dTrack3,dTrack4,dTrack5,biomTrack),from = as.numeric(start_pos), to = as.numeric(end_pos) ,chromosome = chr, scale = 0.5)
dev.off()

#setEPS()
#postscript(paste0(wd,fileAdd,fname,"_wBoxes.eps"),family="sans",colormodel = "rgb",width = 20,height=18,pointsize = 24,paper = "special")
#plotTracks(list(axisTrack,sTrack,dTracks_plusBoxes,biomTrack),from = as.numeric(start_pos), to = as.numeric(end_pos) ,chromosome = chr, scale = 0.5)
#dev.off()