#Choose file name to write the genes to a gmt file
filename<-"genesetFilename.gmt"

##genes_list1 and genes_list2 are character vectors with gene names
#genes_list1_name and genes_list2_name are names of the genes list

#following creates for example the first two lines of the gmt file
l1 = c("genes_list1_name","\t",filename,"\t",paste(genes_list1,collapse="\t"),"\n")

l2 = c("genes_list2_name","\t",filename,"\t",paste(genes_list2,collapse="\t"))

##Add the two lines to the gmt file given by the filename
cat(l1, file=filename, append=F, sep="")
cat(l2, file=filename, append=T, sep="")
