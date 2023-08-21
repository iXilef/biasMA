setwd("C:/Users/Flix/OneDrive/Master/initial/MWU")
# table(BP_clust$term)
# X <- unique(BP_clust[,c("seq","FC.sep")]);colnames(X)<-c("id","stat")
# 
# summary(BP_clust$term)
# 
# X_assign <- X[!(X$id %in% BP_clust$seq[which(is.na(BP_clust$term))]),]
# write.table(X,file="X_noassign",sep=",",dec=".",quote=F,row.names=F,col.names=T)
# write.table(X,file="X_assign",sep=",",dec=".",quote=F,row.names=F,col.names=T)


## for df with NA annots

input="X.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="annotations.BP.go" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
setwd("C:/Users/Flix/OneDrive/Master/initial/MWU")
source("gomwu.functions.R")

# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

mwu_run <- gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=10,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# with_NA_shuffle <- gomwuStats(input, goDatabase, goAnnotations, goDivision,
#                       perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
#                       largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
#                       smallest=10,   # a GO category should contain at least this many genes to be considered
#                       clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#                       adjust.multcomp="shuffle",
#                       #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#                       #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#                       #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
# )
# 
# ## for df with no NA annots
# 
# input="X_assign" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
# goAnnotations="annotations.BP.go" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
# goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
# goDivision="BP" # either MF, or BP, or CC
# setwd("C:/Users/Flix/OneDrive/Master/initial/MWU")
# source("gomwu.functions.R")
# 
# # ------------- Calculating stats
# # It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
# 
# without_NA <- gomwuStats(input, goDatabase, goAnnotations, goDivision,
#            perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
#            largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
#            smallest=10,   # a GO category should contain at least this many genes to be considered
#            clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#            #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#            #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#            #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
# )
# # do not continue if the printout shows that no GO terms pass 10% FDR.
# 
# without_NA_shuffle <- gomwuStats(input, goDatabase, goAnnotations, goDivision,
#                          perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
#                          largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
#                          smallest=10,   # a GO category should contain at least this many genes to be considered
#                          clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#                          adjust.multcomp="shuffle",
#                          #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#                          #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#                          #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
# )
# # do not continue if the printout shows that no GO terms pass 10% FDR.

# for plots
framie <- without_NA[,c("term", "pval", "p.adj")]
framie$pval_2 <-  without_NA_shuffle$pval[match(without_NA$term, without_NA_shuffle$term)]
framie$pval.adj_2 <-  without_NA_shuffle$p.adj[match(without_NA$term, without_NA_shuffle$term)]

save.image("NA_shuffles.RData")
