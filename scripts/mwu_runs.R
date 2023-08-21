###############-------------- initial MWU run------------##############

setwd("C:/Users/iXilef/OneDrive/Master/initial/MWU")
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
setwd("C:/Users/iXilef/OneDrive/Master/initial/MWU")
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

##################----------------Shuffle MWU run -------------#################

setwd("C:/Users/iXilef/OneDrive/Master/initial")

mock <- unique(cbind.data.frame("id" = BP_clust$seq, "stat" = 1))
setwd("C:/Users/iXilef/OneDrive/Master/initial/MWU")
source("gomwu.functions.R")
library(doParallel)
library(foreach)
cores=detectCores()
cl <- makeCluster(cores[1]-2)
registerDoParallel(cl)
#only once
X <- read.table("BP_X.csv", header = T)

# Edit these to match your data file names: 
input="X.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="annotations.BP.go" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# extraOptions=paste("largest=",largest," smallest=",smallest," cutHeight=",clusterCutHeight,sep="")



test <- gomwuStats_cut(input, goDatabase, goAnnotations, goDivision,
                       perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                       largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                       smallest=10,   # a GO category should contain at least this many genes to be considered
                       clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                       iteration=1,
                       #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
                       #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                       #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module
)


resultandseed <- function(result=NULL,seed=NULL)
{
  me <- list(
    result = result,
    seed = seed,
    seq = seq
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"resultandseed")
  return(me)
}

B <- 100

mwu_job <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
  resu <- gomwuStats_cut(input, goDatabase, goAnnotations, goDivision,
                         perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                         largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                         smallest=10,   # a GO category should contain at least this many genes to be considered
                         clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                         iteration=1,
                         #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
                         #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                         #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
  )
  results <- resultandseed()
  results$result <- resu[[1]]
  results$seed <- resu[[2]]
  results$seq <- resu[[3]]
  return(results)
}


## sample from norm
 
B <- 100
mwu_job_norm <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
  resu <- gomwuStats_cut(input, goDatabase, goAnnotations, goDivision,
                         perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                         largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                         smallest=10,   # a GO category should contain at least this many genes to be considered
                         clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                         iteration=3, # integer, 1: random sample from original vector, 2: random shuffle of go terms, 3: sample from N(0,var(log2fc))
                         #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
                         #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                         #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
  )
  results <- resultandseed()
  results$result <- resu[[1]]
  results$seed <- resu[[2]]
  results$seq <- resu[[3]]
  return(results)
}

save(mwu_job, mwu_job_norm, mwu_run, file="./MWU_runs.RData")


# sample from normal 2

B <- 1000
mwu_job_norm_2 <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
  resu <- gomwuStats_cut(input, goDatabase, goAnnotations, goDivision,
                         perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                         largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                         smallest=10,   # a GO category should contain at least this many genes to be considered
                         clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                         iteration=4, # integer, 1: random sample from original vector, 2: random shuffle of go terms, 3: sample from N(0,var(log2fc))
                         #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
                         #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                         #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
  )
  results <- resultandseed()
  results$result <- resu[[1]]
  results$seed <- resu[[2]]
  results$seq <- resu[[3]]
  return(results)
}

# sample from linear model

B <- 100
mwu_job_lin <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
  resu <- gomwuStats_cut(input, goDatabase, goAnnotations, goDivision,
                         perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                         largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                         smallest=10,   # a GO category should contain at least this many genes to be considered
                         clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                         iteration=5, # integer, 1: random sample from original vector, 2: random shuffle of go terms, 3: sample from N(0,var(log2fc))
                         #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
                         #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                         #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
  )
  results <- resultandseed()
  results$result <- resu[[1]]
  results$seed <- resu[[2]]
  results$seq <- resu[[3]]
  return(results)
}

# sample from loess model

B <- 100
mwu_job_loess <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
  resu <- gomwuStats_cut(input, goDatabase, goAnnotations, goDivision,
                         perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                         largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                         smallest=10,   # a GO category should contain at least this many genes to be considered
                         clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                         iteration=6, # integer, 1: random sample from original vector, 2: random shuffle of go terms, 3: sample from N(0,var(log2fc))
                         #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
                         #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                         #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
  )
  results <- resultandseed()
  results$result <- resu[[1]]
  results$seed <- resu[[2]]
  results$seq <- resu[[3]]
  return(results)
}
