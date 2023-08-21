setwd("C:/Users/Flix/OneDrive/Master/initial/MWU")
library(doParallel)
library(foreach)
cores=detectCores()
cl <- makeCluster(cores[1]-2)
registerDoParallel(cl)
# Edit these to match your data file names: 
input="all_topgo.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="all_annots.csv" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")

# run it once to generate the dissim and clustering files

# test <- gomwuStats(input, goDatabase, goAnnotations, goDivision,
#                        perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
#                        largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
#                        smallest=10,   # a GO category should contain at least this many genes to be considered
#                        clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
#                        #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
#                        #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#                        #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module
# )


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

## topgo data run


B <- 100
mwu_job_lin_all <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
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
mwu_job_loess_all <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
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

################ MWU testset


input="heats.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="amil_defog_iso2go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC



test <- gomwuStats(input, goDatabase, goAnnotations, goDivision,
                   perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
                   largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
                   smallest=10,   # a GO category should contain at least this many genes to be considered
                   clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.,
                   #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
                   #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
                   #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module
)



B <- 100
mwu_job_lin_heat <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
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
mwu_job_loess_heat <- foreach(1:B, .export="gomwuStats_cut") %dopar% {
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

mwu_add <- list(heat= list(linear=mwu_job_lin_heat, loess=mwu_job_loess_heat), all= list(linear=mwu_job_lin_all, loess=mwu_job_loess_all))
save(mwu_add, file="add_data_algo4.RData")
