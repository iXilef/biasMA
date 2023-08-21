# make additional GO with more datasets
## preparation from topGO_runs.R

setwd("C:/Users/iXilef/OneDrive/Master/initial")
load("subset.RData")
annotations <- read.table("./MWU_clone/annotations.BP.go")
colnames(annotations) <- c("id", "go")

library(topGO)
library(foreach)
library(dplyr)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-2)
registerDoParallel(cl)
# library(ALL)
# load models
load("./MWU/models.RData")

X <- BP_clust[,c("seq","name","term","lev","FC.sep","p.sep")]
colnames(X) <- c("transcript","GOname","GOterm","GOlevel","log2FC","minuslog10pval")

# Lower case x is just a line for each transcript. 
x <- X[!duplicated(X$transcript),]
# x <- x[,c(1,5,6)]

# make for topGO an named list of character; name: gene, character vector: GO terms

annots <- setNames(annotations$go, annotations$id)
for (i in 1:length(annots)){
  annots[i] <- strsplit(annots[[i]], ";", fixed=T)
  # print(i)
}

annotsrev <- inverseList(annots)

# prune list of unknown mappings
length(which(annots=="unknown"))

annots <- annots[annots!="unknown"]

# clean x from NA pvals

x <- x[!is.na(x$minuslog10pval),]


# custom class for foreach output

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

x <- x[which(x$transcript %in% names(annots)),]
x$n <- sapply(x$transcript, function(x) length(annots[[x]]))
# runs
B <- 100

runall <- function(godata = GOdat, algor, test){
  alltests <- list()
  for(algo in algor){
    for(tests in test){
      alltests <- append(alltests, list(runTest(GOdat, algorithm = algo, statistic = tests), 
                                        name=paste(algo, "-", tests)))
    }
  }
  return(alltests)
}

ks_job_linfun <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  x$log2FC <- rnorm(length(x$log2FC), mean=predict(real_lm, data.frame(nmemberlog10=log10(x$n))),
                    sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(x$n))))
  
  geneList <- setNames(sample(x$log2FC), x$transcript)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}


# prune x for loess run

x <- x[which(x$n < 656),]

ks_job_loess <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  x$log2FC <- rnorm(length(x$log2FC), mean=predict(lo, log10(x$n)),
                    sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(x$n))))
  
  geneList <- setNames(sample(x$log2FC), x$transcript)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}

# get two other data sets ready
# prepare all dataset
library("annotate")
library("hgu95av2.db")
library("GO.db")
library(ALL)
library(genefilter)
data(ALL)
data(geneList)
library(package = affyLib, character.only = TRUE)
# extract mappings
gofile <- hgu95av2GO
mapped_genes <- mappedkeys(gofile)
annots_all <- as.list(gofile[mapped_genes])
annots_all <- setNames(sapply(1:length(annots_all), function(x) names(annots_all[[x]])[which(sapply(annots_all[[x]], function(y) y[["Ontology"]] == "BP"))]),
                       names(annots_all))
# extract all genes from ALL
genes <- rownames(ALL@assayData[["exprs"]])
genes <- as.data.frame(cbind(id=genes, n=sapply(genes, function(x) length(annots_all[[x]]))))
genes$n <- as.numeric(genes$n)
# delete all genes with no GO terms
genes <- genes[!genes$n==0,]

genes$stat <- rnorm(length(genes$id), mean=predict(real_lm, data.frame(nmemberlog10=log10(genes$n))),
                             sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(genes$n))))
# extract genelist


# save for MWU
tosave <- genes[,c(1,3)]
write.table(tosave, file="./MWU/all_topgo.csv", row.names=F, quote=F, sep=",")
tosave <- data.frame("gene"=names(annots_all), "GO"=sapply(1:length(annots_all), 
                                                           function(x) paste0(annots_all[[x]],collapse = ";")))
write.table(tosave, file="./MWU/all_annots.csv", sep="\t", row.names = F, col.names = F, quote=F)




### LM ALL DATA

B <- 100
ks_job_lin_all <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  genes$stat <- rnorm(length(genes$id), mean=predict(real_lm, data.frame(nmemberlog10=log10(genes$n))),
                      sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(genes$n))))
  
  geneList <- setNames(genes$stat, genes$id)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots_all, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}

### LOESS ALL data

ks_job_loess_all <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  genes$stat <- rnorm(length(genes$id), mean=predict(lo, log10(genes$n)),
                    sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(genes$n))))
  
  geneList <- setNames(genes$stat, genes$id)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots_all, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}



################################################### MWU DATA
# get MWU datafiles
mwu_data <- read.table("./MWU/heats.csv", header=T, sep=",")
mwu_annots <- read.table("./MWU/amil_defog_iso2go.tab", header=T)

annots_mwu <- setNames(sapply(mwu_annots$V2, function(x) strsplit(x, ";",fixed=T)), mwu_annots$V1)

mwu_data$n <- sapply(mwu_data$gene, function(x) length(annots_mwu[[x]]))
# delte with no GO term
mwu_data <- mwu_data[mwu_data$n != 0,]


B <- 100
ks_job_lin_heat <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  mwu_data$stat <- rnorm(length(mwu_data$gene), mean=predict(real_lm, data.frame(nmemberlog10=log10(mwu_data$n))),
                      sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(mwu_data$n))))
  
  geneList <- setNames(mwu_data$stat, mwu_data$gene)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots_mwu, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}

### LOESS ALL data

ks_job_loess_heat <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  mwu_data$stat <- rnorm(length(mwu_data$gene), mean=predict(lo, log10(mwu_data$n)),
                    sd=rnorm(1, mean=.9, sd=.01*exp(1.5*log10(mwu_data$n))))
  
  geneList <- setNames(mwu_data$stat, mwu_data$gene)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots_mwu, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}


additional_algo <- list(ks_all = list(loess = ks_job_loess_all, linear = ks_job_lin_all),
                        ks_heat = list(loess = ks_job_loess_heat, linear = ks_job_lin_heat),
                        mwu_all = list(loess = mwu_job_loess_all, linear = mwu_job_lin_all),
                        mwu_heat = list(loess = mwu_job_loess_heat, linear = mwu_job_lin_heat))

save(additional_algo, file="additional_algo_data.RData")
