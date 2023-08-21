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

# create named vector of pvalues
geneList <- setNames(10^(-x$minuslog10pval), x$transcript)


GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
             gene2GO = annots, geneSel = function(x) x < .1, nodeSize = 10)

# runTest(GOdat, algorithm = "elim", statistic = "sum")

# do all the tests
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



# run all algorithms per stat test

# fishrun <- runall(GOdat, algor=whichAlgorithms(), test="fisher") 
ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
# trun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="t")
# gtestrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="globaltest")
# sumrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="sum")


ks_results <- GenTable(GOdat, classic=ksrun[[1]], elim = ksrun[[3]], weight01 = ksrun[[5]], lea = ksrun[[7]], topNodes = 3091)

ks_results <- as.data.frame(type.convert(ks_results, as.is=T)) %>% arrange(GO.ID)


######################--------------Shuffle Mode ------------#################

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

ks_job <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  geneList <- setNames(sample(x$log2FC), x$transcript)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots, geneSel = function(x) x < 0.5 | x > 1, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}

# extract first result
N <- length(ks_job[[1]]$result[[1]]@score)

ks_results_classic <- GenTable(GOdat, ks_job[[1]][["result"]][[1]], topNodes=N) %>% arrange(GO.ID)

ks_results_elim <- GenTable(GOdat, ks_job[[1]][["result"]][[3]], topNodes=N) %>% arrange(GO.ID)

ks_results_weight01 <- GenTable(GOdat, ks_job[[1]][["result"]][[5]], topNodes=N) %>% arrange(GO.ID)

ks_results_lea <- GenTable(GOdat, ks_job[[1]][["result"]][[7]], topNodes=N) %>% arrange(GO.ID)

# function to extract the rest

get_all_its <- function(results, iterations, algo){
  for(i in 1:length(iterations)){
    print(i)
    results <- cbind(results, ks_job[[i]][["result"]][[algo]]@score)
  }
  results <- as.data.frame(type.convert(results, as.is=T))
  results$mean <- apply(results[,7:106], 1, mean)
  results$median <- apply(results[,7:106], 1, median)
  return(results)
}

# make random resampling run
ks_job_norm <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  geneList <- setNames(rnorm(length(x$log2FC), 0, var(x$log2FC)), x$transcript)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots, geneSel = function(x) x < 0.5 | x > 1, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- seed
  results$seq <- geneList
  return(results)
}


# extract first result
N <- length(ks_job[[1]]$result[[1]]@score)

ks_results_classic <- GenTable(GOdat, ks_job[[1]][["result"]][[1]], topNodes=N) %>% arrange(GO.ID)

ks_results_elim <- GenTable(GOdat, ks_job[[1]][["result"]][[3]], topNodes=N) %>% arrange(GO.ID)

ks_results_weight01 <- GenTable(GOdat, ks_job[[1]][["result"]][[5]], topNodes=N) %>% arrange(GO.ID)

ks_results_lea <- GenTable(GOdat, ks_job[[1]][["result"]][[7]], topNodes=N) %>% arrange(GO.ID)

# extract first result
N <- length(ks_job_norm[[1]]$result[[1]]@score)

ks_results_classic_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[1]], topNodes=N) %>% arrange(GO.ID)

ks_results_elim_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[3]], topNodes=N) %>% arrange(GO.ID)

ks_results_weight01_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[5]], topNodes=N) %>% arrange(GO.ID)

ks_results_lea_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[7]], topNodes=N) %>% arrange(GO.ID)


# function to extract the rest

get_all_its <- function(results, iterations, algo){
  for(i in 1:length(iterations)){
    print(i)
    results <- cbind(results, iterations[[i]][["result"]][[algo]]@score)
  }
  results <- as.data.frame(type.convert(results, as.is=T))
  results$mean <- apply(results[,7:106], 1, mean)
  results$median <- apply(results[,7:106], 1, median)
  return(results)
}

# get all pvals into one dataframe

ks_results_classic <- get_all_its(ks_results_classic, ks_job, 1)
ks_results_elim <- get_all_its(ks_results_elim, ks_job, 3)
ks_results_weight01 <- get_all_its(ks_results_weight01, ks_job, 5)
ks_results_lea <- get_all_its(ks_results_lea, ks_job, 7)


ks_results_classic_norm <- get_all_its(ks_results_classic_norm, ks_job_norm, 1)
ks_results_elim_norm <- get_all_its(ks_results_elim_norm, ks_job_norm, 3)
ks_results_weight01_norm <- get_all_its(ks_results_weight01_norm, ks_job_norm, 5)
ks_results_lea_norm <- get_all_its(ks_results_lea_norm, ks_job_norm, 7)


B <- 1000
# make resampling from Normal
ks_job_norm_2 <- foreach(1:B, .packages = c("topGO")) %dopar% {
  seed <- runif(1, 1, 10000000)
  set.seed(seed)
  # get side of special treatment
  side <- sample(1:2, 1)
  # get percentage for base distribution
  percentage <- as.integer(runif(1, 5, 20))
  # get mean
  if(side==1){
    mean.special <- runif(1, .5, 1)
  } else {mean.special <- runif(1, -1, -.5)}
  # split data randomly
  x_indices <- sample(c(1:length(x$log2FC)), as.integer(length(x$log2FC)*percentage/100))
  x_big <- x[-x_indices,]
  x_big$special <- 0
  x_small <- x[x_indices,]
  x_small$special <- 1
  
  # assign norm values
  x_big$log2FC <- rnorm(length(x_big$log2FC), 0, var(x$log2FC))
  x_small$log2FC <- rnorm(length(x_small$log2FC), mean.special, var(x$log2FC))
  
  x <- rbind(x_big, x_small)
  
  geneList <- setNames(x$log2FC, x$transcript)
  GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
               gene2GO = annots, geneSel = function(x) x < 0.5 | x > 1, nodeSize = 10)
  ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
  results <- resultandseed()
  results$result <- ksrun
  results$seed <- list("seed"=seed, "mean" = mean.special, "percentage" = percentage, "values" = x)
  results$seq <- geneList
  
  return(results)
}
save(ks_job_norm_2, file = "./ks_job_norm_2.RData")
