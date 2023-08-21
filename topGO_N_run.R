library(foreach)
library(dplyr)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-4)
registerDoParallel(cl)


library(topGO)

resultandseed <- function(result=NULL,seed=NULL)
{
  me <- list(
    result = result,
    seed = seed
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"resultandseed")
  return(me)
}

B <- 1000

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
  return(results)
}

string <- paste(sapply(1:1000, function(x) paste("ks_job[[",x,"]]$result[[1]],", sep="")), collapse=" ")
# string <- paste("GenTable(GOdat,", string, "topNodes=4081")

ks_results_classic <- GenTable(GOdat, ks_job[[1]][["result"]][[1]], topNodes=4081) %>% arrange(GO.ID)

ks_results_elim <- GenTable(GOdat, ks_job[[1]][["result"]][[3]], topNodes=4081) %>% arrange(GO.ID)

ks_results_weight01 <- GenTable(GOdat, ks_job[[1]][["result"]][[5]], topNodes=4081) %>% arrange(GO.ID)

ks_results_lea <- GenTable(GOdat, ks_job[[1]][["result"]][[7]], topNodes=4081) %>% arrange(GO.ID)

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

ks_results_classic <- get_all_its(ks_results_classic, ks_job, 1)
ks_results_elim <- get_all_its(ks_results_elim, ks_job, 3)
ks_results_weight01 <- get_all_its(ks_results_weight01, ks_job, 5)
ks_results_lea <- get_all_its(ks_results_lea, ks_job, 7)


# eval(sapply(1:10, function(x) paste("ks_job[[",x,"]][['result']][[1]]", sep="")))
# 



# what do we want
# median, mean number of other membership in the go term
# median, mean fc, pval

# per transcript
# number of memmbership in enriched go terms, number of membership

for(i in ks_results$GO.ID){
  transcripts <- 
}


