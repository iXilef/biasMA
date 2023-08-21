# make additional GO with more datasets
## preparation from topGO_runs.R

setwd("C:/Users/iXilef/OneDrive/Master/initial")
load("subset.RData")
annotations <- read.table("./MWU_clone/annotations.BP.go")
colnames(annotations) <- c("id", "go")

library(gageData)
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

B <- 100

runall <- function(godata, algor, test){
  alltests <- list()
  for(algo in algor){
    for(tests in test){
      alltests <- append(alltests, list(runTest(godata, algorithm = algo, statistic = tests), 
                                        name=paste(algo, "-", tests)))
    }
  }
  return(alltests)
}

ks_runs_orig_data <- c()

for(i in seq(0, 0.35, 0.05)){
  print(i)
  B <- 100
  intercept <- 1.24221*i
  ks_job_linfun <- foreach(1:B, .packages = c("topGO")) %dopar% {
    seed <- runif(1, 1, 10000000)
    set.seed(seed)
    x$log2FC <- rnorm(length(x$log2FC), mean=i*log10(x$n)-intercept, sd=0.9192379)
    
    geneList <- setNames(x$log2FC, x$transcript)
    GOdat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
                 gene2GO = annots, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
    ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
    results <- resultandseed()
    results$result <- ksrun
    results$seed <- c(seed, i)
    results$seq <- geneList
    return(results)
  }
  ks_runs_orig_data <- c(ks_runs_orig_data, ks_job_linfun)
}

## prepare other datasets
library(gageData)
# human data 
data(go.sets.hs)
data(go.subs.hs)

human <- go.sets.hs[go.subs.hs$BP]
names(human) <- sapply(names(human), function(x) unlist(strsplit(x, " ", fixed=T))[1])
human <- inverseList(human)

# yeast 
data(go.sets.sc)
data(go.subs.sc)

yeast <- go.sets.sc[go.subs.sc$BP]
names(yeast) <- sapply(names(yeast), function(x) unlist(strsplit(x, " ", fixed=T))[1])
yeast <- inverseList(yeast)


# rat 
data(go.sets.rn)
data(go.subs.rn)

rat <- go.sets.rn[go.subs.rn$BP]
names(rat) <- sapply(names(rat), function(x) unlist(strsplit(x, " ", fixed=T))[1])
rat <- inverseList(rat)

# mouse
data(go.sets.mm)
data(go.subs.mm)

mouse <- go.sets.mm[go.subs.mm$BP]
names(mouse) <- sapply(names(mouse), function(x) unlist(strsplit(x, " ", fixed=T))[1])
mouse <- inverseList(mouse)

par(mfrow=c(2,2))
hist(x$n, nclass=50)
hist(sapply(human, length), nclass=50)
hist(sapply(yeast, length), nclass=50)
hist(sapply(rat, length), nclass=50)
hist(sapply(mouse, length), nclass=50)

# choosing human and rat

human_data <- data.frame(id=names(human), n=sapply(human, length), stat=rnorm(length(human)))

rat_data <- data.frame(id=names(rat), n=sapply(rat, length), stat=rnorm(length(rat)))

# save data for MWU

# tosave <- rat_data[,c(1,3)]
# write.table(tosave, file="./MWU/rat.csv", row.names=F, quote=F, sep=",")
# tosave <- data.frame("gene"=names(rat), "GO"=sapply(1:length(rat), 
#                                                            function(x) paste0(rat[[x]],collapse = ";")))
# write.table(tosave, file="./MWU/rat_annots.csv", sep="\t", row.names = F, col.names = F, quote=F)
# 
# tosave <- human_data[,c(1,3)]
# write.table(tosave, file="./MWU/human.csv", row.names=F, quote=F, sep=",")
# tosave <- data.frame("gene"=names(human), "GO"=sapply(1:length(human), 
#                                                     function(x) paste0(human[[x]],collapse = ";")))
# write.table(tosave, file="./MWU/human_annots.csv", sep="\t", row.names = F, col.names = F, quote=F)



# rat run

ks_runs_rat_data <- c()

for(i in seq(.0, 0.35, 0.05)){
  print(i)
  B <- 100
  ks_job_linfun <- foreach(1:B, .packages = c("topGO")) %dopar% {
    seed <- runif(1, 1, 10000000)
    intercept <- 1.24221*i
    set.seed(seed)
    n <- rat_data$n
    rat_data$stat <- rnorm(length(rat_data$stat), mean=i*log10(rat_data$n)-intercept, sd=0.9192379)
    
    geneList <- setNames(rat_data$stat, rat_data$id)
    GOdat_rat <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
                 gene2GO = rat, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
    ksrun <- runall(godata=GOdat_rat, algor=whichAlgorithms()[-c(3,6)], test="ks")
    results <- resultandseed()
    results$result <- ksrun
    results$seed <- c(seed, i)
    results$seq <- geneList
    return(results)
  }
  ks_runs_rat_data <- c(ks_runs_rat_data, ks_job_linfun)
}

# human run

ks_runs_human_data <- c()

for(i in seq(0, 0.35, 0.05)){
  print(i)
  B <- 100
  intercept <- 1.24221*i
  ks_job_linfun <- foreach(1:B, .packages = c("topGO")) %dopar% {
    seed <- runif(1, 1, 10000000)
    set.seed(seed)
    human_data$stat <- rnorm(length(human_data$stat), mean=i*log10(human_data$n)-intercept, sd=0.9192379)
    
    geneList <- setNames(human_data$stat, human_data$id)
    GOdat_human <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
                 gene2GO = human, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)
    ksrun <- runall(GOdat_human, algor=whichAlgorithms()[-c(3,6)], test="ks")
    results <- resultandseed()
    results$result <- ksrun
    results$seed <- c(seed, i)
    results$seq <- geneList
    return(results)
  }
  ks_runs_human_data <- c(ks_runs_human_data, ks_job_linfun)
}

auxdata <- list(original = list(annotation=annots, data=x, GODat=GOdat),
                human = list(annotation = human, data=human_data, GOdat=GOdat_human),
                rat= list(annotation=rat, data=rat_data, GODat=GOdat_rat))
