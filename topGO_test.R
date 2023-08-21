setwd("C:/Users/Flix/OneDrive/Master/initial")
load("subset.RData")
annotations <- read.table("./MWU_clone/annotations.BP.go")
colnames(annotations) <- c("id", "go")

library(topGO)
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

fishrun <- runall(GOdat, algor=whichAlgorithms(), test="fisher") 
ksrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="ks")
trun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="t")
# gtestrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="globaltest")
# sumrun <- runall(GOdat, algor=whichAlgorithms()[-c(3,6)], test="sum")

fishscores <- as.data.frame(cbind(score(fishrun[[1]]), score(fishrun[[3]]),
                                  score(fishrun[[5]]), score(fishrun[[7]]),
                                  score(fishrun[[9]]), score(fishrun[[11]])))
colnames(fishscores) <- c("classic", "elim", "weight", "weight01", "lea", "parentchild")



ks_results <- GenTable(GOdat, classic=ksrun[[1]], elim = ksrun[[3]], weight01 = ksrun[[5]], lea = ksrun[[7]], topNodes = 3091)

ks_results <- as.data.frame(type.convert(ks_results, as.is=T)) %>% arrange(GO.ID)


par(mfrow=c(2,3))
for(i in c(1,3,5,7,9,11)){hist(fishrun[[i]]@score, nclass=100, main=fishrun[[i+1]])}
par(mfrow=c(2,2))
for(i in c(1,3,5,7)){hist(ksrun[[i]]@score, nclass=100, main=ksrun[[i+1]])}

for(i in c(1,3,5,7)){hist(trun[[i]]@score, nclass=100, main=trun[[i+1]])}

for(i in c(1,3,7,9)){hist(fishrun[[i]]@score, nclass=100, main=fishrun[[i+1]])}





pvalc <- score(fishrun[[1]])
pvalel <- score(fishrun[[3]])[names(pvalc)]

gstat <- termStat(GOdat, names(pvalc))

gSize <- gstat$Annotated/max(gstat$Annotated) * 4
# gCol <- colMap(gstat$Significant)

colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))




# plot color palette
plot(rep(1,378),col=(colfunc(378)), pch=19,cex=2)
plot(pvalc, pvalel, xlab="pval classic", ylab="pval elim", pch=19, cex=gSize, col=(colfunc(length(unique(gstat$Significant)))))
####


## get plot from above with sum of all divided by
annotated <- genesInTerm(GOdat, usedGO(GOdat))

# get number of transcripts as members
storage <- NULL
for(i in fishRes$GO.ID){
  transcripts <- annotated[[i]]
  counts <- unlist(lapply(transcripts, function(x) length(annots[[x]])))
  # pvalcount <- unlist()
  storage <- rbind(storage, cbind(transcripts, counts))
}
storage <- as.data.frame(storage)

storage <- storage[!duplicated(storage$transcripts),]

counts <- c()

for(i in row.names(gstat)){
  transcripts <- annotated[[i]]
  counts <- c(counts, sum(as.numeric(storage$counts[which(storage$transcripts %in% transcripts)]))/length(transcripts))
}


gstat$cummemoftranscripts <- counts

length(unique(gstat$cummemoftranscripts))

colfunc<-colorRampPalette(c("blue","green"))

heatcols <- heat.colors(3866)

# plot color palette
plot(rep(1,3866),col=(gstat$cummemoftranscripts), pch=19,cex=2)
plot(pvalc, pvalel, xlab="pval classic", ylab="pval elim", pch=19, cex=gstat$cummemoftranscripts/30,
     col=unlist(sapply(gstat$cummemoftranscripts, function(x) heatcols[x])), xlim=c(0,0.1), ylim=c(0,0.1))

library(ggplot2)
gstat$pvalc <- pvalc
gstat$pvalel <- pvalel
ggplot(data=gstat, aes(x=pvalc, pvalel, col=log(cummemoftranscripts), size=exp(cummemoftranscripts)))+
  geom_point()+
  scale_colour_gradient(low="red", high="green")

 ####


showGroupDensity(GOdat, row.names(gstat)[1], ranks = TRUE)
# todo: 


# get transcript membership in significant groups

for(i in )



countthemall <- function(x){
  count <- 0
  for(i in 1:length(annots)){
    if(x %in% annots[[i]]){count <- count+1}
  return(count)}}

counties <- c()
for(go in row.names(fishscores)){
  count <- 0
  for(i in 1:length(annotations$id)){
    if(go %in% strsplit(annotations$go[i], ";", fixed=T)){
      count <- count+1
    }
  }
  counties <- c(counties, count)
}

lapply(row.names(fishscores), function(x) return(countthemall(x)))

fishscores$gene_num <- 