library(topGO)
library(dplyr)

setwd("C:/Users/iXilef/OneDrive/Master/initial")

load("additional_algo_data.RData")
load("modelruns_1.RData")
load("./MWU/add_data_algo4.RData")

# annotation original files

annotations <- read.table("./MWU_clone/annotations.BP.go")
colnames(annotations) <- c("id", "go")
annots <- setNames(annotations$go, annotations$id)
for (i in 1:length(annots)){
  annots[i] <- strsplit(annots[[i]], ";", fixed=T)
  # print(i)
}

# annotation ALL data

library(ALL)
library(genefilter)
data(ALL)
data(geneList)
library(package = affyLib, character.only = TRUE)
# extract mappings
gofile <- hgu95av2GO
mapped_genes <- mappedkeys(gofile)
annots_all <- as.list(gofile[mapped_genes])
annots_all <- setNames(sapply(1:length(annots_all), function(x) names(annots_all[[x]])),
                       names(annots_all))
# extract all genes from ALL
genes <- rownames(ALL@assayData[["exprs"]])
genes <- as.data.frame(cbind(id=genes, n=sapply(genes, function(x) length(annots_all[[x]]))))
genes$n <- as.numeric(genes$n)
# delete all genes with no GO terms
genes <- genes[!genes$n==0,]

genes$stat <- 1

# annotation MWU data

mwu_data <- read.table("./MWU/heats.csv", header=T, sep=",")
mwu_annots <- read.table("./MWU/amil_defog_iso2go.tab", header=T)

annots_mwu <- setNames(sapply(mwu_annots$V2, function(x) strsplit(x, ";",fixed=T)), mwu_annots$V1)

mwu_data$n <- sapply(mwu_data$gene, function(x) length(annots_mwu[[x]]))
# delte with no GO term
mwu_data <- mwu_data[mwu_data$n != 0,]


# GO dat
geneList <- setNames(genes$stat, genes$id)
GOdat_all <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
             gene2GO = annots_all, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)



geneList <- setNames(mwu_data$logP, mwu_data$gene)
GOdat_mwu <- new("topGOdata", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO,
             gene2GO = annots_mwu, geneSel = function(x) x < 0.5 | x > .5, nodeSize = 10)


#extract MWU infos 
B <- 100

mwu_lin <- as.data.frame(cbind(model_runs$mwu_lin[[1]]$result$term, as.numeric(model_runs$mwu_lin[[1]]$result$nseqs), model_runs$mwu_lin[[1]]$result$pval))
for(i in 2:length(model_runs$mwu_lin)){
  mwu_lin <- cbind(mwu_lin, model_runs$mwu_lin[[i]]$result$pval)
}

colnames(mwu_lin) <- c("term", "nseqs", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))
mwu_lin$term_single <- sapply(mwu_lin$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])


mwu_lin_trans <- as.data.frame(cbind(unique(model_runs$mwu_lin[[1]]$seq$seq), model_runs$mwu_lin[[1]]$seq$value[!duplicated(model_runs$mwu_lin[[1]]$seq$seq)]))
for(i in 2:length(model_runs$mwu_lin)){
  mwu_lin_trans <- cbind(mwu_lin_trans, model_runs$mwu_lin[[i]]$seq$value[!duplicated(model_runs$mwu_lin[[i]]$seq$seq)])
}

colnames(mwu_lin_trans) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))

mwu_lin_trans_n <- cbind(mwu_lin_trans$transcript, sapply(mwu_lin_trans$transcript, function(x) length(which(mwu_lin$result1[mwu_lin$term_single %in% annots[[x]]] < 0.05))))
for(i in 3:dim(mwu_lin_trans)[2]+1){
  mwu_lin_trans_n <- cbind(mwu_lin_trans_n, sapply(mwu_lin_trans$transcript, function(x) length(which(mwu_lin[,i][mwu_lin$term_single %in% annots[[x]]] < 0.05)))) 
}

mwu_lin_trans_n <- as.data.frame(type.convert(mwu_lin_trans_n, as.is=T))
colnames(mwu_lin_trans_n) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))

  
#


mwu_loess <- as.data.frame(cbind(model_runs$mwu_loess[[1]]$result$term, as.numeric(model_runs$mwu_loess[[1]]$result$nseqs), model_runs$mwu_loess[[1]]$result$pval))
for(i in 2:length(model_runs$mwu_loess)){
  mwu_loess <- cbind(mwu_loess, model_runs$mwu_loess[[i]]$result$pval)
}

colnames(mwu_loess) <- c("term", "nseqs", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))
mwu_loess$term_single <- sapply(mwu_loess$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])


mwu_loess_trans <- as.data.frame(cbind(unique(model_runs$mwu_loess[[1]]$seq$seq), model_runs$mwu_loess[[1]]$seq$value[!duplicated(model_runs$mwu_loess[[1]]$seq$seq)]))
for(i in 2:length(model_runs$mwu_loess)){
  mwu_loess_trans <- cbind(mwu_loess_trans, model_runs$mwu_loess[[i]]$seq$value[!duplicated(model_runs$mwu_loess[[i]]$seq$seq)])
}

colnames(mwu_loess_trans) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))

mwu_loess_trans_n <- cbind(mwu_loess_trans$transcript, sapply(mwu_loess_trans$transcript, function(x) length(which(mwu_loess$result1[mwu_loess$term_single %in% annots[[x]]] < 0.05))))
for(i in 3:dim(mwu_lin_trans)[2]+1){
  mwu_loess_trans_n <- cbind(mwu_loess_trans_n, sapply(mwu_loess_trans$transcript, function(x) length(which(mwu_loess[,i][mwu_loess$term_single %in% annots[[x]]] < 0.05)))) 
}

colnames(mwu_loess_trans_n) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))
mwu_loess_trans_n <- as.data.frame(type.convert(mwu_loess_trans_n, as.is=T))

  
## MWU heat

mwu_lin_heat <- as.data.frame(cbind(mwu_add$heat$linear[[1]]$result$term, as.numeric(mwu_add$heat$linear[[1]]$result$nseqs), mwu_add$heat$linear[[1]]$result$pval))
for(i in 2:length(mwu_add$all$linear)){
  print(i)
  mwu_lin_heat <- cbind(mwu_lin_heat, mwu_add$heat$linear[[i]]$result$pval)
}

colnames(mwu_lin_heat) <- c("term", "nseqs", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))
mwu_lin_heat$term_single <- sapply(mwu_lin_heat$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])


mwu_lin_heat_trans <- as.data.frame(cbind(unique(mwu_add$heat$linear[[1]]$seq$seq), mwu_add$heat$linear[[1]]$seq$value[!duplicated(mwu_add$heat$linear[[1]]$seq$seq)]))
for(i in 2:length(mwu_add$heat$linear)){
  mwu_lin_heat_trans <- cbind(mwu_lin_heat_trans, mwu_add$heat$linear[[i]]$seq$value[!duplicated(mwu_add$heat$linear[[i]]$seq$seq)])
}

colnames(mwu_lin_heat_trans) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))

mwu_lin_heat_trans_n <- cbind(mwu_lin_heat_trans$transcript, sapply(mwu_lin_heat_trans$transcript, function(x) length(which(mwu_lin_heat$result1[mwu_lin_heat$term_single %in% annots_mwu[[x]]] < 0.05))))
for(i in 3:dim(mwu_lin_heat_trans)[2]+1){
  mwu_lin_heat_trans_n <- cbind(mwu_lin_heat_trans_n, sapply(mwu_lin_heat_trans$transcript, function(x) length(which(mwu_lin_heat[,i][mwu_lin_heat$term_single %in% annots_mwu[[x]]] < 0.05))))
}

colnames(mwu_lin_heat_trans_n) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))
mwu_lin_heat_trans_n <- as.data.frame(type.convert(mwu_lin_heat_trans_n, as.is=T))
## MWU ALL

mwu_lin_all <- as.data.frame(cbind(mwu_add$all$linear[[1]]$result$term, as.numeric(mwu_add$all$linear[[1]]$result$nseqs), mwu_add$all$linear[[1]]$result$pval))
for(i in 2:length(mwu_add$all$linear)){
  print(i)
  mwu_lin_all <- cbind(mwu_lin_all, mwu_add$all$linear[[i]]$result$pval)
}

colnames(mwu_lin_all) <- c("term", "nseqs", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))
mwu_lin_all$term_single <- sapply(mwu_lin_all$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])




mwu_lin_all_trans <- as.data.frame(cbind(unique(mwu_add$all$linear[[1]]$seq$seq), mwu_add$all$linear[[1]]$seq$value[!duplicated(mwu_add$all$linear[[1]]$seq$seq)]))
for(i in 2:length(mwu_add$all$linear)){
  mwu_lin_all_trans <- cbind(mwu_lin_all_trans, mwu_add$all$linear[[i]]$seq$value[!duplicated(mwu_add$all$linear[[i]]$seq$seq)])
}

colnames(mwu_lin_all_trans) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))

mwu_lin_all_trans_n <- cbind(mwu_lin_all_trans$transcript, sapply(mwu_lin_all_trans$transcript, function(x) length(which(mwu_lin_all$result1[mwu_lin_all$term_single %in% annots_all[[x]]] < 0.05))))
for(i in 3:dim(mwu_lin_all_trans)[2]+1){
  mwu_lin_all_trans_n <- cbind(mwu_lin_all_trans_n, sapply(mwu_lin_all_trans$transcript, function(x) length(which(mwu_lin_all[,i+1][mwu_lin_all$term_single %in% annots_all[[x]]] < 0.05)))) 
}

colnames(mwu_lin_all_trans_n) <- c("transcript", unlist(lapply(1:B, function(x) paste("result", x, sep=""))))
mwu_lin_all_trans_n <- as.data.frame(type.convert(mwu_lin_all_trans_n, as.is=T))

######################### GO runs
get_all_its <- function(results, iterations, algo){
  for(i in 1:length(iterations)){
    print(i)
    results <- cbind(results, iterations[[i]][["result"]][[algo]]@score)
  }
  results <- as.data.frame(type.convert(results, as.is=T))
  results$mean <- apply(results[,7:(length(iterations)+6)], 1, mean)
  results$median <- apply(results[,7:(length(iterations)+6)], 1, median)
  return(results)
}
## extract first from ks


N <- length(model_runs$ks_lin[[1]]$result[[1]]@score)

ks_results_classic_lin <- GenTable(GOdat, model_runs$ks_lin[[1]]$result[[1]], topNodes=N) %>% arrange(GO.ID)
ks_results_classic_lin <- get_all_its(ks_results_classic_lin, model_runs$ks_lin, 1) 

ks_results_elim_lin <- GenTable(GOdat, model_runs$ks_lin[[1]]$result[[3]], topNodes=N) %>% arrange(GO.ID)
ks_results_elim_lin <- get_all_its(ks_results_elim_lin, model_runs$ks_lin, 3) 


ks_results_weight01_lin <- GenTable(GOdat, model_runs$ks_lin[[1]]$result[[5]], topNodes=N) %>% arrange(GO.ID)
ks_results_weight01_lin <- get_all_its(ks_results_weight01_lin, model_runs$ks_lin, 5) 


ks_results_lea_lin <- GenTable(GOdat, model_runs$ks_lin[[1]]$result[[7]], topNodes=N) %>% arrange(GO.ID)
ks_results_lea_lin <- get_all_its(ks_results_lea_lin, model_runs$ks_lin, 7) 

ks_trans <- as.data.frame(cbind(names(model_runs$ks_lin[[1]]$seq),
                                model_runs$ks_lin[[1]]$seq))
for(i in 2:length(model_runs$ks_lin)){
  ks_trans <- cbind(ks_trans, model_runs$ks_lin[[i]]$seq)
}
colnames(ks_trans) <- c("transcript", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

### extract first results from all

N <- length(additional_algo$ks_all$linear[[1]]$result[[1]]@score)

ks_results_classic_lin_all <- GenTable(GOdat_all, additional_algo$ks_all$linear[[1]]$result[[1]], topNodes=N) %>% arrange(GO.ID)
ks_results_classic_lin_all <- get_all_its(ks_results_classic_lin_all, additional_algo$ks_all$linear, 1) 


ks_results_elim_lin_all <- GenTable(GOdat_all, additional_algo$ks_all$linear[[1]]$result[[3]], topNodes=N) %>% arrange(GO.ID)
ks_results_elim_lin_all <- get_all_its(ks_results_elim_lin_all, additional_algo$ks_all$linear, 3) 


ks_results_weight01_lin_all <- GenTable(GOdat_all, additional_algo$ks_all$linear[[1]]$result[[5]], topNodes=N) %>% arrange(GO.ID)
ks_results_weight01_lin_all <- get_all_its(ks_results_weight01_lin_all, additional_algo$ks_all$linear, 5) 


ks_results_lea_lin_all <- GenTable(GOdat_all, additional_algo$ks_all$linear[[1]]$result[[7]], topNodes=N) %>% arrange(GO.ID)
ks_results_lea_lin_all <- get_all_its(ks_results_lea_lin_all, additional_algo$ks_all$linear, 7) 

all_trans <- as.data.frame(cbind(names(additional_algo$ks_all$linear[[1]]$seq),
                                 additional_algo$ks_all$linear[[1]]$seq))
for(i in 2:length(additional_algo$ks_all$linear)){
  all_trans <- cbind(all_trans, additional_algo$ks_all$linear[[i]]$seq)
}
colnames(all_trans) <- c("transcript", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

### extract first results from heat

N <- length(additional_algo$ks_heat$linear[[1]]$result[[1]]@score)

ks_results_classic_lin_heat<- GenTable(GOdat_mwu, additional_algo$ks_heat$linear[[1]]$result[[1]], topNodes=N) %>% arrange(GO.ID)
ks_results_classic_lin_heat <- get_all_its(ks_results_classic_lin_heat,additional_algo$ks_heat$linear, 1) 

ks_results_elim_lin_heat <- GenTable(GOdat_mwu, additional_algo$ks_heat$linear[[1]]$result[[3]], topNodes=N) %>% arrange(GO.ID)
ks_results_elim_lin_heat <- get_all_its(ks_results_elim_lin_heat, additional_algo$ks_heat$linear, 1) 

ks_results_weight01_lin_heat <- GenTable(GOdat_mwu, additional_algo$ks_heat$linear[[1]]$result[[5]], topNodes=N) %>% arrange(GO.ID)
ks_results_weight01_lin_heat <- get_all_its(ks_results_weight01_lin_heat, additional_algo$ks_heat$linear, 1) 

ks_results_lea_lin_heat <- GenTable(GOdat_mwu, additional_algo$ks_heat$linear[[1]]$result[[7]], topNodes=N) %>% arrange(GO.ID)
ks_results_lea_lin_heat <- get_all_its(ks_results_lea_lin_heat, additional_algo$ks_heat$linear, 1) 

heat_trans <- as.data.frame(cbind(names(additional_algo$ks_heat$linear[[1]]$seq),
                                  additional_algo$ks_heat$linear[[1]]$seq))
for(i in 2:length(additional_algo$ks_heat$linear)){
  heat_trans <- cbind(heat_trans, additional_algo$ks_heat$linear[[i]]$seq)
}
colnames(heat_trans) <- c("transcript", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

### get N from lin runs

get_trans_go <- function(results, iterations, algo, annotation){
  res <- data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq)
  for(i in 2:iterations){
    res <- cbind(res, results[[i]]$seq)
  }
  res_n <- data.frame(transcript=names(annotation))
  res_n <- left_join(res_n, data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq), by="transcript")
  # classic
  res_classic_n <- res_n
  res_classic_n$first <- sapply(annotation, function(x) length(which(ks_results_classic_lin$result1[which(ks_results_classic_lin$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_classic_n <- cbind(res_classic_n, sapply(annotation, function(x) length(which(ks_results_classic_lin[,i][which(ks_results_classic_lin$GO.ID %in% x)] < 0.05))))
  }
  res_classic_n <- na.omit(res_classic_n)
  
  # elim 
  res_elim_n <- res_n
  res_elim_n$first <- sapply(annotation, function(x) length(which(ks_results_elim_lin$result1[which(ks_results_elim_lin$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_elim_n <- cbind(res_elim_n, sapply(annotation, function(x) length(which(ks_results_elim_lin[,i][which(ks_results_elim_lin$GO.ID %in% x)] < 0.05))))
  }
  res_elim_n <- na.omit(res_elim_n)
  
  # weight
  res_weight_n <- res_n
  res_weight_n$first <- sapply(annotation, function(x) length(which(ks_results_weight01_lin$result1[which(ks_results_weight01_lin$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_weight_n <- cbind(res_weight_n, sapply(annotation, function(x) length(which(ks_results_weight01_lin[,i][which(ks_results_weight01_lin$GO.ID %in% x)] < 0.05))))
  }
  res_weight_n <- na.omit(res_weight_n)
  
  # lea
  res_lea_n <- res_n
  res_lea_n$first <- sapply(annotation, function(x) length(which(ks_results_lea_lin$result1[which(ks_results_lea_lin$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_lea_n <- cbind(res_lea_n, sapply(annotation, function(x) length(which(ks_results_lea_lin[,i][which(ks_results_lea_lin$GO.ID %in% x)] < 0.05))))
  }
  res_lea_n <- na.omit(res_lea_n)
  
  res <- res[which(res$transcript %in% res_classic_n$transcript),]
  
  return(list(res, res_classic_n, res_elim_n, res_weight_n, res_lea_n))
}

res_lin_n <- get_trans_go(model_runs$ks_lin, 100, 0, annots)

#### get N from ALL runs

get_trans_go <- function(results, iterations, algo, annotation){
  res <- data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq)
  for(i in 2:iterations){
    res <- cbind(res, results[[i]]$seq)
  }
  res_n <- data.frame(transcript=names(annotation))
  res_n <- left_join(res_n, data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq), by="transcript")
  # classic
  res_classic_n <- res_n
  res_classic_n$first <- sapply(annotation, function(x) length(which(ks_results_classic_lin_all$result1[which(ks_results_classic_lin_all$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_classic_n <- cbind(res_classic_n, sapply(annotation, function(x) length(which(ks_results_classic_lin_all[,i][which(ks_results_classic_lin_all$GO.ID %in% x)] < 0.05))))
  }
  res_classic_n <- na.omit(res_classic_n)
  
  # elim 
  res_elim_n <- res_n
  res_elim_n$first <- sapply(annotation, function(x) length(which(ks_results_elim_lin_all$result1[which(ks_results_elim_lin_all$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_elim_n <- cbind(res_elim_n, sapply(annotation, function(x) length(which(ks_results_elim_lin_all[,i][which(ks_results_elim_lin_all$GO.ID %in% x)] < 0.05))))
  }
  res_elim_n <- na.omit(res_elim_n)
  
  # weight
  res_weight_n <- res_n
  res_weight_n$first <- sapply(annotation, function(x) length(which(ks_results_weight01_lin_all$result1[which(ks_results_weight01_lin_all$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_weight_n <- cbind(res_weight_n, sapply(annotation, function(x) length(which(ks_results_weight01_lin_all[,i][which(ks_results_weight01_lin_all$GO.ID %in% x)] < 0.05))))
  }
  res_weight_n <- na.omit(res_weight_n)
  
  # lea
  res_lea_n <- res_n
  res_lea_n$first <- sapply(annotation, function(x) length(which(ks_results_lea_lin_all$result1[which(ks_results_lea_lin_all$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_lea_n <- cbind(res_lea_n, sapply(annotation, function(x) length(which(ks_results_lea_lin_all[,i][which(ks_results_lea_lin_all$GO.ID %in% x)] < 0.05))))
  }
  res_lea_n <- na.omit(res_lea_n)
  
  res <- res[which(res$transcript %in% res_classic_n$transcript),]
  
  return(list(res, res_classic_n, res_elim_n, res_weight_n, res_lea_n))
}

res_lin_all_n <- get_trans_go(additional_algo$ks_all$linear, 100, 0, annots_all)

### get N from HEAT runs

get_trans_go <- function(results, iterations, algo, annotation){
  res <- data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq)
  for(i in 2:iterations){
    res <- cbind(res, results[[i]]$seq)
  }
  res_n <- data.frame(transcript=names(annotation))
  res_n <- left_join(res_n, data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq), by="transcript")
  # classic
  res_classic_n <- res_n
  res_classic_n$first <- sapply(annotation, function(x) length(which(ks_results_classic_lin_heat$result1[which(ks_results_classic_lin_heat$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_classic_n <- cbind(res_classic_n, sapply(annotation, function(x) length(which(ks_results_classic_lin_heat[,i][which(ks_results_classic_lin_heat$GO.ID %in% x)] < 0.05))))
  }
  res_classic_n <- na.omit(res_classic_n)
  
  # elim 
  res_elim_n <- res_n
  res_elim_n$first <- sapply(annotation, function(x) length(which(ks_results_elim_lin_heat$result1[which(ks_results_elim_lin_heat$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_elim_n <- cbind(res_elim_n, sapply(annotation, function(x) length(which(ks_results_elim_lin_heat[,i][which(ks_results_elim_lin_heat$GO.ID %in% x)] < 0.05))))
  }
  res_elim_n <- na.omit(res_elim_n)
  
  # weight
  res_weight_n <- res_n
  res_weight_n$first <- sapply(annotation, function(x) length(which(ks_results_weight01_lin_heat$result1[which(ks_results_weight01_lin_heat$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_weight_n <- cbind(res_weight_n, sapply(annotation, function(x) length(which(ks_results_weight01_lin_heat[,i][which(ks_results_weight01_lin_heat$GO.ID %in% x)] < 0.05))))
  }
  res_weight_n <- na.omit(res_weight_n)
  
  # lea
  res_lea_n <- res_n
  res_lea_n$first <- sapply(annotation, function(x) length(which(ks_results_lea_lin_heat$result1[which(ks_results_lea_lin_heat$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_lea_n <- cbind(res_lea_n, sapply(annotation, function(x) length(which(ks_results_lea_lin_heat[,i][which(ks_results_lea_lin_heat$GO.ID %in% x)] < 0.05))))
  }
  res_lea_n <- na.omit(res_lea_n)
  
  res <- res[which(res$transcript %in% res_classic_n$transcript),]
  
  return(list(res, res_classic_n, res_elim_n, res_weight_n, res_lea_n))
}

res_lin_heat_n <- get_trans_go(additional_algo$ks_heat$linear, 100, 0, annots_mwu)

## save into one object

results_lin <- list(classic=list(pval=ks_results_classic_lin,
                                       seq=ks_trans,
                                       seq_n=res_lin_n[[2]]),
                          
                          weight=list(pval=ks_results_weight01_lin,
                                      seq=ks_trans,
                                      seq_n=res_lin_n[[4]]),
                          
                          elim=list(pval=ks_results_elim_lin,
                                    seq=ks_trans,
                                    seq_n=res_lin_n[[3]]),
                          
                          lea=list(pval=ks_results_lea_lin,
                                   seq=ks_trans,
                                   seq_n=res_lin_n[[5]]),
                    
                          mwu=list(pval=as.data.frame(type.convert(mwu_lin, as.is=T)),
                                   seq=mwu_lin_trans,
                                   seq_n=mwu_lin_trans_n))

results_lin_all <- list(classic=list(pval=ks_results_classic_lin_all,
                                 seq=all_trans,
                                 seq_n=res_lin_all_n[[2]]),
                    
                    weight=list(pval=ks_results_weight01_lin_all,
                                seq=all_trans,
                                seq_n=res_lin_all_n[[4]]),
                    
                    elim=list(pval=ks_results_elim_lin_all,
                              seq=all_trans,
                              seq_n=res_lin_all_n[[3]]),
                    
                    lea=list(pval=ks_results_lea_lin_all,
                             seq=all_trans,
                             seq_n=res_lin_all_n[[5]]),
                    
                    mwu=list(pval=as.data.frame(type.convert(mwu_lin_all, as.is=T)),
                             seq=mwu_lin_all_trans,
                             seq_n=mwu_lin_all_trans_n))

results_lin_heat <- list(classic=list(pval=ks_results_classic_lin_heat,
                                 seq=heat_trans,
                                 seq_n=res_lin_heat_n[[2]]),
                    
                    weight=list(pval=ks_results_weight01_lin_heat,
                                seq=heat_trans,
                                seq_n=res_lin_heat_n[[4]]),
                    
                    elim=list(pval=ks_results_elim_lin_heat,
                              seq=heat_trans,
                              seq_n=res_lin_heat_n[[3]]),
                    
                    lea=list(pval=ks_results_lea_lin_heat,
                             seq=heat_trans,
                             seq_n=res_lin_heat_n[[5]]),
                    
                    mwu=list(pval=as.data.frame(type.convert(mwu_lin_heat, as.is=T)),
                             seq=mwu_lin_heat_trans,
                             seq_n=mwu_lin_heat_trans_n))

# get line objects 
# calculate moving median for all random runs
# mwu
get_lines_mwu <- function(fcvals, nvals, iterations){
  int.data <- data.frame(N = seq(0,300,1))
  for(i in 2:(iterations+1)){
    print(i)
    int.data <- cbind(int.data, sapply(0:300, function(x) median(as.numeric(fcvals[,i][which(nvals[,i] > x)]))))
  }
  int.data <- as.data.frame(int.data)
  return(int.data)                     
}

get_lines_ks <- function(fcvals, nvals, iterations){
  int.data <- data.frame(N = seq(0,300,1))
  for(i in 2:(iterations+1)){
    print(i)
    int.data <- cbind(int.data, sapply(0:300, function(x) median(as.numeric(fcvals[,i][which(nvals[,i+1] > x)]), na.rm=T)))
  }
  int.data <- as.data.frame(int.data)
  return(int.data)                     
}

mwu_lines <- get_lines_mwu(mwu_lin_trans, mwu_lin_trans_n, 100)
classic_lines <- get_lines_ks(results_lin$classic$seq, results_lin$classic$seq_n, 100)
elim_lines <- get_lines_ks(results_lin$elim$seq, results_lin$elim$seq_n, 100)
weight_lines <- get_lines_ks(results_lin$weight$seq, results_lin$weight$seq_n, 100)
lea_lines <- get_lines_ks(results_lin$lea$seq, results_lin$lea$seq_n, 100)

mwu_lines_all <- get_lines_mwu(mwu_lin_all_trans, mwu_lin_all_trans_n, 100)
classic_lines_all <- get_lines_ks(results_lin_all$classic$seq, results_lin_all$classic$seq_n, 100)
elim_lines_all <- get_lines_ks(results_lin_all$elim$seq, results_lin_all$elim$seq_n, 100)
weight_lines_all <- get_lines_ks(results_lin_all$weight$seq, results_lin_all$weight$seq_n, 100)
lea_lines_all <- get_lines_ks(results_lin_all$lea$seq, results_lin_all$lea$seq_n, 100)

mwu_lines_heat <- get_lines_mwu(mwu_lin_heat_trans, mwu_lin_heat_trans_n, 100)
classic_lines_heat <- get_lines_ks(results_lin_heat$classic$seq, results_lin_heat$classic$seq_n, 100)
elim_lines_heat <- get_lines_ks(results_lin_heat$elim$seq, results_lin_heat$elim$seq_n, 100)
weight_lines_heat <- get_lines_ks(results_lin_heat$weight$seq, results_lin_heat$weight$seq_n, 100)
lea_lines_heat <- get_lines_ks(results_lin_heat$lea$seq, results_lin_heat$lea$seq_n, 100)


