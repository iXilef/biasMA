library(dplyr)
library(ontologyIndex)
library(topGO)

setwd("C:/Users/iXilef/OneDrive/Master/initial")
load("./original_objects.RData")
load("./ks_runs.RData")
load("./N_KS_run.RData")
load("./ks_job_norm_2.RData")
load("./MWU/MWU_runs.RData")
load("./MWU/mwu_job_norm_2.RData")
annotations <- read.table("./MWU_clone/annotations.BP.go")
colnames(annotations) <- c("id", "go")

annots <- setNames(annotations$go, annotations$id)
for (i in 1:length(annots)){
  annots[i] <- strsplit(annots[[i]], ";", fixed=T)
  # print(i)
}

annotsrev <- inverseList(annots)


# prepare MWU objects
B <- 100

mwu_run$term_single <- sapply(mwu_run$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])

# MWU resample

mwu_results_rand <- as.data.frame(cbind(mwu_job[[1]]$result$term, as.numeric(mwu_job[[1]]$result$nseqs), mwu_job[[1]]$result$pval)) 
for(i in 2:length(mwu_job)){
  mwu_results_rand <- cbind(mwu_results_rand, mwu_job[[i]]$result$pval)
}

colnames(mwu_results_rand) <- c("term", "nseqs", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))
mwu_results_rand$term_single <- sapply(mwu_results_rand$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])


mwu_trans_rand <- as.data.frame(cbind(mwu_job[[1]]$seq$id, as.numeric(mwu_job[[1]]$seq$stat), mwu_job[[1]]$seq$random_fc)) 
for(i in 2:length(mwu_job)){
  mwu_trans_rand <- cbind(mwu_trans_rand, mwu_job[[i]]$seq$random_fc)
}

colnames(mwu_trans_rand) <- c("transcript", "original FC", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

# MWU norm 

mwu_results_norm <- as.data.frame(cbind(mwu_job_norm[[1]]$result$term, as.numeric(mwu_job_norm[[1]]$result$nseqs), mwu_job_norm[[1]]$result$pval)) 
for(i in 2:length(mwu_job_norm)){
  mwu_results_norm <- cbind(mwu_results_norm, mwu_job_norm[[i]]$result$pval)
}

colnames(mwu_results_norm) <- c("term", "nseqs", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

mwu_results_norm$term_single <- sapply(mwu_results_norm$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])

mwu_trans_norm <- as.data.frame(cbind(mwu_job_norm[[1]]$seq$id, as.numeric(mwu_job_norm[[1]]$seq$stat), mwu_job_norm[[1]]$seq$random_fc)) 
for(i in 2:length(mwu_job_norm)){
  mwu_trans_norm <- cbind(mwu_trans_norm, mwu_job_norm[[i]]$seq$random_fc)
}

colnames(mwu_trans_norm) <- c("transcript", "original FC", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

# MWU norm 2

mwu_results_norm_2 <- as.data.frame(cbind(mwu_job_norm_2[[1]]$result$term, as.numeric(mwu_job_norm_2[[1]]$result$nseqs), mwu_job_norm_2[[1]]$result$pval)) 
for(i in 2:length(mwu_job_norm_2)){
  mwu_results_norm_2 <- cbind(mwu_results_norm_2, mwu_job_norm_2[[i]]$result$pval)
}

colnames(mwu_results_norm_2) <- c("term", "nseqs", unlist(lapply(1:1000,function(x) paste("result", x, sep=""))))

mwu_results_norm_2$term_single <- sapply(mwu_results_norm_2$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])

mwu_trans_norm_2 <- as.data.frame(cbind(mwu_job_norm_2[[1]]$seq$id, as.numeric(mwu_job_norm_2[[1]]$seq$stat), mwu_job_norm_2[[1]]$seq$stat)) 
for(i in 2:length(mwu_job_norm_2)){
  mwu_trans_norm_2 <- cbind(mwu_trans_norm_2, mwu_job_norm_2[[i]]$seq$stat)
}

colnames(mwu_trans_norm_2) <- c("transcript", "original FC", unlist(lapply(1:1000,function(x) paste("result", x, sep=""))))

# correction
mwu_trans_norm_2$`original FC` <- mwu_trans_norm$`original FC`

### prepare go stat runs 

# function to extract the rest

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

# extract first result - reshuffle
N <- length(ks_job[[1]]$result[[1]]@score)

ks_results_classic <- GenTable(GOdat, ks_job[[1]][["result"]][[1]], topNodes=N) %>% arrange(GO.ID)

ks_results_elim <- GenTable(GOdat, ks_job[[1]][["result"]][[3]], topNodes=N) %>% arrange(GO.ID)

ks_results_weight01 <- GenTable(GOdat, ks_job[[1]][["result"]][[5]], topNodes=N) %>% arrange(GO.ID)

ks_results_lea <- GenTable(GOdat, ks_job[[1]][["result"]][[7]], topNodes=N) %>% arrange(GO.ID)

# extract first result - norm
N <- length(ks_job_norm[[1]]$result[[1]]@score)

ks_results_classic_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[1]], topNodes=N) %>% arrange(GO.ID)

ks_results_elim_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[3]], topNodes=N) %>% arrange(GO.ID)

ks_results_weight01_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[5]], topNodes=N) %>% arrange(GO.ID)

ks_results_lea_norm <- GenTable(GOdat, ks_job_norm[[1]][["result"]][[7]], topNodes=N) %>% arrange(GO.ID)

# extract first result - norm2 
N <- length(ks_job_norm_2[[1]]$result[[1]]@score)

ks_results_classic_norm_2 <- GenTable(GOdat, ks_job_norm_2[[1]][["result"]][[1]], topNodes=N) %>% arrange(GO.ID)

ks_results_elim_norm_2 <- GenTable(GOdat, ks_job_norm_2[[1]][["result"]][[3]], topNodes=N) %>% arrange(GO.ID)

ks_results_weight01_norm_2 <- GenTable(GOdat, ks_job_norm_2[[1]][["result"]][[5]], topNodes=N) %>% arrange(GO.ID)

ks_results_lea_norm_2 <- GenTable(GOdat, ks_job_norm_2[[1]][["result"]][[7]], topNodes=N) %>% arrange(GO.ID)

# get the rest

ks_results_classic <- get_all_its(ks_results_classic, ks_job, 1)
ks_results_elim <- get_all_its(ks_results_elim, ks_job, 3)
ks_results_weight01 <- get_all_its(ks_results_weight01, ks_job, 5)
ks_results_lea <- get_all_its(ks_results_lea, ks_job, 7)
ks_trans <- as.data.frame(cbind(names(ks_job[[1]]$seq),
                                ks_job[[1]]$seq))
for(i in 2:length(ks_job)){
  ks_trans <- cbind(ks_trans, ks_job[[i]]$seq)
}
colnames(ks_trans) <- c("transcript", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))


ks_results_classic_norm <- get_all_its(ks_results_classic_norm, ks_job_norm, 1)
ks_results_elim_norm <- get_all_its(ks_results_elim_norm, ks_job_norm, 3)
ks_results_weight01_norm <- get_all_its(ks_results_weight01_norm, ks_job_norm, 5)
ks_results_lea_norm <- get_all_its(ks_results_lea_norm, ks_job_norm, 7)
ks_trans_norm <- as.data.frame(cbind(names(ks_job_norm[[1]]$seq),
                                     ks_job_norm[[1]]$seq))

for(i in 2:length(ks_job_norm)){
  ks_trans_norm <- cbind(ks_trans_norm, ks_job_norm[[i]]$seq)
}
colnames(ks_trans_norm) <- c("transcript", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

# get everything form norm2 run

ks_results_classic_norm_2 <- get_all_its(ks_results_classic_norm_2, ks_job_norm_2, 1)
ks_results_elim_norm_2 <- get_all_its(ks_results_elim_norm_2, ks_job_norm_2, 3)
ks_results_weight01_norm_2 <- get_all_its(ks_results_weight01_norm_2, ks_job_norm_2, 5)
ks_results_lea_norm_2 <- get_all_its(ks_results_lea_norm_2, ks_job_norm_2, 7)
ks_trans_norm_2 <- as.data.frame(cbind(names(ks_job_norm_2[[1]]$seq),
                                     ks_job_norm_2[[1]]$seq))

for(i in 2:length(ks_job_norm_2)){
  ks_trans_norm_2 <- cbind(ks_trans_norm_2, ks_job_norm_2[[i]]$seq)
}
colnames(ks_trans_norm_2) <- c("transcript", unlist(lapply(1:B,function(x) paste("result", x, sep=""))))

# get run infos

ks_run_metadata_norm_2 <- as.data.frame(cbind(run = c(1:1000),
                                              seed = sapply(ks_job_norm_2, function(x) x$seed$seed),
                                mean = sapply(ks_job_norm_2, function(x) x$seed$mean),
                                percentage = sapply(ks_job_norm_2, function(x) x$seed$percentage)))

mwu_run_metadata_norm_2 <- as.data.frame(cbind(run = c(1:1000),
                                               seed = sapply(mwu_job_norm_2, function(x) x$seed$seed),
                                               mean = sapply(mwu_job_norm_2, function(x) x$seed$mean),
                                               percentage = sapply(mwu_job_norm_2, function(x) x$seed$percentage)))


######### get results from single runs

ks_results <- GenTable(GOdat, classic=ksrun[[1]], elim = ksrun[[3]], weight01 = ksrun[[5]], lea = ksrun[[7]], topNodes = 3191)

ks_results <- as.data.frame(type.convert(ks_results, as.is=T)) %>% arrange(GO.ID)





######################## PADJUST FOR CLASSIC and MWU
# FDR and Bonferroni

############
correct <- "bonferroni"

correct.it <- function(method, pvals, cols){
  for(i in cols){
    pvals[,i] <- p.adjust(pvals[,i], method=method)
  }
  return(pvals)
  }

mwu_results_rand_bonf <- correct.it("fdr", mwu_results_rand, 3:103)
mwu_results_norm_bonf <- correct.it("fdr", mwu_results_norm, 3:103)
mwu_results_norm_2_bonf <- correct.it("fdr", mwu_results_norm_2, 3:103)

ks_results_classic_bh
ks_results_classic_norm_bh
ks_results_classic_norm_2_bh







########### function how many times goterm is over threshold




threshold <- 0.05

over_threshold <- function(data, pvalcols, threshold = 0.05){
  return(apply(data[,pvalcols], 1, function(x) length(which(x < threshold))))
}

############# Infos about GO terms

## GET infos about GO TERMS
transinfo <- BP_clust[!duplicated(BP_clust$seq),]

GO_info <- data.frame(term = names(annotsrev))
GO_info$ntrans_overall <- sapply(annotsrev, length) # number of members in the dataset
GO_info$ntrans_real <- sapply(annotsrev, function(x) length(which(x %in% transinfo$seq))) # number of members in the dataset
GO_info$meanFC <- sapply(annotsrev, function(x) mean(transinfo$FC.sep[which(transinfo$seq %in% x)], na.rm = T)) # mean fold change of all members in the go term
GO_info$medianFC <- sapply(annotsrev, function(x) median(transinfo$FC.sep[which(transinfo$seq %in% x)], na.rm = T)) # median fold change of all members in the go term
GO_info$n_up <- sapply(annotsrev, function(x) length(transinfo$FC.sep[which(transinfo$seq %in% x & transinfo$FC.sep >= 1)])) # number of members being "upregulated"
GO_info$n_down <- sapply(annotsrev, function(x) length(transinfo$FC.sep[which(transinfo$seq %in% x & transinfo$FC.sep <= -1)])) # number of members being "downregulated"
GO_info$n_interesting <- sapply(annotsrev, function(x) length(transinfo$FC.sep[which(transinfo$seq %in% x & (transinfo$FC.sep >= 1 | transinfo$FC.sep <= -1))])) # number of all interesting members (sum of down and upregulated)

GO_info <- GO_info[which(GO_info$term %in% ks_results$GO.ID),] %>% arrange(ntrans_real)


######### Infos about Transcripts

seq_info <- data.frame(seq=names(annots))
seq_info$nmember_overall <- sapply(annots, length) # how many memberships in the ontology
seq_info$nmember_real <- sapply(annots, function(x) length(which(x %in% GO_info$term))) # how many memberships in the dataset
# seq_info$ncomember <- sapply(annots, function(x) length(unique(unlist(annotsrev[which(names(annotsrev) %in% x)])))) # number of transcripts which share the same go terms on the first level -> needs work because rn all transcripts are comembers
seq_info$mean_member <- sapply(annots, function(x) mean(sapply(annotsrev[which(names(annotsrev) %in% x)], length))) # mean number of members in the go term
seq_info$median_member <- sapply(annots, function(x) median(sapply(annotsrev[which(names(annotsrev) %in% x)], length))) # median number of members in the go terms 

# how many go terms where the transcript has a membership are enriched in each method
seq_info$nmember_classic_sig <- sapply(annots, function(x) length(which(ks_results$classic[which(ks_results$GO.ID %in% x)] < 0.05)))
seq_info$nmember_elim_sig <- sapply(annots, function(x) length(which(ks_results$elim[which(ks_results$GO.ID %in% x)] < 0.05)))
seq_info$nmember_weight_sig <- sapply(annots, function(x) length(which(ks_results$weight01[which(ks_results$GO.ID %in% x)] < 0.05)))
seq_info$nmember_lea_sig <- sapply(annots, function(x) length(which(ks_results$lea[which(ks_results$GO.ID %in% x)] < 0.05)))
seq_info$nmember_mwu <- sapply(annots, function(x) length(which(mwu_run$pval[which(mwu_run$term_single %in% x)] < 0.05)))

seq_info$nmember_classic_sig_fdr <- sapply(annots, function(x) length(which(p.adjust(ks_results$classic[which(ks_results$GO.ID %in% x)], method="fdr") < 0.05)))
seq_info$nmember_elim_sig_fdr <- sapply(annots, function(x) length(which(p.adjust(ks_results$elim[which(ks_results$GO.ID %in% x)], method="fdr") < 0.05)))
seq_info$nmember_weight_sig_fdr <- sapply(annots, function(x) length(which(p.adjust(ks_results$weight01[which(ks_results$GO.ID %in% x)], method="fdr") < 0.05)))
seq_info$nmember_mwu_fdr <- sapply(annots, function(x) length(which(p.adjust(mwu_run$pval[which(mwu_run$term_single %in% x)], method="fdr") < 0.05)))

# seq_info$nmember_wilcox <- sapply(annots, function(x) length(which(wilcox_results$pval[which(wilcox_results$term_single %in% x)] < 0.05)))

seq_info <- left_join(seq_info, transinfo[, c(1,10)], by="seq")

seq_info <- na.omit(seq_info)

# get number of times for mwu

# random resample

mwu_trans_rand_n <- data.frame(transcript=names(annots))
mwu_trans_rand_n <- left_join(mwu_trans_rand_n, mwu_trans_rand[,c(1,2)], by="transcript")

mwu_trans_rand_n$first <- sapply(annots, function(x) length(which(mwu_trans_rand$result1[which(mwu_results_rand$term_single %in% x)] < 0.05)))
for(i in 4:dim(mwu_trans_rand)[2]){
  print(i)
  mwu_trans_rand_n <- cbind(mwu_trans_rand_n, sapply(annots, function(x) length(which(mwu_trans_rand[,i][which(mwu_results_rand$term_single %in% x)] < 0.05))))
}
mwu_trans_rand_n <- na.omit(mwu_trans_rand_n)

# sample from N

mwu_trans_norm_n <- data.frame(transcript=names(annots))
mwu_trans_norm_n <- left_join(mwu_trans_norm_n, mwu_trans_norm[,c(1,2)], by="transcript")

mwu_trans_norm_n$first <- sapply(annots, function(x) length(which(mwu_trans_norm$result1[which(mwu_results_norm$term_single %in% x)] < 0.05)))
for(i in 4:dim(mwu_trans_norm)[2]){
  print(i)
  mwu_trans_norm_n <- cbind(mwu_trans_norm_n, sapply(annots, function(x) length(which(mwu_trans_norm[,i][which(mwu_results_norm$term_single %in% x)] < 0.05))))
}
mwu_trans_norm_n <- na.omit(mwu_trans_norm_n)

# sample from norm 2


mwu_trans_norm_n_2 <- data.frame(transcript=names(annots))
mwu_trans_norm_n_2 <- left_join(mwu_trans_norm_n_2, mwu_trans_norm_2[,c(1,2)], by="transcript")

mwu_trans_norm_n_2$first <- sapply(annots, function(x) length(which(mwu_trans_norm_2$result1[which(mwu_results_norm_2$term_single %in% x)] < 0.05)))
for(i in 4:dim(mwu_trans_norm_2)[2]){
  print(i)
  mwu_trans_norm_n_2 <- cbind(mwu_trans_norm_n_2, sapply(annots, function(x) length(which(mwu_trans_norm_2[,i][which(mwu_results_norm_2$term_single %in% x)] < 0.05))))
}
mwu_trans_norm_n_2 <- na.omit(mwu_trans_norm_n_2)


################### prepare go packages
get_trans_go <- function(results, iterations, algo){
  res <- data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq)
  for(i in 2:iterations){
    res <- cbind(res, results[[i]]$seq)
  }
  res_n <- data.frame(transcript=names(annots))
  res_n <- left_join(res_n, data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq), by="transcript")
  # classic
  res_classic_n <- res_n
  res_classic_n$first <- sapply(annots, function(x) length(which(ks_results_classic$result1[which(ks_results_classic$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_classic_n <- cbind(res_classic_n, sapply(annots, function(x) length(which(ks_results_classic[,i][which(ks_results_classic$GO.ID %in% x)] < 0.05))))
  }
  res_classic_n <- na.omit(res_classic_n)
  
  # elim 
  res_elim_n <- res_n
  res_elim_n$first <- sapply(annots, function(x) length(which(ks_results_elim$result1[which(ks_results_elim$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_elim_n <- cbind(res_elim_n, sapply(annots, function(x) length(which(ks_results_elim[,i][which(ks_results_elim$GO.ID %in% x)] < 0.05))))
  }
  res_elim_n <- na.omit(res_elim_n)
  
  # weight
  res_weight_n <- res_n
  res_weight_n$first <- sapply(annots, function(x) length(which(ks_results_weight01$result1[which(ks_results_weight01$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_weight_n <- cbind(res_weight_n, sapply(annots, function(x) length(which(ks_results_weight01[,i][which(ks_results_weight01$GO.ID %in% x)] < 0.05))))
  }
  res_weight_n <- na.omit(res_weight_n)
  
  # lea
  res_lea_n <- res_n
  res_lea_n$first <- sapply(annots, function(x) length(which(ks_results_lea$result1[which(ks_results_lea$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_lea_n <- cbind(res_lea_n, sapply(annots, function(x) length(which(ks_results_lea[,i][which(ks_results_lea$GO.ID %in% x)] < 0.05))))
  }
  res_lea_n <- na.omit(res_lea_n)
  
  res <- res[which(res$transcript %in% res_classic_n$transcript),]
  
  return(list(res, res_classic_n, res_elim_n, res_weight_n, res_lea_n))
}

res_rand_n <- get_trans_go(ks_job, 100, )

get_norm_go <- function(results, iterations, algo){
  res <- data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq)
  for(i in 2:iterations){
    res <- cbind(res, results[[i]]$seq)
  }
  res_n <- data.frame(transcript=names(annots))
  res_n <- left_join(res_n, data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq), by="transcript")
  # classic
  res_classic_n <- res_n
  res_classic_n$first <- sapply(annots, function(x) length(which(ks_results_classic_norm$result1[which(ks_results_classic_norm$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_classic_n <- cbind(res_classic_n, sapply(annots, function(x) length(which(ks_results_classic_norm[,i][which(ks_results_classic_norm$GO.ID %in% x)] < 0.05))))
  }
  res_classic_n <- na.omit(res_classic_n)
  
  # elim 
  res_elim_n <- res_n
  res_elim_n$first <- sapply(annots, function(x) length(which(ks_results_elim_norm$result1[which(ks_results_elim_norm$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_elim_n <- cbind(res_elim_n, sapply(annots, function(x) length(which(ks_results_elim_norm[,i][which(ks_results_elim_norm$GO.ID %in% x)] < 0.05))))
  }
  res_elim_n <- na.omit(res_elim_n)
  
  # weight
  res_weight_n <- res_n
  res_weight_n$first <- sapply(annots, function(x) length(which(ks_results_weight01_norm$result1[which(ks_results_weight01_norm$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_weight_n <- cbind(res_weight_n, sapply(annots, function(x) length(which(ks_results_weight01_norm[,i][which(ks_results_weight01_norm$GO.ID %in% x)] < 0.05))))
  }
  res_weight_n <- na.omit(res_weight_n)
  
  # lea
  res_lea_n <- res_n
  res_lea_n$first <- sapply(annots, function(x) length(which(ks_results_lea_norm$result1[which(ks_results_lea_norm$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_lea_n <- cbind(res_lea_n, sapply(annots, function(x) length(which(ks_results_lea_norm[,i][which(ks_results_lea_norm$GO.ID %in% x)] < 0.05))))
  }
  res_lea_n <- na.omit(res_lea_n)
  
  res <- res[which(res$transcript %in% res_classic_n$transcript),]
  
  return(list(res, res_classic_n, res_elim_n, res_weight_n, res_lea_n))
}

res_norm_n <- get_norm_go(ks_job_norm, 100, )

## get from ks run norm 2

get_norm_go_2 <- function(results, iterations, algo){
  res <- data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq)
  for(i in 2:iterations){
    res <- cbind(res, results[[i]]$seq)
  }
  res_n <- data.frame(transcript=names(annots))
  res_n <- left_join(res_n, data.frame(transcript = names(results[[1]]$seq), fc = results[[1]]$seq), by="transcript")
  # classic
  res_classic_n <- res_n
  res_classic_n$first <- sapply(annots, function(x) length(which(ks_results_classic_norm_2$result1[which(ks_results_classic_norm_2$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_classic_n <- cbind(res_classic_n, sapply(annots, function(x) length(which(ks_results_classic_norm_2[,i][which(ks_results_classic_norm_2$GO.ID %in% x)] < 0.05))))
  }
  res_classic_n <- na.omit(res_classic_n)
  
  # elim 
  res_elim_n <- res_n
  res_elim_n$first <- sapply(annots, function(x) length(which(ks_results_elim_norm_2$result1[which(ks_results_elim_norm_2$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_elim_n <- cbind(res_elim_n, sapply(annots, function(x) length(which(ks_results_elim_norm_2[,i][which(ks_results_elim_norm_2$GO.ID %in% x)] < 0.05))))
  }
  res_elim_n <- na.omit(res_elim_n)
  
  # weight
  res_weight_n <- res_n
  res_weight_n$first <- sapply(annots, function(x) length(which(ks_results_weight01_norm_2$result1[which(ks_results_weight01_norm_2$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_weight_n <- cbind(res_weight_n, sapply(annots, function(x) length(which(ks_results_weight01_norm_2[,i][which(ks_results_weight01_norm_2$GO.ID %in% x)] < 0.05))))
  }
  res_weight_n <- na.omit(res_weight_n)
  
  # lea
  res_lea_n <- res_n
  res_lea_n$first <- sapply(annots, function(x) length(which(ks_results_lea_norm_2$result1[which(ks_results_lea_norm_2$GO.ID %in% x)] < 0.05)))
  
  for(i in 7:(iterations+6)){
    res_lea_n <- cbind(res_lea_n, sapply(annots, function(x) length(which(ks_results_lea_norm_2[,i][which(ks_results_lea_norm_2$GO.ID %in% x)] < 0.05))))
  }
  res_lea_n <- na.omit(res_lea_n)
  
  res <- res[which(res$transcript %in% res_classic_n$transcript),]
  
  return(list(res, res_classic_n, res_elim_n, res_weight_n, res_lea_n))
}

res_norm_n_2 <- get_norm_go_2(ks_job_norm_2, 1000, )


testruns <- list(ks=list(original=ksrun,
                           reshuffle=ks_job,
                           normal=ks_job_norm,
                           normal2=ks_job_norm_2),
                  mwu=list(original=mwu_run,
                           reshuffle=mwu_job,
                           normal=mwu_job_norm,
                           normal2=mwu_job_norm_2))
save(testruns, file="./dataframes/testruns.RData")

results_reshuffle <- list(original=list(ks=ks_results,
                                        mwu=mwu_run,
                                        seq_info=seq_info),
                  classic=list(pval=ks_results_classic,
                                       seq=ks_trans,
                                       seq_n=res_rand_n[[2]]),
                  
                         weight=list(pval=ks_results_weight01,
                                     seq=ks_trans,
                                     seq_n=res_rand_n[[4]]),
                  
                         elim=list(pval=ks_results_elim,
                                   seq=ks_trans,
                                   seq_n=res_rand_n[[3]]),
                  
                         lea=list(pval=ks_results_lea,
                                  seq=ks_trans,
                                  seq_n=res_rand_n[[5]]),
                        mwu=list(pval=mwu_results_rand,
                                 seq=mwu_trans_rand,
                                 seq_n=mwu_trans_rand_n))

results_normal <- list(original=list(ks=ks_results,
                                        mwu=mwu_run,
                                        seq_info=seq_info),
                          classic=list(pval=ks_results_classic_norm,
                                       seq=ks_trans_norm,
                                       seq_n=res_norm_n[[2]]),
                          
                          weight=list(pval=ks_results_weight01_norm,
                                      seq=ks_trans_norm,
                                      seq_n=res_norm_n[[4]]),
                          
                          elim=list(pval=ks_results_elim_norm,
                                    seq=ks_trans_norm,
                                    seq_n=res_norm_n[[3]]),
                          
                          lea=list(pval=ks_results_lea_norm,
                                   seq=ks_trans_norm,
                                   seq_n=res_norm_n[[5]]),
                       
                          mwu=list(pval=mwu_results_norm,
                                   seq=mwu_trans_norm,
                                   seq_n=mwu_trans_norm_n))

results_normal2 <- list(original=list(ks=ks_results,
                                     mwu=mwu_run,
                                     seq_info=seq_info),
                       classic=list(pval=ks_results_classic_norm_2,
                                    seq=ks_trans_norm_2,
                                    seq_n=res_norm_n_2[[2]]),
                       
                       weight=list(pval=ks_results_weight01_norm_2,
                                   seq=ks_trans_norm_2,
                                   seq_n=res_norm_n_2[[4]]),
                       
                       elim=list(pval=ks_results_elim_norm_2,
                                 seq=ks_trans_norm_2,
                                 seq_n=res_norm_n_2[[3]]),
                       
                       lea=list(pval=ks_results_lea_norm_2,
                                seq=ks_trans_norm_2,
                                seq_n=res_norm_n_2[[5]]),
                       
                       mwu=list(pval=mwu_results_norm_2,
                                seq=mwu_trans_norm_2,
                                seq_n=mwu_trans_norm_n_2),
                       metadata=list(ks=ks_run_metadata_norm_2,
                                     mwu=mwu_run_metadata_norm_2))

misc <- list(GO=GO_info,
             trans=seq_info,
             annots=annots,
             annotsrev=annotsrev,
             origannots=annotations,
             originaldata=BP_clust)

save(misc, file="./dataframes/misc.RData")
save(results_reshuffle, file="./dataframes/results_reshuffle.RData")
save(results_normal, file="./dataframes/results_normal.RData")
save(results_normal2, file="./dataframes/results_normal2.RData")
save(misc, results_normal, results_normal2, results_reshuffle, file="./dataframes/all_working_files.RData")
save.image(file="./20230707")


results_orig <-  list(ks=ks_results,
              mwu=mwu_run,
              seq_info=seq_info,
              go =GO_info,
              medi_lines = orig_fc_N)

save(results_orig, file="./data/orig_results.RData")
