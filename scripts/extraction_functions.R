# 

get.data <- function(data=F, sequ=F, what=F, correct="none", iterations=F, GoDat=F, annotfile=F){
  # extract transcripts
  if(what=="MWU"){
    fin.data <- setNames(rep(list(NaN), 8), sequ)
    
    for(slope in 1:length(sequ)){
      
      int.data <- data[ifelse(slope==1, slope, (iterations*(slope-1))+1):(slope*iterations)]
      annot.int <- annotfile[names(annotfile) %in% unique(int.data[[1]]$seq$seq)]
      transcripts <- data.frame(matrix(nrow = length(unique(int.data[[1]]$seq$seq)), ncol = iterations), row.names = unique(int.data[[1]]$seq$seq)) ## DEFINE LENGTH
      GOterms <- data.frame(matrix(nrow = length(int.data[[1]]$result$term), ncol = iterations), 
                            row.names = sapply(int.data[[1]]$result$term, function(x) strsplit(x, ";", fixed=T)[[1]][1]))
      transcripts_n <- data.frame(matrix(nrow = length(unique(int.data[[1]]$seq$seq)), ncol = iterations), row.names = unique(int.data[[1]]$seq$seq))
      
      for(iteration in 1:iterations){
        print(paste("iteration:", iteration, "in slope: ", sequ[slope]))
        transcripts[,iteration] <- int.data[[iteration]]$seq$value[!duplicated(int.data[[iteration]]$seq$seq)] # DEFINE things
        
        GOterms[,iteration] <- p.adjust(int.data[[iteration]]$result$pval, method=correct)
        
        transcripts_n[,iteration] <- sapply(annot.int, function(x) length(which(GOterms[row.names(GOterms) %in% x ,iteration] < .05)))
      }
      
      fin.data[[slope]] <- list(GOterms = GOterms, transcripts = transcripts, transcripts_n)
    }
  } else {
    
    
    fin.data <- setNames(rep(list(NA), 5), c("classic", "elim", "weight", "lea", "misc"))
    fin.data$classic <- setNames(rep(list(NA), 8), sequ)
    fin.data$elim <- setNames(rep(list(NA), 8), sequ)
    fin.data$weight <- setNames(rep(list(NA), 8), sequ)
    fin.data$lea <- setNames(rep(list(NA), 8), sequ)
    
    for(slope in 1:length(sequ)){
      
      print(paste("slope: ", sequ[slope]))
      
      int.data <- data[ifelse(slope==1, slope, (iterations*(slope-1))+1):(slope*iterations)]
      annot.int <- annotfile[names(annotfile) %in% names(int.data[[1]]$seq)]
      all(names(annot.int)==names(int.data[[1]]$seq))
      
      ### CLASSIC
      
      transcripts <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq)) ## DEFINE LENGTH
      GOterms_c <- data.frame(matrix(nrow = length(int.data[[1]]$result[[1]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_c_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[1]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_c <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_c_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      ### ELIM
      GOterms_e <- data.frame(matrix(nrow = length(int.data[[1]]$result[[3]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_e_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[3]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_e <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_e_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      ### WEIGHT
      GOterms_w <- data.frame(matrix(nrow = length(int.data[[1]]$result[[5]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_w_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[5]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_w <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_w_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
            
      ### LEA
      GOterms_l <- data.frame(matrix(nrow = length(int.data[[1]]$result[[7]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_l_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[7]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_l <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_l_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      for(iteration in 1:iterations){
        print(paste("iteration:", iteration, "in slope: ", sequ[slope]))
        
        transcripts[,iteration] <- int.data[[iteration]]$seq # DEFINE things
        
        ### CLASSIC

        GOterms_c[,iteration] <- p.adjust(int.data[[iteration]]$result[[1]]@score, method="none")
        
        GOterms_c_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[1]]@score, method="fdr")

        transcripts_n_c[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_c[row.names(GOterms_c) %in% x ,iteration] < .05)))
        
        transcripts_n_c_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_c_fdr[row.names(GOterms_c_fdr) %in% x ,iteration] < .05)))
        

        ### ELIM

        GOterms_e[,iteration] <-  p.adjust(int.data[[iteration]]$result[[3]]@score, method="none")
        
        GOterms_e_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[3]]@score, method="fdr")

        transcripts_n_e[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_e[row.names(GOterms_e) %in% x ,iteration] < .05)))
        
        transcripts_n_e_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_e_fdr[row.names(GOterms_e_fdr) %in% x ,iteration] < .05)))
        

        ### WEIGHT

        GOterms_w[,iteration] <-  p.adjust(int.data[[iteration]]$result[[5]]@score, method="none")
        
        GOterms_w_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[5]]@score, method="fdr")
        
        transcripts_n_w[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_w[row.names(GOterms_w) %in% x ,iteration] < .05)))
                transcripts_n_w_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_w_fdr[row.names(GOterms_w_fdr) %in% x ,iteration] < .05)))
        

        ### LEA

        GOterms_l[,iteration] <-  p.adjust(int.data[[iteration]]$result[[7]]@score, method="none")
        
        GOterms_l_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[7]]@score, method="fdr")

        transcripts_n_l[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_l[row.names(GOterms_l) %in% x ,iteration] < .05)))
        
        transcripts_n_l_fdr <- sapply(annot.int, function(x) length(which(GOterms_l_fdr[row.names(GOterms_l_fdr) %in% x ,iteration] < .05)))
      }
      
      trans <- data.frame(trans = names(int.data[[1]]$seq),
                          GO_N=sapply(annot.int, function(x) length(which(names(int.data[[1]]$result[[1]]@score) %in% x))))
      
      annot.int.rev <- inverseList(annot.int)
      # prune list to present GO terms 
      annot.int.rev <- annot.int.rev[names(annot.int.rev) %in% names(int.data[[1]]$result[[1]]@score)]
      
      annotfile[names(annotfile) %in% names(int.data[[1]]$seq)]
      
      go <- GenTable(auxdata$rat$GODat, int.data[[1]]$result[[1]], topNodes=length(int.data[[1]]$result[[1]]@score))[,c(1,3)] %>% arrange(GO.ID)
      
      fin.data$classic[[slope]] <- list(GOterms = GOterms_c, trans = transcripts, trans_n = transcripts_n_c, trans_n_fdr = transcripts_n_c_fdr)
      fin.data$elim[[slope]] <- list(GOterms = GOterms_e, trans = transcripts, trans_n = transcripts_n_e, trans_n_fdr = transcripts_n_e_fdr)
      fin.data$weight[[slope]] <- list(GOterms = GOterms_w, trans = transcripts, trans_n = transcripts_n_w, trans_n_fdr = transcripts_n_w_fdr)
      fin.data$lea[[slope]] <- list(GOterms = GOterms_l, trans = transcripts, trans_n = transcripts_n_l, trans_n_fdr = transcripts_n_l_fdr)
      fin.data$misc <- list(GOterms = go, transcripts = trans)
      
      
    }

  }
  # extract GOterm pvals
  # extract transcript N times sig
  # extract lines for plots
  return(fin.data)
}


get.data.parallel <- function(data=F, sequ=F, what=F, correct="none", iterations=F, GoDat=F, annotfile=F){
  library(foreach)
  library(doParallel)
  library(topGO)
  library(dplyr)
  cores=detectCores()
  cl <- makeCluster(cores[1]-2)
  registerDoParallel(cl)
  

  
  # extract transcripts
  if(what=="MWU"){
    
    resultandseed <- function(result=NULL,seed=NULL)
    {
      me <- list(
        mwu = mwu
      )
      
      # Set the name for the class
      class(me) <- append(class(me),"resultandseed")
      return(me)
    }
    
    fin.data <- foreach(slope = 1:length(sequ)) %dopar% {
      int.data <- data[ifelse(slope==1, slope, (iterations*(slope-1))+1):(slope*iterations)]
      annot.int <- annotfile[names(annotfile) %in% unique(int.data[[1]]$seq$seq)]
      transcripts <- data.frame(matrix(nrow = length(unique(int.data[[1]]$seq$seq)), ncol = iterations), row.names = unique(int.data[[1]]$seq$seq)) ## DEFINE LENGTH
      GOterms <- data.frame(matrix(nrow = length(int.data[[1]]$result$term), ncol = iterations), 
                            row.names = sapply(int.data[[1]]$result$term, function(x) strsplit(x, ";", fixed=T)[[1]][1]))
      
      GOterms_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result$term), ncol = iterations), 
                            row.names = sapply(int.data[[1]]$result$term, function(x) strsplit(x, ";", fixed=T)[[1]][1]))
      
      transcripts_n <- data.frame(matrix(nrow = length(unique(int.data[[1]]$seq$seq)), ncol = iterations), row.names = unique(int.data[[1]]$seq$seq))
      
      transcripts_n_fdr <- data.frame(matrix(nrow = length(unique(int.data[[1]]$seq$seq)), ncol = iterations), row.names = unique(int.data[[1]]$seq$seq))
      
      for(iteration in 1:iterations){
        print(iteration)
        transcripts[,iteration] <- int.data[[iteration]]$seq$value[!duplicated(int.data[[iteration]]$seq$seq)] # DEFINE things
        
        GOterms[,iteration] <- p.adjust(int.data[[iteration]]$result$pval, method="none")
        
        GOterms_fdr[,iteration] <- p.adjust(int.data[[iteration]]$result$pval, method="fdr")
        
        transcripts_n[,iteration] <- sapply(annot.int, function(x) length(which(GOterms[row.names(GOterms) %in% x ,iteration] < .05)))
        
        transcripts_n_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_fdr[row.names(GOterms_fdr) %in% x ,iteration] < .05)))
      }
      
      #################################
      
      trans <- data.frame(table(int.data[[1]]$seq$seq))
      colnames(trans) <- c("trans", "GO_N")
      trans <- trans[order(match(trans$trans,row.names(transcripts))),]

      go <- data.frame(table(int.data[[1]]$seq$term))
      colnames(go) <- c("GO", "trans_N")
      go$GO <- sapply(go$GO, function(x) strsplit(as.character(x), ";", fixed=T)[[1]][1])
      go <- go[order(match(go$GO, row.names(GOterms))),]
      
      #################################
      
      # fin.data[[slope]] <- list(GOterms = GOterms, transcripts = transcripts, transcripts_n)
      results <- list(GOterms = GOterms, trans = transcripts, trans_n = transcripts_n, trans_n_fdr = transcripts_n_fdr, misc=list(slope=sequ[slope], go=go, trans=trans))
      return(results)
    }
  } else {
    
    resultandseed <- function(result=NULL,seed=NULL)
    {
      me <- list(
        classic = classic,
        elim = elim,
        weight=weight,
        lea=lea
      )
      
      # Set the name for the class
      class(me) <- append(class(me),"resultandseed")
      return(me)
    }
    
    # 
    # fin.data <- setNames(rep(list(NA), 4), c("classic", "elim", "weight", "lea"))
    # fin.data$classic <- setNames(rep(list(NA), 8), sequ)
    # fin.data$elim <- setNames(rep(list(NA), 8), sequ)
    # fin.data$weight <- setNames(rep(list(NA), 8), sequ)
    # fin.data$lea <- setNames(rep(list(NA), 8), sequ)
    
    fin.data <- foreach(slope = 1:length(sequ), .packages = c("topGO", "dplyr")) %dopar%  {
      
      print(paste("slope: ", sequ[slope]))
      
      int.data <- data[ifelse(slope==1, slope, (iterations*(slope-1))+1):(slope*iterations)]
      annot.int <- annotfile[names(annotfile) %in% names(int.data[[1]]$seq)]
      all(names(annot.int)==names(int.data[[1]]$seq))
      
      ### CLASSIC
      
      transcripts <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq)) ## DEFINE LENGTH
      GOterms_c <- data.frame(matrix(nrow = length(int.data[[1]]$result[[1]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_c_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[1]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_c <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_c_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      ### ELIM
      GOterms_e <- data.frame(matrix(nrow = length(int.data[[1]]$result[[3]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_e_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[3]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_e <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_e_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      ### WEIGHT
      GOterms_w <- data.frame(matrix(nrow = length(int.data[[1]]$result[[5]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_w_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[5]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_w <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_w_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      ### LEA
      GOterms_l <- data.frame(matrix(nrow = length(int.data[[1]]$result[[7]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      GOterms_l_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$result[[7]]@score), ncol = iterations), row.names = names(int.data[[1]]$result[[1]]@score))
      transcripts_n_l <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      transcripts_n_l_fdr <- data.frame(matrix(nrow = length(int.data[[1]]$seq), ncol = iterations), row.names = names(int.data[[1]]$seq))
      
      
      for(iteration in 1:iterations){
        print(paste("iteration:", iteration, "in slope: ", sequ[slope]))
        
        transcripts[,iteration] <- int.data[[iteration]]$seq # DEFINE things
        
        ### CLASSIC
        
        GOterms_c[,iteration] <- p.adjust(int.data[[iteration]]$result[[1]]@score, method="none")
        
        GOterms_c_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[1]]@score, method="fdr")
        
        transcripts_n_c[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_c[row.names(GOterms_c) %in% x ,iteration] < .05)))
        
        transcripts_n_c_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_c_fdr[row.names(GOterms_c_fdr) %in% x ,iteration] < .05)))
        
        
        ### ELIM
        
        GOterms_e[,iteration] <-  p.adjust(int.data[[iteration]]$result[[3]]@score, method="none")
        
        GOterms_e_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[3]]@score, method="fdr")
        
        transcripts_n_e[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_e[row.names(GOterms_e) %in% x ,iteration] < .05)))
        
        transcripts_n_e_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_e_fdr[row.names(GOterms_e_fdr) %in% x ,iteration] < .05)))
        
        
        ### WEIGHT
        
        GOterms_w[,iteration] <-  p.adjust(int.data[[iteration]]$result[[5]]@score, method="none")
        
        GOterms_w_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[5]]@score, method="fdr")
        
        transcripts_n_w[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_w[row.names(GOterms_w) %in% x ,iteration] < .05)))
        
        transcripts_n_w_fdr[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_w_fdr[row.names(GOterms_w_fdr) %in% x ,iteration] < .05)))
        
        
        ### LEA
        
        GOterms_l[,iteration] <-  p.adjust(int.data[[iteration]]$result[[7]]@score, method="none")
        
        GOterms_l_fdr[,iteration] <-  p.adjust(int.data[[iteration]]$result[[7]]@score, method="fdr")
        
        transcripts_n_l[,iteration] <- sapply(annot.int, function(x) length(which(GOterms_l[row.names(GOterms_l) %in% x ,iteration] < .05)))
        
        transcripts_n_l_fdr <- sapply(annot.int, function(x) length(which(GOterms_l_fdr[row.names(GOterms_l_fdr) %in% x ,iteration] < .05)))
      }
      
      trans <- data.frame(trans = names(int.data[[1]]$seq),
                          GO_N=sapply(annot.int, function(x) length(which(names(int.data[[1]]$result[[1]]@score) %in% x))))
      
      annot.int.rev <- inverseList(annot.int)
      # prune list to present GO terms 
      annot.int.rev <- annot.int.rev[names(annot.int.rev) %in% names(int.data[[1]]$result[[1]]@score)]
      
      annotfile[names(annotfile) %in% names(int.data[[1]]$seq)]
      
      go <- GenTable(GoDat, int.data[[1]]$result[[1]], topNodes=length(int.data[[1]]$result[[1]]@score))[,c(1,3)] %>% arrange(GO.ID)
      # 
      # fin.data$classic[[slope]] <- list(GOterms = GOterms_c, trans = transcripts, trans_n = transcripts_n_c)
      # fin.data$elim[[slope]] <- list(GOterms = GOterms_e, trans = transcripts, trans_n = transcripts_n_e)
      # fin.data$weight[[slope]] <- list(GOterms = GOterms_w, trans = transcripts, trans_n = transcripts_n_w)
      # fin.data$lea[[slope]] <- list(GOterms = GOterms_l, trans = transcripts, trans_n = transcripts_n_l)
      # 
      
      results <- list(classic = list(GOterms = GOterms_c, trans = transcripts, trans_n = transcripts_n_c, trans_n_fdr = transcripts_n_c_fdr),
                      elim = list(GOterms = GOterms_e, trans = transcripts, trans_n = transcripts_n_e, trans_n_fdr = transcripts_n_e_fdr),
                      weight = list(GOterms = GOterms_w, trans = transcripts, trans_n = transcripts_n_w, trans_n_fdr = transcripts_n_w_fdr),
                      lea = list(GOterms = GOterms_l, trans = transcripts, trans_n = transcripts_n_l, trans_n_fdr = transcripts_n_l_fdr),
                      misc = list(slope =sequ[slope], go=go, trans=trans))
                      
      # results$classic <- list(GOterms = GOterms_c, trans = transcripts, trans_n = transcripts_n_c)
      # results$elim <- list(GOterms = GOterms_e, trans = transcripts, trans_n = transcripts_n_e)
      # results$weight <- list(GOterms = GOterms_w, trans = transcripts, trans_n = transcripts_n_w)
      # results$lea <- list(GOterms = GOterms_l, trans = transcripts, trans_n = transcripts_n_l)
      
      return(results)
      
      
    }
    
  }
  # extract GOterm pvals
  # extract transcript N times sig
  # extract lines for plots
  return(fin.data)
}

sapply(annot.int, function(x) length(which(transcripts[, iterations][which(names(GOterms_c) %in% x)])))

cbind(mwu_trans_rand_n, sapply(annots, function(x) length(which(mwu_trans_rand[,i][which(mwu_results_rand$term_single %in% x)] < 0.05))))


term_single <- sapply(mwu_results_norm_2$term, function(x) strsplit(x, ";", fixed=T)[[1]][1])