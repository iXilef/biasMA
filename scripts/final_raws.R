library(gageData)
library(topGO)
library(foreach)
library(dplyr)
# library(doParallel)
# cores=detectCores()
# cl <- makeCluster(cores[1]-2)
# registerDoParallel(cl)


data(go.sets.hs)
data(go.subs.hs)

human <- go.sets.hs[go.subs.hs$BP]
names(human) <- sapply(names(human), function(x) unlist(strsplit(x, " ", fixed=T))[1])

# delete GO sets below 10
# new_human <- human[which(sapply(human, function(x) length(x) > 9))]
trans_hum <- inverseList(human)
trans_hum <- sample(trans_hum, 1000)
go_human <- inverseList(trans_hum) 
go_human <- go_human[which(sapply(go_human, function(x) length(x) > 9 ))]

results <- list()

B <- 100
res <- list()

for(b in 1:B){
  print(b)
  results <- list()
    for(i in seq(.0, 0.35, 0.1)){
      print(i)
      set.seed(runif(1, 0, 10000000))
      intercept <- 1.24221*i
      int.data <- data.frame(trans = names(trans_hum))
      int.data$gomem <- sapply(trans_hum, length)
      int.data$FC <- rnorm(1000, mean=i*log10(int.data$gomem)-intercept, sd=0.9192379)
      
      mwu.p <- 1:1063
      ks.p <- 1:1063
       
      mwu.p_dupl <- 1:1063
      ks.p_dupl <- 1:1063
      bg2 <- rep(int.data$FC, sapply(trans_hum, length))
      
      for(j in 1:length(go_human)){
        # without duplication of transcript in background
        int.data$FC[which(go_human[[j]] %in% int.data$trans)]
        
      
        mwu.p[j] <- wilcox.test(int.data$FC[which(go_human[[j]] %in% int.data$trans)], int.data$FC)$p.value
        ks.p[j] <- ks.test(int.data$FC[which(go_human[[j]] %in% int.data$trans)], int.data$FC)$p.value
        
        mwu.p_dupl[j] <- wilcox.test(int.data$FC[which(go_human[[j]] %in% int.data$trans)], bg2)$p.value
        ks.p_dupl[j] <- ks.test(int.data$FC[which(go_human[[j]] %in% int.data$trans)], bg2)$p.value
        
      }
      pvals <- cbind(mwu.p = mwu.p, ks.p = ks.p, mwu.p_dupl = mwu.p_dupl, ks.p_dupl = ks.p_dupl)
      results <- append(results, list(pvals))
    }
    
    counts <- c()
    
    for(i in 1:4){
      for(j in 1:4){
        counts <- c(counts, unlist(length(which(results[[i]][,j] < 0.05))))
      }
    }
    counts <- unlist(counts)
# 
res <- append(res, results)
}

ergebnis <- NULL
for(i in 1:length(res)){
  ergebnis <- rbind(ergebnis, apply(res[[i]], 2, function(x) length(which(x < 0.05))))
}
