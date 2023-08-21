loess50 <- loess(FC.sep~log10(nmember_real+0.001), data=na.omit(results_orig$seq_info), span=.75)
smooth50 <- predict(loess50)
this <- as.data.frame(cbind(smooth50,fc=results_orig$seq_info$nmember_real))
this <- this %>% arrange(fc)
plot(results_orig$seq_info$FC.sep~log10(results_orig$seq_info$nmember_real))
lines(this$smooth50, x=log10(this$fc), col="red")




fc_N <- as.data.frame(cbind(0:max(results_orig$seq_info$nmember_real), sapply(0:max(results_orig$seq_info$nmember_real),
                                                                              function(x) median(results_orig$seq_info$FC.sep[which(results_orig$seq_info$nmember_real > x)])))); colnames(fc_N) <- c("N", "medianFC")

slope <- seq(0,0.35,0.05)
buchstab <- c("A", "B", "C", "D", "E", "F", "G", "H")

par(mfrow=c(2,4))
for(j in 1:8){
  plot(x=1:500, y=sapply(1:500, function(y) median(orig_results[[j]]$classic$trans[,1][which(x$n >= y)], na.rm=T)), type="l", ylim=c(-1, 1), xlab="Number of annotated GO terms", ylab="Cumulative Median log2FC", main=paste0(buchstab[j], ") ", "slope: ", slope[j]))
  for(i in 2:100){
    lines(x=1:500, y=sapply(1:500, function(y) median(orig_results[[j]]$classic$trans[,i][which(x$n >= y)], na.rm=T)))
  }
  if(j == 6){
    lines(fc_N$N, fc_N$medianFC, col="blue") 
  }
  lines(x=1:500, sapply(1:500, function(y) median(unlist(orig_results[[j]]$classic$trans[which(x$n >= y),]), na.rm=T)), col="red")}



slope <- seq(0,0.35,0.05)
buchstab <- c("A", "B", "C", "D", "E", "F", "G", "H")


par(mfrow=c(1,1))
for(j in 8:8){
  lo <- loess(orig_results[[j]]$classic$trans[,1]~log10(x$n))
  lo <- as.data.frame(cbind(n=log10(x$n), fc=predict(lo)))
  lo2 <- lo
  lo <- lo %>% arrange(n)
  
  plot(x=lo$n, y=lo$fc, type="l", ylim=c(-1, 1),
       xlab="Number of annotated GO terms", ylab="Cumulative Median log2FC", main=paste0(buchstab[j], ") ", "slope: ", slope[j]))
  
  for(i in 2:100){
    print(i)
    lo <- loess(orig_results[[j]]$classic$trans[,i]~log10(x$n))
    lo <- as.data.frame(cbind(fc=predict(lo), n=log10(x$n)))
    lo2 <- rbind(lo2, lo)
    lo <- lo %>% arrange(n)
    lines(x=lo$n, y=lo$fc)
  }
  print("predict")
  lo_fin <- loess(lo2$fc~lo2$n)
  lo2 <- as.data.frame(cbind(n=log10(x$n), predict(lo_fin, newdata = log10(x$n))))
  lo2 <- lo2 %>% arrange(n)
  lines(x=lo2$n, y=lo2$fc, col="blue", lwd=2)
  
}
