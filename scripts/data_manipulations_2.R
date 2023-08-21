library(dplyr)
setwd("C:/Users/iXilef/OneDrive/Master/initial")
load("./dataframes/results_reshuffle.RData")
load("./dataframes/results_normal.RData")
load("./dataframes/results_normal2.RData")
load("./dataframes/misc.RData")

# calculate moving median/mean(fc) for # of members on original runs
seq_info <- results_normal$original$seq_info

orig_fc_N <- data.frame(N = seq(0,max(seq_info$nmember_mwu),1))
orig_fc_N$median_mwu <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_mwu > x)]))
orig_fc_N$median_classic <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_classic_sig > x)]))
orig_fc_N$median_elim <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_elim_sig > x)]))
orig_fc_N$median_weight01 <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_weight_sig > x)]))
orig_fc_N$median_lea <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_lea_sig > x)]))

orig_fc_N$median_mwu_fdr <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_mwu_fdr > x)]))
orig_fc_N$median_classic_fdr <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_classic_sig_fdr > x)]))
orig_fc_N$median_elim_fdr <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_elim_sig_fdr > x)]))
orig_fc_N$median_weight01_fdr <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_weight_sig_fdr > x)]))
orig_fc_N$median_lea_fdr <- sapply(orig_fc_N$N, function(x) median(seq_info$FC.sep[which(seq_info$nmember_lea_sig_fdr > x)]))

orig_fc_N$N_mwu <- sapply(orig_fc_N$N, function(x) length(seq_info$FC.sep[which(seq_info$nmember_mwu > x)]))
orig_fc_N$N_classic <- sapply(orig_fc_N$N, function(x) length(seq_info$FC.sep[which(seq_info$nmember_classic_sig > x)]))
orig_fc_N$N_elim <- sapply(orig_fc_N$N, function(x) length(seq_info$FC.sep[which(seq_info$nmember_elim_sig > x)]))
orig_fc_N$N_weight01 <- sapply(orig_fc_N$N, function(x) length(seq_info$FC.sep[which(seq_info$nmember_weight_sig > x)]))
orig_fc_N$N_lea <- sapply(orig_fc_N$N, function(x) length(seq_info$FC.sep[which(seq_info$nmember_lea_sig > x)]))






# calculate moving median for all random runs
# mwu
get_lines_mwu <- function(fcvals, nvals, iterations){
  int.data <- data.frame(N = seq(0,288,1))
  for(i in 3:(iterations+2)){
    print(i)
    int.data <- cbind(int.data, sapply(0:288, function(x) median(as.numeric(fcvals[,i][which(nvals[,i] > x)]))))
  }
  int.data <- as.data.frame(int.data)
  return(int.data)                     
}

get_lines_ks <- function(fcvals, nvals, iterations){
  int.data <- data.frame(N = seq(0,288,1))
  for(i in 3:(iterations+2)){
    print(i)
    int.data <- cbind(int.data, sapply(0:288, function(x) median(as.numeric(fcvals[,i-1][which(nvals[,i] > x)]), na.rm=T)))
  }
  int.data <- as.data.frame(int.data)
  return(int.data)                     
}

mwu_shuffle <- get_lines_mwu(results_reshuffle$mwu$seq, results_reshuffle$mwu$seq_n,100)
classic_shuffle <- get_lines_ks(results_reshuffle$classic$seq, results_reshuffle$classic$seq_n,100)
elim_shuffle <- get_lines_ks(results_reshuffle$elim$seq, results_reshuffle$elim$seq_n,100)
weight_shuffle <- get_lines_ks(results_reshuffle$weight$seq, results_reshuffle$weight$seq_n,100)
lea_shuffle <- get_lines_ks(results_reshuffle$lea$seq, results_reshuffle$lea$seq_n,100)


# get it from normal run

mwu_norm <- get_lines_mwu(results_normal$mwu$seq, results_normal$mwu$seq_n,100)
classic_norm <- get_lines_ks(results_normal$classic$seq, results_normal$classic$seq_n,100)
elim_norm <- get_lines_ks(results_normal$elim$seq, results_normal$elim$seq_n,100)
weight_norm <- get_lines_ks(results_normal$weight$seq, results_normal$weight$seq_n,100)
lea_norm <- get_lines_ks(results_normal$lea$seq, results_normal$lea$seq_n,100)

# get it from normal2 run

mwu_norm_2 <- get_lines_mwu(results_normal2$mwu$seq, results_normal2$mwu$seq_n,1000)
classic_norm_2 <- get_lines_ks(results_normal2$classic$seq, results_normal2$classic$seq_n,1000)
elim_norm_2 <- get_lines_ks(results_normal2$elim$seq, results_normal2$elim$seq_n,1000)
weight_norm_2 <- get_lines_ks(results_normal2$weight$seq, results_normal2$weight$seq_n,1000)
lea_norm_2 <- get_lines_ks(results_normal2$lea$seq, results_normal2$lea$seq_n,1000)

### at least some pictures

# mwu 

plot(mwu_trans_rand[,3], mwu_trans_rand_n[,3], ylim=c(0, 250), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (MWU FC Shuffle) \n red = real data, grey = random resampling of FC")
for(i in 4:102){
  points(mwu_trans_rand[,i], mwu_trans_rand_n[,i], col="grey")
}

points(seq_info$FC.sep, seq_info$nmember_mwu, col="red")


plot(mwu_trans_norm[,3], mwu_trans_norm_n[,3], ylim=c(0, 250), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (MWU FC N(0, sd(realFC)) ) \n red = real data, grey = random resampling of FC")
for(i in 4:102){
  points(mwu_trans_norm[,i], mwu_trans_norm_n[,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_mwu, col="red")

# classic 

plot(res_rand_n[[1]][,3], res_rand_n[[2]][,3], ylim=c(0, 250), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (classic FC Shuffle) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_rand_n[[1]][,i], res_rand_n[[2]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_classic_sig, col="red")
abline(v=median(seq_info$FC.sep), col="red")



plot(res_norm_n[[1]][,3], res_norm_n[[2]][,3], ylim=c(0, 250), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (Classic FC N(0, sd(realFC)) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_norm_n[[1]][,i], res_norm_n[[2]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_classic_sig, col="red")
abline(v=median(seq_info$FC.sep), col="red")

# elim

plot(res_rand_n[[1]][,3], res_rand_n[[3]][,3], ylim=c(0, 100), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (elim FC Shuffle) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_rand_n[[1]][,i], res_rand_n[[3]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_elim_sig, col="red")
abline(v=median(seq_info$FC.sep), col="red")



plot(res_norm_n[[1]][,3], res_norm_n[[3]][,3], ylim=c(0, 100), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (elim FC N(0, sd(realFC)) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_norm_n[[1]][,i], res_norm_n[[3]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_elim_sig, col="red")

# weight01

plot(res_rand_n[[1]][,3], res_rand_n[[4]][,3], ylim=c(0, 100), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (weight01 FC Shuffle) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_rand_n[[1]][,i], res_rand_n[[4]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_weight_sig, col="red")
abline(v=median(seq_info$FC.sep), col="red")



plot(res_norm_n[[1]][,3], res_norm_n[[4]][,3], ylim=c(0, 100), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (weight01 FC N(0, sd(realFC)) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_norm_n[[1]][,i], res_norm_n[[4]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_weight_sig, col="red")

# lea

plot(res_rand_n[[1]][,3], res_rand_n[[5]][,3], ylim=c(0, 100), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (lea FC Shuffle) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_rand_n[[1]][,i], res_rand_n[[5]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_lea_sig, col="red")
abline(v=median(seq_info$FC.sep), col="red")



plot(res_norm_n[[1]][,3], res_norm_n[[5]][,3], ylim=c(0, 100), xlim=c(-6.5, 6), col="grey", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms (lea FC N(0, sd(realFC)) ) \n red = real data, grey = random resampling of FC")
for(i in 4:101){
  points(res_norm_n[[1]][,i], res_norm_n[[5]][,i], col="grey")
}
points(seq_info$FC.sep, seq_info$nmember_lea_sig, col="red")
abline(v=median(seq_info$FC.sep), col="red")

plot(seq_info$FC.sep, seq_info$nmember_mwu, col="black", ylab="Number of times Transcript was member of significant GO term",
     xlab="log2 Fold change", main="Fold change vs number of membership in significant GO terms
     black=MWU; blue=classic; orange=elim; green=weight; lea=grey
     Vorizontal lines: Mean(FC) | N > 20")
abline(v=mean(seq_info$FC.sep[which(seq_info$nmember_mwu > 20)]), col="black")
points(seq_info$FC.sep, seq_info$nmember_classic_sig, col="blue")
abline(v=mean(seq_info$FC.sep[which(seq_info$nmember_classic_sig > 20)]), col="blue")
points(seq_info$FC.sep, seq_info$nmember_elim_sig, col="orange")
abline(v=mean(seq_info$FC.sep[which(seq_info$nmember_elim_sig > 20)]), col="orange")
points(seq_info$FC.sep, seq_info$nmember_lea_sig, col="green")
abline(v=mean(seq_info$FC.sep[which(seq_info$nmember_weight_sig > 20)]), col="green")
points(seq_info$FC.sep, seq_info$nmember_weight_sig, col="grey")
abline(v=mean(seq_info$FC.sep[which(seq_info$nmember_lea_sig > 20)]), col="grey")