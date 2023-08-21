library(fasano.franceschini.test)

# median(fc) over N original runs

plot(orig_fc_N$N, orig_fc_N$median_mwu, col="black", type="l")
text(orig_fc_N$N, orig_fc_N$median_mwu, labels=orig_fc_N$N_mwu)
lines(orig_fc_N$N, orig_fc_N$median_classic, col="blue")
text(orig_fc_N$N, orig_fc_N$median_classic, labels=orig_fc_N$N_classic)
lines(orig_fc_N$N, orig_fc_N$median_elim, col="orange")
text(orig_fc_N$N, orig_fc_N$median_elim, labels=orig_fc_N$N_elim)
lines(orig_fc_N$N, orig_fc_N$median_weight01, col="green")
text(orig_fc_N$N, orig_fc_N$median_weight01, labels=orig_fc_N$N_weight01)
lines(orig_fc_N$N, orig_fc_N$median_lea, col="gray")
text(orig_fc_N$N, orig_fc_N$median_lea, labels=orig_fc_N$N_lea)


# lines for shuffle runs

plot(mwu_shuffle$N, mwu_shuffle[,2], type="l", ylim=c(-1.6,2.4))
for(i in 2:101){
  lines(mwu_shuffle$N, mwu_shuffle[,i])
}
lines(1:289, apply(mwu_shuffle[,c(2:101)],1, mean), col="red")

plot(classic_shuffle$N, classic_shuffle[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,120))
for(i in 3:101){
  lines(classic_shuffle$N, classic_shuffle[,i])
}
lines(1:289, apply(classic_shuffle[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

plot(elim_shuffle$N, elim_shuffle[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:101){
  lines(elim_shuffle$N, elim_shuffle[,i])
}
lines(1:289, apply(elim_shuffle[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

plot(weight_shuffle$N, weight_shuffle[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:101){
  lines(weight_shuffle$N, weight_shuffle[,i])
}
lines(1:289, apply(weight_shuffle[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

# lines for normal runs

plot(mwu_norm$N, mwu_norm[,2], type="l", ylim=c(-1.6,2.4))
for(i in 2:101){
  lines(mwu_norm$N, mwu_norm[,i])
}
lines(1:289, apply(mwu_norm[,c(2:101)],1, mean), col="red")

plot(classic_norm$N, classic_norm[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,120))
for(i in 3:101){
  lines(classic_norm$N, classic_norm[,i])
}
lines(1:289, apply(classic_norm[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

plot(elim_norm$N, elim_norm[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:101){
  lines(elim_norm$N, elim_norm[,i])
}
lines(1:289, apply(elim_norm[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

plot(weight_norm$N, weight_norm[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:101){
  lines(weight_norm$N, weight_norm[,i])
}
lines(1:289, apply(weight_norm[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

plot(lea_norm$N, lea_norm[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:101){
  lines(lea_norm$N, lea_norm[,i])
}
lines(1:289, apply(lea_norm[,c(2:101)],1, function(x) mean(x, na.rm=T)), col="red")

# lines for normal2 runs

plot(mwu_norm_2$N, mwu_norm_2[,2], type="l", ylim=c(-1.6,2.4))
for(i in 2:1001){
  lines(mwu_norm_2$N, mwu_norm_2[,i])
}
lines(1:289, apply(mwu_norm_2[,c(2:1001)],1, mean), col="red")

plot(classic_norm_2$N, classic_norm_2[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,120))
for(i in 3:1001){
  lines(classic_norm_2$N, classic_norm_2[,i])
}
lines(1:289, apply(classic_norm_2[,c(2:1001)],1, function(x) mean(x, na.rm=T)), col="red")

plot(elim_norm_2$N, elim_norm_2[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:1001){
  lines(elim_norm_2$N, elim_norm_2[,i])
}
lines(1:289, apply(elim_norm_2[,c(2:1001)],1, function(x) mean(x, na.rm=T)), col="red")

plot(weight_norm_2$N, weight_norm_2[,2], type="l", ylim=c(-2.7,3.4), xlim=c(0,50))
for(i in 3:1001){
  lines(weight_norm_2$N, weight_norm_2[,i])
}
lines(1:289, apply(weight_norm_2[,c(2:1001)],1, function(x) mean(x, na.rm=T)), col="red")

plot(lea_norm_2$N, lea_norm_2[,2], type="l", ylim=c(-1,1), xlim=c(0,50))
for(i in 3:1001){
  lines(lea_norm_2$N, lea_norm_2[,i], col=ifelse(results_normal2$metadata$ks$mean[i-1] > 0.5 & results_normal2$metadata$ks$percentage[i-1] > 15, "red", "grey"))
}
lines(1:289, apply(lea_norm_2[,c(2:1001)],1, function(x) mean(x, na.rm=T)), col="red")


# small peacock test try
# DONT RUN TAKES HOURS
# library(fasano.franceschini.test)
# 
# mwu_file <- results_reshuffle$original$seq_info[,c("nmember_mwu", "FC.sep")]
# reshuffle_file <- cbind(as.numeric(results_reshuffle$mwu$seq[,3]), as.numeric(results_reshuffle$mwu$seq_n[,3]))
# for(i in 2:100){
#   print(i)
#   reshuffle_file <- rbind(reshuffle_file, cbind(as.numeric(results_reshuffle$mwu$seq[,i+2]), as.numeric(results_reshuffle$mwu$seq_n[,i+2])))
# }
# 
# fasano.franceschini.test(mwu_file, reshuffle_file, threads=12)

# linear models 
mwu_lin <- lm(log(seq_info$nmember_mwu+1)~seq_info$FC.sep)
summary(mwu_lin)
plot(mwu_lin)

classic_lin <- lm(log(seq_info$nmember_classic_sig+1)~seq_info$FC.sep)
summary(classic_lin)
plot(classic_lin)

elim_lin <- lm(log(seq_info$nmember_elim_sig+1)~seq_info$FC.sep)
summary(elim_lin)
plot(elim_lin)

weight_lin <- lm(log(seq_info$nmember_weight_sig+1)~seq_info$FC.sep)
summary(weight_lin)
plot(weight_lin)


# pvals <- c(fasano.franceschini.test(seq_info[1, c(11,12)], cbind(as.numeric(t(results_reshuffle$mwu$seq[1,])), as.numeric(t(results_reshuffle$mwu$seq_n[1,])))[-c(1,2),], nPermute = 100, seed=12345)$p.value)
pvals <- c()
for(i in seq_info$seq){
  if(!(i %in% results_reshuffle$mwu$seq$transcript)){next}
  # print(i)
  pvals <- c(pvals, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                             cbind(as.numeric(t(results_reshuffle$mwu$seq[which(results_reshuffle$mwu$seq$transcript==i),])), 
                                                   as.integer(t(results_reshuffle$mwu$seq_n[which(results_reshuffle$mwu$seq$transcript==i),])))[-c(1,2),], 
                                             nPermute = 100, verbose=FALSE)$p.value)
}

pvals_clas <- c()
for(i in seq_info$seq){
  if(!(i %in% results_reshuffle$classic$seq$transcript)){next}
  # print(i)
  pvals_clas <- c(pvals_clas, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                             cbind(as.numeric(t(results_reshuffle$classic$seq[which(results_reshuffle$classic$seq$transcript==i),])), 
                                                   as.integer(t(results_reshuffle$classic$seq_n[which(results_reshuffle$classic$seq$transcript==i),])))[-c(1,2),], 
                                             nPermute = 100, verbose=FALSE)$p.value)
}

pvals_elim <- c()
for(i in seq_info$seq){
  if(!(i %in% results_reshuffle$elim$seq$transcript)){next}
  # print(i)
  pvals_elim <- c(pvals_elim, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                       cbind(as.numeric(t(results_reshuffle$elim$seq[which(results_reshuffle$elim$seq$transcript==i),])), 
                                                             as.integer(t(results_reshuffle$elim$seq_n[which(results_reshuffle$elim$seq$transcript==i),])))[-c(1,2),], 
                                                       nPermute = 100, verbose=FALSE)$p.value)
}

pvals_weight <- c()
for(i in seq_info$seq){
  if(!(i %in% results_reshuffle$weight$seq$transcript)){next}
  # print(i)
  pvals_weight <- c(pvals_weight, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                       cbind(as.numeric(t(results_reshuffle$weight$seq[which(results_reshuffle$weight$seq$transcript==i),])), 
                                                             as.integer(t(results_reshuffle$weight$seq_n[which(results_reshuffle$weight$seq$transcript==i),])))[-c(1,2),], 
                                                       nPermute = 100, verbose=FALSE)$p.value)
}

pvals_lea <- c()
for(i in seq_info$seq){
  if(!(i %in% results_reshuffle$lea$seq$transcript)){next}
  # print(i)
  pvals_lea <- c(pvals_lea, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                       cbind(as.numeric(t(results_reshuffle$lea$seq[which(results_reshuffle$lea$seq$transcript==i),])), 
                                                             as.integer(t(results_reshuffle$lea$seq_n[which(results_reshuffle$lea$seq$transcript==i),])))[-c(1,2),], 
                                                       nPermute = 100, verbose=FALSE)$p.value)
}


pvals_norm <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal$mwu$seq$transcript)){next}
  # print(i)
  pvals_norm <- c(pvals_norm, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                             cbind(as.numeric(t(results_normal$mwu$seq[which(results_normal$mwu$seq$transcript==i),])), 
                                                   as.integer(t(results_normal$mwu$seq_n[which(results_normal$mwu$seq$transcript==i),])))[-c(1,2),], 
                                             nPermute = 100, verbose=FALSE)$p.value)
}

pvals_clas_norm <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal$classic$seq$transcript)){next}
  # print(i)
  pvals_clas_norm <- c(pvals_clas_norm, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                       cbind(as.numeric(t(results_normal$classic$seq[which(results_normal$classic$seq$transcript==i),])), 
                                                             as.integer(t(results_normal$classic$seq_n[which(results_normal$classic$seq$transcript==i),])))[-c(1,2),], 
                                                       nPermute = 100, verbose=FALSE)$p.value)
}

pvals_elim_norm <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal$elim$seq$transcript)){next}
  # print(i)
  pvals_elim_norm <- c(pvals_elim_norm, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                       cbind(as.numeric(t(results_normal$elim$seq[which(results_normal$elim$seq$transcript==i),])), 
                                                             as.integer(t(results_normal$elim$seq_n[which(results_normal$elim$seq$transcript==i),])))[-c(1,2),], 
                                                       nPermute = 100, verbose=FALSE)$p.value)
}

pvals_weight_norm <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal$weight$seq$transcript)){next}
  # print(i)
  pvals_weight_norm <- c(pvals_weight_norm, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                           cbind(as.numeric(t(results_normal$weight$seq[which(results_normal$weight$seq$transcript==i),])), 
                                                                 as.integer(t(results_normal$weight$seq_n[which(results_normal$weight$seq$transcript==i),])))[-c(1,2),], 
                                                           nPermute = 100, verbose=FALSE)$p.value)
}

pvals_lea_norm <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal$lea$seq$transcript)){next}
  # print(i)
  pvals_lea_norm <- c(pvals_lea_norm, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                     cbind(as.numeric(t(results_normal$lea$seq[which(results_normal$lea$seq$transcript==i),])), 
                                                           as.integer(t(results_normal$lea$seq_n[which(results_normal$lea$seq$transcript==i),])))[-c(1,2),], 
                                                     nPermute = 100, verbose=FALSE)$p.value)
}

pvals_norm_2 <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal2$mwu$seq$transcript)){next}
  # print(i)
  pvals_norm_2 <- c(pvals_norm_2, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                       cbind(as.numeric(t(results_normal2$mwu$seq[which(results_normal2$mwu$seq$transcript==i),])), 
                                                             as.integer(t(results_normal2$mwu$seq_n[which(results_normal2$mwu$seq$transcript==i),])))[-c(1,2),], 
                                                       nPermute = 100, verbose=FALSE)$p.value)
}

pvals_clas_norm_2 <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal2$classic$seq$transcript)){next}
  # print(i)
  pvals_clas_norm_2 <- c(pvals_clas_norm_2, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                                 cbind(as.numeric(t(results_normal2$classic$seq[which(results_normal2$classic$seq$transcript==i),])), 
                                                                       as.integer(t(results_normal2$classic$seq_n[which(results_normal2$classic$seq$transcript==i),])))[-c(1,2),], 
                                                                 nPermute = 100, verbose=FALSE)$p.value)
}

pvals_elim_norm_2 <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal2$elim$seq$transcript)){next}
  # print(i)
  pvals_elim_norm_2 <- c(pvals_elim_norm_2, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                                 cbind(as.numeric(t(results_normal2$elim$seq[which(results_normal2$elim$seq$transcript==i),])), 
                                                                       as.integer(t(results_normal2$elim$seq_n[which(results_normal2$elim$seq$transcript==i),])))[-c(1,2),], 
                                                                 nPermute = 100, verbose=FALSE)$p.value)
}

pvals_weight_norm_2 <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal2$weight$seq$transcript)){next}
  # print(i)
  pvals_weight_norm_2 <- c(pvals_weight_norm_2, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                                     cbind(as.numeric(t(results_normal2$weight$seq[which(results_normal2$weight$seq$transcript==i),])), 
                                                                           as.integer(t(results_normal2$weight$seq_n[which(results_normal2$weight$seq$transcript==i),])))[-c(1,2),], 
                                                                     nPermute = 100, verbose=FALSE)$p.value)
}

pvals_lea_norm_2 <- c()
for(i in seq_info$seq){
  if(!(i %in% results_normal2$lea$seq$transcript)){next}
  # print(i)
  pvals_lea_norm_2 <- c(pvals_lea_norm_2, fasano.franceschini.test(seq_info[which(seq_info$seq==i), c(11,12)],
                                                               cbind(as.numeric(t(results_normal2$lea$seq[which(results_normal2$lea$seq$transcript==i),])), 
                                                                     as.integer(t(results_normal2$lea$seq_n[which(results_normal2$lea$seq$transcript==i),])))[-c(1,2),], 
                                                               nPermute = 100, verbose=FALSE)$p.value)
}