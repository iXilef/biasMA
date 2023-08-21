# I first generate a simple mock data set with some bias:
# 1000 features (genes, transcripts, whatever), 100 terms, each feature member of some terms. 
# The way I set this up is as follows:
# I draw from a normal distribution with mean = 0 a set of fold changes.
# I make a membership matrix with 1000 rows (one for each feature) and 100 columns (one for each GO term).
# In this matrix, each row corresponds to the memberships (logical) of the corresponding feature in the "values" vector
# Now in this matrix, I increase the number of memberships, with only about 20 memberships for the first features, and up to 60 memberships for the last ones
# If now we sort the values vector (and add some noise), we have recreated our bias. 

for(i in seq(0,0.3,.1)){
# I define log2FC values by adding a sorted random distribution to an unsorted one. By playing around with the standard deviation of either of them, I can make the bias stronger or weaker.
values          <- sort(rnorm(1000,0,i)) + rnorm(1000,0,1)
values.control  <- sample(values)

# Quick and dirty way of generating a matrix with increasing numbers of memberships. 
start_prob <- 0.2
end_prob   <- 0.6
probabilities <- seq(start_prob, end_prob, length.out = 100000)
memberships <- rep(F,100000)
for (i in 1:100000)
{
 memberships[i] <- sample(c(TRUE,FALSE),1,prob=c(probabilities[i],1-probabilities[i]))
}
	
memberships <- matrix (memberships,ncol=100,byrow=T)

value_multiplied <- rep(values, apply(memberships, 1, function(x) length(which(x==T))))

# Now let's do some statistical tests for each term:
pvals.t   <- rep(0,100)
pvals.mwu <- rep(0,100)
pvals.ks  <- rep(0,100) 

pvals.control.t   <- rep(0,100)
pvals.control.mwu <- rep(0,100)
pvals.control.ks  <- rep(0,100) 

pvals.incl.t   <- rep(0,100)
pvals.incl.mwu <- rep(0,100)
pvals.incl.ks  <- rep(0,100) 

pvals.all.t   <- rep(0,100)
pvals.all.mwu <- rep(0,100)
pvals.all.ks  <- rep(0,100) 

for (i in 1:100)
{
	pvals.t[i]   <- -log10(t.test(values~memberships[,i])$p.value)	
	pvals.mwu[i] <- -log10(wilcox.test(values~memberships[,i])$p.value)
	pvals.ks[i]  <- -log10(ks.test(values~memberships[,i])$p.value)

	pvals.control.t[i]   <- -log10(t.test(values.control~memberships[,i])$p.value)	
	pvals.control.mwu[i] <- -log10(wilcox.test(values.control~memberships[,i])$p.value)
	pvals.control.ks[i]  <- -log10(ks.test(values.control~memberships[,i])$p.value)	
	
	pvals.incl.t[i]   <- -log10(t.test(values[memberships[,i]], values)$p.value)	
	pvals.incl.mwu[i] <- -log10(wilcox.test(values[memberships[,i]], values)$p.value)
	pvals.incl.ks[i]  <- -log10(ks.test(values[memberships[,i]], values)$p.value)
	
	pvals.all.t[i]   <- -log10(t.test(values[memberships[,i]], value_multiplied)$p.value)
	pvals.all.mwu[i] <- -log10(wilcox.test(values[memberships[,i]], value_multiplied)$p.value)
	pvals.all.ks[i]  <- -log10(ks.test(values[memberships[,i]], value_multiplied)$p.value)
	

}
MAX <- max(c(pvals.t,pvals.mwu,pvals.ks)) + 1

par(mfrow=c(2,4))
# hist(pvals.t,xlim=c(0,MAX),breaks=20);abline(v=-log10(0.05),col="red",main="t test, with bias")
# hist(pvals.incl.t,xlim=c(0,MAX),breaks=20);abline(v=-log10(0.05),col="red",main="t test, with bias")
# hist(pvals.control.t,xlim=c(0,MAX),breaks=20);abline(v=-log10(0.05),col="red",main="t test, no bias")

hist(pvals.ks,xlim=c(0,MAX),breaks=20,main="ks test, with bias");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.ks < -log10(0.05)))))
hist(pvals.incl.ks,xlim=c(0,MAX),breaks=20,main="ks test, with bias, whole universe");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.incl.ks < -log10(0.05)))))
hist(pvals.all.ks, xlim=c(0,MAX), breaks=20, main="ks test, with bias, multiplied bg");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.all.ks < -log10(0.05)))))
hist(pvals.control.ks,xlim=c(0,MAX),breaks=20,main="ks test, no bias");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.control.ks < -log10(0.05)))))


hist(pvals.mwu,xlim=c(0,MAX),breaks=20,main="mwu test, with bias");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.mwu < -log10(0.05)))))
hist(pvals.incl.mwu,xlim=c(0,MAX),breaks=20,main="mwu test, with bias, whole universe");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.incl.mwu < -log10(0.05)))))
hist(pvals.all.mwu, xlim=c(0,MAX), breaks=20, main="mwu test, with bias, multiplied bg");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.all.mwu < -log10(0.05)))))
hist(pvals.control.mwu,xlim=c(0,MAX),breaks=20,main="mwu test, no bias");abline(v=-log10(0.05),col="red")
text(2, 8, paste0("N (< 0.05)= ", length(which(pvals.control.mwu < -log10(0.05)))))

length(which(pvals.t < -log10(0.05)))

length(which(pvals.mwu < -log10(0.05)))

length(which(pvals.ks < -log10(0.05)))

length(which(pvals.incl.t < -log10(0.05)))

length(which(pvals.incl.mwu < -log10(0.05)))

length(which(pvals.incl.ks < -log10(0.05)))


# Signal obviously correlated between tests
# par(mfrow=c(1,1))
# plot(pvals.mwu,pvals.ks,xlim=c(0,MAX),ylim=c(0,MAX))
}

load("this.RData")
resis2 <- resis[,c(1,5,9,13,
                   3,7,11,15,
                   2,6,10,14,
                   4,8,12,16)]
test <- as.list(as.data.frame(resis2))
test2 <- sapply(test, function(x) x[which(x != 0)])
test3 <- sapply(test2, log10)

rep(seq(0,.3,0.1), each=4)
rep(c("MWU", "KS", "MWU dupl", "KS dupl"), times=4)
names(test3) <- paste(rep(seq(0.0,0.3, 0.1), times=4))
text1 <- paste0("N(0)=", sapply(test2, function(x) 100-length(x)))

par(mfrow=c(1,1))

boxplot(x=test3, las=2, yaxt="n")
axis(2, at=log10(c(1,2,3,5,7,10,15,20,30,50,75,100,150,250,500)), labels=c(1,2,3,5,7,10,15,20,30,50,75,100,150,250,500), las=1)
lines(x=c(1:16),y=log10(sapply(test2, function(x) 100-length(x))), col="blue")
text(x=c(1:16)-0.1, y=log10(.9), text1)
abline(v=c(4.5, 8.5, 12.5), h=log10(75))
text(x=c(2,6,10,14), y=log10(800), labels = c("MWU", "MWU dupl", "KS", "KS dupl"))
text(x=0.8, y= log10(65), labels=c("N(0)"), col="blue")

