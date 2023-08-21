setwd("C:/Users/iXilef/OneDrive/Master/initial")
load("subset.RData")
annotations <- read.table("./MWU_clone/annotations.BP.go")
colnames(annotations) <- c("id", "go")

# BP_clust is the object containing all DEG stats for the ontology "biological process"
# For now, we only work with fold changes for the contrast called "sep", which is a comparison of migratory and non-migratory butterflies in September.
# Later on, we can look into other contrasts. 

# Extract data for MWU
# exclude non annotated transcripts

X <- unique(BP_clust[,c("seq","FC.sep")]); colnames(X)<-c("id","stat")
not.annot <- as.vector(BP_clust$seq[which(is.na(BP_clust$term))])
X <- X[!X$id %in% not.annot,]
write.table(X,file="X.csv",sep=",",dec=".",quote=F,row.names=F,col.names=T)

# Capital X is the object that contains a line for each combination of GO-term and transcript.
# It's a long format, with a lot of redundancy, but you'll see how it's going to be important to have it like this.

X <- BP_clust[,c("seq","name","term","lev","FC.sep","p.sep")]
colnames(X) <- c("transcript","GOname","GOterm","GOlevel","log2FC","minuslog10pval")
X <- na.omit(X)
# Lower case x is just a line for each transcript. 
x <- X[!duplicated(X$transcript),]
x <- x[,c(1,5,6)]

# Now, first quickly look at the distribution of fold changes:
par(mfrow=c(2,1))
hist(x$log2FC,main=paste("mean: ",round(mean(x$log2FC,na.rm=T),3), ", median: ", round(median(x$log2FC,na.rm=T),3)),breaks=100)
hist(X$log2FC,main=paste("mean: ",round(mean(X$log2FC,na.rm=T),3), ", median: ", round(median(X$log2FC,na.rm=T),3)),breaks=100)

# Here you can start to see that the transcripts that show up in a larger number of GO terms, have a different distribution.
# Another way of looking at this:

plot(density(x$log2FC,na.rm=T),lwd=3, main="Distribution of log2 Fold change values by how many membership in GO terms it has")
# Which transcripts show up more than N times:
par(mfrow=c(1,1))
colcounter=1
for (N in seq(10,300,10))
{
overN<-names(table(X$transcript))[table(X$transcript)>N] # select transcripts over N
lines(density(x$log2FC[x$transcript %in% overN],na.rm=T),col=rgb(0,0+colcounter/30,0.5,1)) # selected transcripts plotted
colcounter=colcounter+1
}
# Each density is a selection of fold changes of fewer transcripts, cutting off at belonging to "N" Go terms. 


# Now the wilcoxon test:
# Just a random sample of 100 GO terms for which we do a simple rank test
# When using x, we calculate the rank based on each transcript featuring one time
# When using X, we calculate the rank based on each transcript featuring the number of times it was assignedhttp://127.0.0.1:31281/graphics/plot_zoom_png?width=1200&height=900 to a GO term

set.seed(3)


pvals.x <- NULL
pvals.X <- NULL
fc      <- NULL
goterms <- NULL
lengths <- NULL

for ( i in 1:length(unique(X$GOterm)))
{

goterm  <- unique(X$GOterm)[i]
targets <- na.omit(X$transcript[X$GOterm==goterm])
print(goterm)
wcX <- wilcox.test(X$log2FC ~ c(X$GOterm==goterm), alternative ="two.sided")
wcx <- wilcox.test(x$log2FC ~ c(x$transcript %in% targets), alternative ="two.sided")

goterms <- c(goterms,as.character(goterm))
lengths <- c(lengths,length(targets))
fc      <- c(fc, median(x$log2FC[x$transcript %in% targets]))
pvals.x <- c(pvals.x,wcx$p.value)
pvals.X <- c(pvals.X,wcX$p.value)
}

wilcox_results <- as.data.frame(cbind(goterms, lengths, fc, pvals.X))
# Here you see two things: when calculating ranks of log2FC in x, you get: 
# * way more significant GO terms
# * lower p-values. So after correcting for multiple testing, the effect will be even more severe!


plot(-log10(pvals.x),-log10(pvals.X));abline(h=-log10(0.01));abline(v=-log10(0.01))


