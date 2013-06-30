In a typical differential expression analysis, the majority of the genes are assumed to be non-differentially expressed. For this reason, it is useful to examine boxplots of counts across samples in each dataset, both before and after normalization. An effective normalization should result in a stabilization of read count distributions across replicates.

```rconsole
# create a list to simplify the plotting step
l=list(
        split(counts, colnames(counts)),
        split(rpkm, colnames(rpkm)),
        split(deseq_ncounts, colnames(deseq_ncounts))
)
names(l)=c("raw counts", "normalised counts - RPKMs", "normalised counts - DESeq")

# visualise the data as a set of boxplots
pdf(file="./normalisingg_ex2.pdf", width=10)
par(mfrow=c(1,3), las=2, mar=c(7,5,3,3), cex=1)
boxplot(l[[1]], outline=F, ylab="raw counts", ylim=c(0, 1100))
boxplot(l[[2]], outline=F, ylab="RPKMs", ylim=c(0, 60))
boxplot(l[[3]], outline=F, ylab="normalised counts - DESeq", ylim=c(0, 1100))
dev.off()
```
