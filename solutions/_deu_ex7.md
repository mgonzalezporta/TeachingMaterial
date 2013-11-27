```rscript
r=result[result$padjust<0.1,]
r=r[complete.cases(r),]
dim(r)[1]     # number of exons
length(unique(r$geneID))    # # number of genes
```
