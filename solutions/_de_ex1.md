```rconsole
resSig=res[ res$padj < 0.1, ]
head( resSig[ order(resSig$pval), ] )
```