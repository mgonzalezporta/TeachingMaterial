```rconsole
plot( res$baseMean, res$log2FoldChange, log="x", 
    pch=19, cex=.3, col = ifelse( res$padj < .1, "red", "black" ), 
    ylim=c(-5,5) )
abline(h=0)
```