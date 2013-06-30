```rconsole
pdf(file="./de_ex2.pdf")
plot( de$baseMean, de$log2FoldChange, log="x", 
    pch=19, cex=.3, col = ifelse( de$padj < .1, "red", "black" ), 
    ylim=c(-5,5) )
abline(h=0)
dev.off()
```
