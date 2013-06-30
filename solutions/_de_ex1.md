```rconsole
de_significant=de[de$padj < 0.1, ]
head(de_significant[order(de_significant$pval), ])
```
