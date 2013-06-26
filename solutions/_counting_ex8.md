```rconsole
split_counts=split(exon_counts_chr4, names(exon_counts_chr4) ) 
head(split_counts)
gene_counts_chr4=sapply( split_counts, function(x) sum(x) ) 
head(gene_counts_chr4)
```