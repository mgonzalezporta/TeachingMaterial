```rconsole
# retrieve exon lengths 
exon_lengths=width(annot) 
names(exon_lengths)=elementMetadata(annot)$gene
head(exon_lengths)

# calculate gene lengts 
exon_lengths_by_gene=split( exon_lengths, names(exon_lengths) )
head(exon_lengths_by_gene, n=1)
gene_lengths=sapply( exon_lengths_by_gene, sum ) 
head(gene_lengths)

# normalise counts by library size 
lib_size=colSums(counts) 
ncounts=t(t(counts)/lib_size)

# normalise counts by gene length 
common_genes=intersect(row.names(ncounts), names(gene_lengths)) 
subset_ncounts=ncounts[row.names(ncounts) %in% common_genes,] 
gene_lengths=gene_lengths[names(gene_lengths) %in% common_genes] 
ncounts=subset_ncounts/gene_lengths

# obtain RPKMs 
rpkm=ncounts*1e9 
head(rpkm)
```