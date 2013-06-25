## Normalising counts
### With RPKMs
While in the previous sections the data was derived from a single sample, in this exercise we will work with the precomputed counts for all the samples in [our experiment](http://bioconductor.org/packages/2.11/data/experiment/html/pasilla.html):

```R
library("pasilla")
data("pasillaGenes")
counts=counts(pasillaGenes)
head(counts)
```

A common way to normalise reads is to convert them to RPKM. This implies normalising the read counts depending on the feature size (exon, transcript, gene model...) and on the total number of reads sequenced for that library:
![rpkms](../img/rpkms/png)

**Exercise:** This exercise may be challenging but if you find you need help the solution is provided at the end of this document. Obtain RPKMs for the table `counts` following these hints:

* find out what the length of each gene is by summing its exon lengths: 
    * try to use the accessors of the `annot` object to get the individual exon lengths and store it in an object called `exon_lengths`        
    * get the list of genes from `annot` with the *elementMetadata* accessor and store it in an object called `gene_names`; notice that some genes are repeated since they correspond to the different exons in `annot` and that they correspond to the exon lengths you just determined 
    * use the *split* function on the exon lengths and the gene names to create a list called `length_by_gene` in which each element 
    corresponds to a gene and contains a vector of the lengths of each of its exons 
    * use the function *sapply* on the list you have just created to sum the exon length vectors for each gene with 
    `sapply( length_by_gene, sum )` and store it in `gene_lengths`
* normalise counts by library size:
    * use *colSums* to find out the total number of reads per library and store it in an object called `lib_size`
    * divide the table `counts` by `lib_size` and save the result to an object called `ncounts`
* normalise counts by gene length, keeping in mind that we were working with a subset of the annotation, therefore we just want to keep genes in chromosome 4
* obtain RPKM values by multiplying by a factor of 10^9

Such a count normalisation is suited for visualisation, but sub-optimal for further analyses. A better way of normalising the data is to use either the *edgeR* or *DESeq* packages.

### With DESeq
*Note: This section is based and includes text from the [DESeq package vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf) to which you can refer for more details.*

Despite of the easy interpretation, RPKM normalisation is not the most adequate for certain types of downstream analysis (eg differential gene expression), given that it is susceptible to library normalisation biases. There are many other normalisation methods that should be considered with that goal in mind (see \cite{Dillies:2012gz} for a comparison). In this section we are going to explore the one offered within the DESeq package:

```R
\begin{verbatim}
library(DESeq)

# create an object of class CountDataSet, 
# which is the data container used by DESeq
# remember we have already loaded the count data in the previous section
cds = newCountDataSet(counts(pasillaGenes), phenoData(pasillaGenes)$condition)
cds

# the CountDataSet class in a container for high-throughput 
# assays and experimental metadata
# the counts can be accessed with the counts function
head( counts(cds) )
# and the phenotypical data with the pData function
pData(cds)
```

In order to normalise the raw counts we start by determining the relative library sizes, or *size factors* of each library. For example, if the counts of the expressed genes in one sample are, on average, twice as high as in another, the size factor for the first sample should be twice as large as the one for the other sample. These size factors can be obtained with the function *estimateSizeFactors*:

```R
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
```

Once we have this information, the normalised data is obtained by dividing each column of the count table by the corresponding size factor. We can perform this calculation by calling the function counts with a specific argument as follows:

```R
norm_counts=counts( cds, normalized=TRUE )
```

**Exercise:** We have now accumulated three different versions of the same dataset: the raw counts (`counts`), the RPKM quantifications (`rpkm`) and the DESeq normalised counts (`norm_counts`). How would you visualise the performance of each in getting rid of the variation that does not associate to the experimental conditions that are being evaluated?