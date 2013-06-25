## Differential gene expression
*Note: This section is based and includes text from the [DESeq package vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf) to which you can refer for more details.*

A basic task in the analysis of expression data is the detection of differentially expressed genes. The *DESeq* package provides a method to test for differential expression by using the negative binomial distribution and a shrinkage estimator for the distribution variance. *DESeq* expects a matrix of count values where each column corresponds to a sample and each line to a feature such as a gene, transcript or exon. *DESeq* uses its own normalisation procedure so the input values must be the raw counts of sequencing reads rather than normalised counts.

In this section we are going to use the pasilla dataset to find genes differentially expressed across the two conditions represented there.

```R
library(DESeq)
library("pasilla")
data("pasillaGenes")

# create an object of class CountDataSet, 
# which is the data container used by DESeq
cds = newCountDataSet(counts(pasillaGenes), phenoData(pasillaGenes)$condition)
cds

# the CountDataSet class in a container for high-throughput 
# assays and experimental metadata
# the counts can be accessed with the counts function
head( counts(cds) )
# and the phenotypical data with the pData function
pData(cds)
```

The first processing step is to normalise the read counts by library size. For this we start by determining the relative library sizes, or *size factors* of each library. For example, if the counts of non-differentially expressed genes in one sample are, on average, twice as high as in another, the size factor for the first sample should be twice as large as the one for the other sample. These size factors can be obtained with the function *estimateSizeFactors*:

```R
cds <- estimateSizeFactors( cds )
sizeFactors( cds )

# the normalised data is obtained by dividing each column of 
# the count table by the size factor for this column
# this calculation is performed by calling the function counts 
# with normalized=TRUE
head( counts( cds, normalized=TRUE ) )
```

#### Variance inference
DESeq relies on an estimation of the typical relationship between the data dispersion and their mean. The dispersion can be understood as the square of the coefficient of biological variation. For example, if the expression of a gene  typically differs from replicate to replicate sample by 20%, the dispersion is 0.20^2 = 0.04. Note that the variance seen between counts is the sum of two components: the sample-to-sample variation just mentioned, and the uncertainty in measuring a concentration by counting reads. The latter, known as Poisson noise, is the dominating noise source for lowly expressed genes. The sum of both Poisson noise and dispersion is considered in the differential expression inference. To estimate the dispersions we can use the function *estimateDispersions*, which performs three steps. First, it estimates, for each gene, a dispersion value, then, it fits a curve through the estimates. Finally, it assigns to each gene a dispersion value, using either the estimated or the fitted value:

```R
cds <- estimateDispersions( cds )
str( fitInfo( cds ) )
head( fData( cds ) )
```

There are several ways in which the empirical dispersion can be computed (by setting the argument *method* of *estimateDispersions*). By default, *estimateDispersions* uses the samples from all conditions to estimate a single pooled empirical dispersion value and assigns it to all samples. The list returned by the *fitInfo* function contains the empirical dispersion values for each gene, the curve fitted through the dispersions, and the fitted values. This can be visualised by plotting the per-gene estimates against the normalised mean expressions per gene, and overlaying the fitted curve (in red):

```R
plot( rowMeans( counts( cds, normalized=TRUE ) ), 
  fitInfo(cds)$perGeneDispEsts, log="xy" )
xg <- 10^seq( -.5, 5, length.out=300 ) 
lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
```

You now have to decide whether to use the empirical dispersions derived only from the data or the fitted values in the differential expression test. By default, *estimateDispersions* takes the maximum of the two values. This is the conservative choice, recommended once you have at least three or four replicates. If you had used less replicates you would have estimated the dispersion from less samples and should expect the estimates to scatter with quite some sampling variance around their true values. In this case you should consider using only the fitted values by using  `estimateDispersions(cds, sharingMode = "fit-only")`. Find out about other *sharingMode* values by looking in the manual pages:

```R
?estimateDispersions
```
              
Finally, to contrast two conditions, in our case to find whether there is differential expression between conditions untreated and treated, we use the function *nbinomTest*:

```R
res <- nbinomTest( cds, "treated", "untreated" )
head(res)
```

The padj column in table `res` contains the p-values, adjusted for multiple testing with the Benjamini-Hochberg procedure, which controls false discovery rate (FDR). We can now plot the normalised mean versus log_2 fold change and colour the genes that are significantly differentially expressed at a 10% FDR in red:
 
```R
plot( res$baseMean, res$log2FoldChange, log="x", 
  pch=19, cex=.3, col = ifelse( res$padj < .1, "red", "black" ), ylim=c(-5,5) )
abline(h=0)

# filter significant genes with a threshold for 
# false dicovery rate (FDR) of 10%
resSig <- res[ res$padj < 0.1, ]

# list the most significant
head( resSig[ order(resSig$pval), ] )
```