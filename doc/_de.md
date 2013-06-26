## Differential gene expression
A basic task in the analysis of expression data is the detection of differentially expressed genes. The *DESeq* package provides a method to test for differential expression by using the negative binomial distribution. It expects a matrix of count values where each column corresponds to a sample and each line to a feature (e.g. a gene). Typically, a *DESeq* analysis is performed in three steps: count normalisation, dispersion estimation and differential expression test.

### Count normalisation
Since we have already generated a matrix with the normalised counts in the previous section, we will use it directly as input for the next step.

### Dispersion estimation
An important step in differential expression analysis is to figure out how much variability we can expect in the expression measurements within the same condition. Unless this is known, we cannot make inferences about whether the change in expression observed for a given gene is big enough to be considered significant, or whether it corresponds to the variability that we would expect by chance. This is why it is so important to have replicates: they show how much variation occurs without a difference in the condition.


In *DESeq*, in order to estimate the dispersion for each gene, we can use the function *estimateDispersions*:

```rconsole
cds=estimateDispersions(cds)
str(fitInfo(cds))
head(fData(cds))
```

The list returned by the *fitInfo* function contains the empirical dispersion values for each gene, the curve fitted through the dispersions, and the fitted values. This can be visualised by plotting the per-gene estimates against the normalised mean expressions per gene, and overlaying the fitted curve:

```rconsole
plot(rowMeans(counts(cds, normalized=TRUE)), fitInfo(cds)$perGeneDispEsts, log="xy")
xg=10^seq(-.5, 5, length.out=300) 
lines(xg, fitInfo(cds)$dispFun(xg), col="red")
```

### Differential expression test
              
Finally, we will use the function *nbinomTest* to contrast the two studied conditions:

```rconsole
de=nbinomTest(cds, "treated", "untreated")
head(de)
```

The *padj* column in the table `de` contains the p-values adjusted for multiple testing with the Benjamini-Hochberg procedure (i.e. FDR). This is the information that we will use to decide whether the expression of a given gene differs significantly across conditions (e.g. we can arbitrarily decide that genes with an FDR<0.10 are differentially expressed).

**Exercise:** How would you select those genes that pass a given FDR threshold (e.g. FDR<0.10)? Which are the most significant?
[Solution](../solutions/_de_ex1.md)

**Exercise:** Let us generate an MA plot to evaluate the results of the differential expression analysis. Follow this steps:

* Obtain the mean expression values from the `de` object.

* Obtain the log2 fold change values from the `de` object.

* Create a scatterplot with the mean in the *x* axis and the log2 fold change in the *y* axis.

* How would we highlight those genes that are differentially expressed with a different color?
  *Hint:* check the *col* argument for the *plot* function and the function *ifelse*.

[Solution](../solutions/_de_ex2.md)