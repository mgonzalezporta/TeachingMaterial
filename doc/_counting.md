## Counting reads overlapping annotated genes
Following the read mapping step, we can proceed working with BAM files with standalone tools or load them directly in R. These two worfkflows are not exclusive and we will cover both of them for illustrative purposes.

### With htseq-count
[htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) is a simple but yet powerful tool to overlap a BAM file with the genome annotation and thus obtain the number of reads that overlap with our features of interest. As usual, we can obtain information on the tool by typing `htseq-count -h` and by referring to its website (see above). 

**Exercise:** 0ne of the input files required by htseq-count is a GTF file. For this practical, you will find this file under the directory `reference`. Which information does it contain?

**Exercise:** As we have already mentioned, the other required input file is a BAM file. Can you spot any specific requirement regarding this file?
*Hint:* check the website of the tool.

In addition to the input file requirements, special care must be taken in dealing with reads that overlap more than one feature (e.g. overlapping genes), and thus might be counted several times in different features. To deal with this, htseq-count offers three different counting modes: union, intersection-strict and intersection-nonempty.

**Exercise:** What are the differences between these three counting modes?
*Hint:* check the website of the tool.

Now that we have a good understanding of the input files and options, we can proceed to execute htseq-count:

```bash
# the following command takes a while to execute
# the output file is already provided in the data/mapped directory
samtools view untreated3_paired.bam | htseq-count \
    --mode=intersection-nonempty \
    --stranded=no \
    --type=exon \
    --idattr=gene_id - \
    ../../reference/Drosophila_melanogaster.BDGP5.25.62.gtf > htseq_count.out
```

**Exercise:** In addition to the counts that overlap known genes, the output file also contains some extra information on reads that could not be assigned to any of those; can you find it?

### With R
Computing gene counts in R is very similar to what we have done so far with htseq-count. However, it requires some extra steps, since we first need to load the necessary files (i.e. BAM files and annotation).

*Note: All the commands provided in this section have to be executed in R. Make sure to specify the working directory properly before starting (e.g. `setwd("./data")`).*

#### Importing BAM files
There are three main functions to load BAM files into R:

* *scanBam*: this function is part of the *Rsamtools* package and is the low level function used by the other two. It potentially reads *all* fields (including CIGAR strings and user defined tags) of a BAM file into a list structure, but allows you to select specific fields and records to import.
* *readAligned*: a higher-level function defined in the *ShortRead* package which imports some of the data (query names, sequences, quality, strand, reference name, position, mapping quality and flag) into an *AlignedRead* object (which an extension of a *ShortReadQ( object, and thus has the same functions). *ShortRead( was the first package developed to read in NGS data and is able to read almost sequencer every manufacturer proprietary formats (with the notable exception of ABI color-space), so you could for example also use it to read an Illumina export file produced by a GenomeAnalyzer GAIIx.
* *readGappedAlignments* and *readGappedAlignmentPairs*: two functions from the *GenomicRanges* package that create an object intended for operations such as searching for overlaps or coverage. Each alignment is described by its position and strand on the reference and read ids, sequences and base qualities are discarded for the sake of memory usage and speed.

In this section we will import the data using the *readGappedAlignmentPairs}* function, intended for paired-end data. This will create a *GappedAlignmentPairs* object that contains only the aligned reads. In order to spead up the process of importing the data, we will use the function *ScanBamParam* to load only the reads that map to chromosome 4: 

```
library(GenomicRanges)
library(Rsamtools) 

# define a filter
which <- RangesList(IRanges(1, 1351857)) 
names(which) <- "chr4"
which
param <- ScanBamParam(which=which)

# import the data
aln_chr4 <- readGappedAlignmentPairs("untreated3.bam", use.names=T, param=param)
aln_chr4
```

We have now stored our data in an object of the class *GappedAlignmentPairs*, which has been defined in the *GenomicRanges* package and does not correspond to the standard R classes. For this reason, it is useful to check the [documentation for this package](http://www.bioconductor.org/packages/release/bioc/manuals/GenomicRanges/man/GenomicRanges.pdf) in order to learn how to access our data.

**Exercise:** After having a look at the *GenomicRanges* documentation, try to answer the following questions:

* How many reads have been loaded?
* How can we access the read names? What about the strand information?
* How can we access the information for the first reads in the pair? Try to print a vector with their start coordinates.
* How many reads are properly paired?
* What is the percentage of reads that map to multiple locations?
* What information does the following command provide?
    `seqlevels(aln\_chr4)`
*Hint:* look for the *GappedAlignmentPairs* class

Since we have loaded only the reads that map to chromosome 4, we can proceed to modify the `aln_chr4` object accordingly:
```R
seqlevels(aln_chr4)="chr4"
aln_chr4
```

#### Importing the annotation
To link the alignments to their respective features, we need access to the genome annotation for the studied organism, in our case *Drosophila melanogaster*, which contains information on the coordinates of known exons, genes and transcripts. Similarly to what we encountered when loading BAM files, there is more than one way to load the annotation in R (see XX for further details). It is extremely important to pay attention to overlapping features (e.g. exons shared by multiple transcripts within the same gene), since they might end up complicating the downstream analysis (e.g. we need to make sure not to count the same read multiple times). In order to circumvent this limitation, in this practical we will use the *biomaRt* package to query Ensembl directly from R and retrieve only the necessary information:

```R
library(biomaRt)

ensembl62=useMart(host="apr2011.archive.ensembl.org", 
    biomart="ENSEMBL_MART_ENSEMBL",
    dataset="dmelanogaster_gene_ensembl")

fields=c("chromosome_name", 
    "strand", 
    "ensembl_gene_id", 
    "ensembl_exon_id", 
    "start_position", 
    "end_position", 
    "exon_chrom_start", 
    "exon_chrom_end")
annot=getBM( fields, mart=ensembl62)
```

**Exercise:** Have a look at the newly created `annot` object. What type of object is it?
*Hint:* use the function *class*

**Exercise:** In the next subsection we will calculate the overlap between the loaded BAM file and the annotation with the function *summarizeOverlaps* from the *GenomicRanges* package. What is the input required? Do we have all the necessary objects ready?
*Hint:* type `?summarizeOverlaps`

Before we proceed to calculate the counts, we need to store the annotation information in an object of the proper class:
```R
annot=GRanges(
    seqnames = Rle(paste("chr", annot$chromosome_name, sep="")),
    ranges = IRanges(start=annot$exon_chrom_start, 
                   end=annot$exon_chrom_end), 
    strand = Rle(annot$strand), 
    exon=annot$ensembl_exon_id, 
    gene=annot$ensembl_gene_id
  )
annot
class(annot)
```

#### Counting reads over known genes in R
Now that we have the alignment locations (`aln_unique` object) and the genome annotation (`annot` object), we can quantify gene expression by counting reads over all exons of a gene and summing them together. Similarly to what we encountered with htseq-count, we need to pay attention to those reads that overlap with several features.

```R
counts_chr4=summarizeOverlaps(
    annot_chr4, aln_chr4, ignore.strand=T, mode="IntersectionNotEmpty")
exon_counts_chr4=assays(counts_chr4)$counts[,1]
names(exon_counts_chr4)=elementMetadata(annot_chr4)$gene

head(exon_counts_chr4, n=15)
```

**Exercise:** So far e have obtained the number of reads overlapping each exon. How can we combine this information to obtain gene counts?\\
*Hint:* use the functions *split* and *sapply*

### Alternative approaches
In this section of the practical we have seen how to calculate the number of reads that overlap known gene models. In the two approaches evaluated here, those reads that mapped to multiple features were not considered. This is a simplification we may not want to pursue, and alternatively, there are several methods to probabilistically estimate the expression of overlapping features [Turro:2011p4448, Li:2010p4264, Trapnell:2010p3907].