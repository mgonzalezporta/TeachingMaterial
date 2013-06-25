# RNA-seq data analysis practical
This tutorial will illustrate how to use standalone tools, together with R and Bioconductor for the analysis of expression data obtained from high-throughput sequencing (HTS) assays. Starting from raw files as they come from a sequencing machine, we will make diagnostic plots to assess the quality of the run and filter out bad quality reads. We will then assume that these raw data have been mapped to a reference genome and we will use the aligned reads to obtain a matrix of normalised expression estimates per gene and to define novel expressed regions.


Throughout this document you will find commands that you should type or copy paste into the UNIX terminal or the R command line. There are also exercises (some of them might be challenging) that you should try to solve (the solutions are provided at the end). Whenever there is an output to your command please take some time to understand what it means.


Under the several subdirectories of `2013_RNAseq_course` there is a copy of this document, the data and the annotation files we will be using in this tutorial. We recommend that you copy all the files into your Desktop.


Along this tutorial, we will be working with a subset of a publicly available RNA-seq experiment from *Drosophila melanogaster*. This dataset is provided as a Bioconductor package ([pasilla](http://www.bioconductor.org/packages/release/data/experiment/html/pasilla.html)) and the raw sequence data can be found in the Short Read Archive ([SRP001537](http://www.ebi.ac.uk/ena/data/view/SRP001537)). More information on the experiment can be found in the original publication ([Brooks et al. 2010](http://genome.cshlp.org/content/early/2010/10/04/gr.108662.110)).

## Table of contents

1. **Dealing with raw data**
    1. [The FASTQ format](doc/_fastq.md)
    2. [Quality assessment (QA)](doc/_qa.md)
    3. [Filtering the FASTQ files](doc/_filtering_fastq.md)
    4. [De-multiplexing samples](doc/_demultiplexing.md)
    3. [Aligning reads to the genome](doc/_aligning.md)
2. **Dealing with aligned data**
    1. [The SAM/BAM format](doc/_bam.md)
    1. [Visualising aligned reads](doc/_visualising.md)
    1. [Filtering BAM files using samtools](doc/_filtering_bam.md)
    1. [Counting reads overlapping annotated genes](doc/_counting.md)
        1. With htseq-count
        1. With R
        1. Alternative approaches
    1. [Normalising counts](doc/_normalising.md)
        1. With RPKMs
        1. With DESeq
    1. [Differential gene expression](doc/_de.md)

## Software requirements
Be aware...

## Aknowledgments
This tutorial is partially based on materials developed by Angela Gonçalves, Nicolas Delhomme, Patrick Aboyoun, Simon Anders and Martin Morgan, who we would like to thank and acknowledge.