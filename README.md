# RNA-seq data analysis practical
This tutorial will illustrate how to use standalone tools, together with R and Bioconductor for the analysis of RNA-seq data. Keep in mind that this is a rapidly evolving field and that this document is not intended as a review of the many tools available to perform each step; instead, we will cover one of the many existing workflows to analyse this type of data.

We will be working with a subset of a publicly available dataset from *Drosophila melanogaster*, which is available both in the Short Read archive ([SRP001537](http://www.ebi.ac.uk/ena/data/view/SRP001537) - raw data) and in Bioconductor ([pasilla package](http://www.bioconductor.org/packages/release/data/experiment/html/pasilla.html) - processed data). For more information about this dataset please refer to the original publication ([Brooks et al. 2010](http://genome.cshlp.org/content/early/2010/10/04/gr.108662.110)).

The tools and R packages that we will be using during the practical are listed below (see [Software requirements](https://github.com/mgonzalezporta/TeachingMaterial#software-requirements)) and the necessary data files can be found [here](http://www.ebi.ac.uk/~mar/courses/RNAseq.tar.gz). After dowloading and uncompressing the `tar.gz` file, you should have the following directory structure in yuor computer:

```
RNAseq
|-- reference               # reference info (e.g. genome sequence and annotation)
`-- data
    |-- raw                 # raw data: fastq files
    |-- mapped              # mapped data: BAM files
    `-- demultiplexing      # extra fastq files for the demultiplexing section

```

## Table of contents

1. **Dealing with raw data**
    1. [The FASTQ format](doc/_fastq.md)
    2. [Quality assessment (QA)](doc/_qa.md)
    3. [Filtering FASTQ files](doc/_filtering_fastq.md)
    4. [De-multiplexing samples](doc/_demultiplexing.md)
    3. [Aligning reads to the genome](doc/_aligning.md)
2. **Dealing with aligned data**
    1. [The SAM/BAM format](doc/_bam.md)
    1. [Visualising aligned reads](doc/_visualising.md)
    1. [Filtering BAM files](doc/_filtering_bam.md)
    1. [Counting reads overlapping annotated genes](doc/_counting.md)
        * With htseq-count
        * With R
        * Alternative approaches
    1. [Normalising counts](doc/_normalising.md)
        * With RPKMs
        * With DESeq
    1. [Differential gene expression](doc/_de.md)

## Software requirements
* Standalone tools:
  * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [PRINSEQ](http://prinseq.sourceforge.net/)
  * [eautils](https://code.google.com/p/ea-utils/)
  * [samtools](http://sourceforge.net/projects/samtools/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)
  * [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

* Bioconductor packages:
  * [GenomicRanges](http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
  * [Rsamtools](http://www.bioconductor.org/packages/release/bioc/html/Rsamtools.html)
  * [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html)
  * [pasilla](http://www.bioconductor.org/packages/release/data/experiment/html/pasilla.html)
  * [DESeq](http://www.bioconductor.org/packages/release/bioc/html/DESeq.html)

*Note: depending on the topics covered in the course we might end up not using some of the above listed tools.*

## Other resources

### Course data
* [Complete course data, including command outputs and R sessions](http://www.ebi.ac.uk/~mar/courses/RNAseq_all.tar.gz)

### Tutorials
* [Course materials available at the Bioconductor website](http://www.bioconductor.org/help/course-materials/)
* [Online training resources at the EBI website](http://www.ebi.ac.uk/training/online/course-list?topic%5B%5D=13&views_exposed_form_focused_field=)
* [R and Bioconductor tutorial by Thomas Girke](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual)
* Do not forget to check the documentation for the packages used in the practical!

### Cheat sheets
* [R reference card](http://cran.r-project.org/doc/contrib/Short-refcard.pdf)
* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)


## Aknowledgments
This tutorial has been inspired on material developed by Ângela Gonçalves, Nicolas Delhomme, Simon Anders and Martin Morgan, who I would like to thank and acknowledge. Special thanks must go to Ângela Gonçalves, with whom I started teaching, and Gabriella Rustici, for always finding a way to organise a new course.