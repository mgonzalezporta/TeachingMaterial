## Filtering BAM
Samtools can also be used to further modify and/or subset BAM files. For example, some tools will require that the reads are sorted by coordinate or by name. In addition, and similarly to what we did with the fastq files, one might consider to discard reads with a low alignment quality (including reads that align to several locations in the genome), or in the case of paired-end data, discard reads that are not properly paired.

**Exercise:** Try to answer to the following questions using samtools and the information provided in [the documentation](http://samtools.sourceforge.net/samtools.shtml):

* How many reads are properly paired?
* By default, TopHat creates BAM files where the reads are sorted by coordinate. How would you sort the properly paired reads by name instead? Save the output in a new BAM file called `untreated3_paired.bam`, which we will use later on during the practical.
* Which percentage of those properly paired reads map uniquely?
  *Hint:* Have a look at the options for `samtools sort` and `samtools view`.

[Solution](../solutions/_filtering_bam.md)