## The FASTQ format
The nucleotide sequences and qualities of the short reads produced in a sequencing experiment are commonly stored in a plain text file using the FASTQ format. In the `data/raw` directory, you will find two fastq files, which contain information about the short reads obtained from one of the samples in the *Drosophila melanogaster* experiment.

**Exercise:** Why do we have two fastq files for this given sample?
[Solution](../solutions/_fastq_ex1.md)

To confirm that we are working with a fastq file and to get an idea of how this format looks like we can print the first lines of our files by typing this into the terminal:

```bash
zcat SRR031714_1.fastq.gz | head
zcat SRR031714_2.fastq.gz | head
```

**Exercise:** How many lines are used to represent a read in the fastq file? Which information do they contain?

[Solution](../solutions/_fastq_ex2.md)

**Exercise:** How many reads are there in each file? Do both files contain the same number of reads? Is that what we would expect?

[Solution](../solutions/_fastq_ex3.md)