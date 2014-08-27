We can get the number of reads in each fastq file from the first table in the FastQC report.

We can also execute the following command if we want to automate the process:

```bash
grep '^@SRR' SRR031714_1_filt1.fastq | wc -l
    # find the lines starting with "@SRR" and count how many there are
    # 5,150,415
```

Note that being too strict in some filtering steps might lead to an important loss of information. Discarding too many reads since the beginning is generally not a good option, and one can also rely on the mapping step to discard low quality reads.
