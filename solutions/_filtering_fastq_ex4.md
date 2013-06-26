We can get the number of reads in each fastq file from the first table in the FastQC report:

```
SRR031709.fastq: 3,812,809		
SRR031709_filt1.fastq: 3,668,330
SRR031709_filt2.fastq: 2,771,725
SRR031709_filt3.fastq: 3,034,517
```

 We can also execute the following command if we want to automate the process:

```bash
grep '^@SRR' SRR031709.fastq | wc -l
    # find the lines starting with "@" and count how many there are
    # 3,812,809
```

Note that being too strict in some filtering steps might lead to an important loss of information. Discarding too many reads since the beginning is generally not a good option, and one can also rely on the mapping step to discard low quality reads.