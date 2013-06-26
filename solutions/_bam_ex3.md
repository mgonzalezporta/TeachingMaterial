
* How many reads map in total?
```bash
samtools view untreated3.bam | wc -l
	# 33,672,534 /2 = 16,836,267
```

* How many reads map to each chromosome?
```bash
samtools view untreated3.bam | awk '{print $3}' | sort | uniq -c | sort -nr > chr.txt
less chr.txt
```

* How many different mapping qualities are represented in the BAM file. What do they represent?
```bash
samtools view untreated3.bam | awk '{print $5}' | sort | uniq -c | sort -nr
```

  TopHat (and Bowtie) does not provide as much information encoded in the mapping quality as other software (e.g. BWA). Still, those are usually the values reported:
  * `255`: unique mapping
  * `3`: the read maps to 2 locations in the target
  * `2`: the read maps to 3 locations
  * `1`: the read maps to 4-9 locations
  * `0`: the read maps to 10 or more locations
  
* How many different alignment flags can you find in the BAM file?
```bash
samtools view tophat_out/accepted_hits.bam | awk '{print $2}' | sort | uniq -c | sort -nr
```

* Try to print the unique CIGAR strings for the first 300 reads. What is their meaning?
```bash
samtools view untreated3.bam | head -n 300 | awk '{print $6}' | sort -u
```

  * `36M` - 36 mapped positions
  * `19M260N11M` - 19 mapped positions, followed by a gap of 260 nucleotides, followed by 11 mapped positions
  * etc.

  Note that `M` indicates that the position could be mapped, but it does not specify if it was an exact match or a mismatch.