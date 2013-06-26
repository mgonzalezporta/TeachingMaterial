* Try the following accessor methods:
```rconsole
length(aln_chr4)
head(names(aln_chr4))
seqnames(aln_chr4)
strand(aln_chr4)
first(aln_chr4)
last(aln_chr4)
head(start(first(aln_chr4)))
```

* How many reads are properly paired?
```rconsole
table(isProperPair(aln_chr4))
```

* What is the percentage of reads that map to multiple locations?

	```rconsole
	t=table(names(aln_chr4))
	head(t[t>1])
	length(t[t>1])/length(t)*100
	
	# let us inspect one of the multireads
	aln_chr4[names(aln_chr4)=="SRR031714.1029411"]
	```

* What information does the command `seqlevels(aln_chr4)` provide?

  It contains information on the chromosome names available in the BAM file.

