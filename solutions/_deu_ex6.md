As we have seen in the previous section this exon had been labeled as non-testable by *DEXSeq*:
```
counts_subset=head(counts(ecs))
testable_subset=head(fData(ecs)$testable)
cbind(counts_subset, testable_subset)
```
