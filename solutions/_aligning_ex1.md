This is certainly a topic for discussion with no definitive answers, but here are some suggestions:

* *De novo assembly*
  * *Pros:* useful if there is no reference genome for the species in question and/or if the annotation is of poor quality.
  * *Cons:* the assembly if difficult and only the most abundant transcripts are likely to be fully assembled.

* Mapping to the genome
  * *Pros:* discovery of novel transcribed regions; *de novo* assembly of gene models for species with no annotation.
  * *Cons:* limitations intrinstic to gapped mapping; analysis of the results in order to obtain transcript expression estimates slightly more complex.

* Mapping to the transcriptome
  * *Pros:* less complexity in the mapping (reads map contiguously); output easily interpretable; faster than mapping to the genome.
  * *Cons:* only for annotated regions; multiple transcripts for the same gene (reads might map equally well to these because of shared sequence).