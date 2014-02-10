#!/bin/bash

# IMPORTANT: always leave two empty spaces at the end of each md file,
#			 otherwise pandonc doesn't distinguish the sections properly
#
# human readable command:
# pandoc *.md 
#     -o ../pdf/practical.pdf 
#     --toc 
#     --variable title:"RNA-seq data analysis practical - San Michele all'Adige (Trento), Italy"
#     --variable date:"2013/07/02"
#     --variable author:"Mar Gonzàlez-Porta"
#     --variable links-as-notes 
#     --variable linkcolor:black 
#     --variable urlcolor:black 
#     --variable geometry:margin=3cm

pandoc *.md -o ../pdf/practical.pdf --toc --variable title:"RNA-seq data analysis practical - Agricultural Omics, EBI" --variable date:"2014/02/19" --variable author:"Mar Gonzàlez-Porta" --variable links-as-notes --variable linkcolor:black --variable urlcolor:black --variable geometry:margin=3cm
