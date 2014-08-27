#!/bin/bash

# IMPORTANT: always leave two empty lines at the end of each md file,
#			 otherwise pandonc doesn't distinguish the sections properly
#
# human readable command:
# pandoc *.md 
#     -o ../pdf/practical.pdf 
#     --toc 
#     --variable title:"RNA-seq data analysis practical"
#     --variable date:"2013/07/02"
#     --variable author:"Mar Gonzàlez-Porta"
#     --variable links-as-notes 
#     --variable linkcolor:black 
#     --variable urlcolor:black 
#     --variable geometry:margin=3cm

pandoc *.md -o ../pdf/practical.pdf --toc --variable title:"RNA-seq data analysis practical" --variable date:"2014/08/27" --variable author:"Mar Gonzàlez-Porta" --variable links-as-notes --variable linkcolor:black --variable urlcolor:black --variable geometry:margin=3cm
