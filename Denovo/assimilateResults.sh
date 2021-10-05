#!/bin/bash

#### Input:
## 1. Path to blat folder with different complexity forders within
## 2. Output File Path {print substr(FILENAME,1,4),

file=$2.txt
fileComplexity_Recall_Ladder=$2"CompRe_ladder".txt
fileComplexity_Recall_original=$2"CompRe_original".txt


echo "Assimilating results from blat."
echo "Complexity Recall Recall_num PercentageLength  bandProb" > $fileComplexity_Recall_Ladder
echo "Complexity Recall Recall_num PercentageLength  bandProb" > $fileComplexity_Recall_original

## Assimilating results for ladder assembly
for f in $1complexity[0-9]/results.txt ; do
   awk 'NR%2==0 {print FILENAME,"\t",$3,"\t",$6,"\t",$4,"\t",$5," "}' $f | awk '{print substr($1,length($1)-28,28),"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall_Ladder
done
awk 'NR%2==0 {print FILENAME"\t",$3,"\t",$6,"\t",$4,"\t",$5," "}' $1complexity10/results.txt | awk '{print substr($1,length($1)-28,28),"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall_Ladder
awk 'NR%2==0 {print FILENAME"\t",$3,"\t",$6,"\t",$4,"\t",$5," "}' $1complexity11/results.txt | awk '{print substr($1,length($1)-28,28),"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall_Ladder

## Assimilating results for original trinity assembly
for f in $1complexity[0-9]/results.txt ; do
   awk 'NR!=1 && NR%2==1 {print FILENAME"\t",$3,"\t",$6,"\t",$4,"\t",$5," "}' $f | awk '{print substr($1,length($1)-28,28),"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall_original
done
awk 'NR!=1 && NR%2==1 {print FILENAME"\t",$3,"\t",$6,"\t",$4,"\t",$5," "}' $1complexity10/results.txt | awk '{print substr($1,length($1)-28,28),"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall_original
awk 'NR!=1 && NR%2==1 {print FILENAME"\t",$3,"\t",$6,"\t",$4,"\t",$5," "}' $1complexity11/results.txt | awk '{print substr($1,length($1)-28,28),"\t",$2,"\t",$3,"\t",$4,"\t",$5}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall_original
