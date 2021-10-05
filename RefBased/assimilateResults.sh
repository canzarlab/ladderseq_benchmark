#!/bin/bash

#### Input:
## 1. gffCompare results or gtfCompare results (input values : gff ot gtf)
## 2. concat, original or merged
## 3. output file name

file=$3.txt
fileComplexity_Recall=$3"CompRe".txt
filePrecision=$3"Pre".txt
echo $file
if [ $1 = "gff" ];then
echo "Assimilating results from gffCompare."
echo "Complexity Recall" > $fileComplexity_Recall
echo "Precision bandProb" > $filePrecision
for f in $2gffComp*.stats ; do
   awk 'NR == 14 {print FILENAME"	",$4"	"}' $f | awk '{print substr($1,length($1)-28,28),"\t",$2}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1'  >> $fileComplexity_Recall
done
for f in $2CompleteGround*.stats ; do
   awk 'NR == 14 {printf $6"	"}' $f >> $filePrecision
   echo $4 >> $filePrecision

done
pr -m -t $fileComplexity_Recall $filePrecision >> $file
rm $fileComplexity_Recall
rm $filePrecision
else
echo "Assimilating results from gtfCompare."
echo "Complexity Recall" > $fileComplexity_Recall
echo "Precision bandProb" > $filePrecision
for f in $2/*gtfComp.stat ; do
   awk '{if(NR == 3)printf FILENAME"	"};{if(NR == 3)printf $4"	"}' $f  | awk '{print substr($1,length($1)-28,28),"\t",$4}' | awk '{gsub(/[^[:digit:]]/, "", $1)}1' >> $fileComplexity_Recall
done
for f in $2/*CompleteGround.stat ; do
   awk '{if(NR == 8)printf $4"	"}' $f > $filePrecision
   echo $4 > $filePrecision
done
pr -m -t $fileComplexity_Recall $filePrecision >> $file
rm $fileComplexity_Recall
rm $filePrecision
fi
