# $1 complexity
# $2 complete ground truth path
# $3 assembled gtf name
# $4 folder path
# Software binaries
# $5 exonrefine binary
# $6 groupGenes groupGenes binary
# $7 gtfFilter binary

cat $2/$1_iso_ground_truth.gtf $3 > $4singleExon_pred_$1.gtf

$5 -p $4pred_$1 $4singleExon_pred_$1.gtf

$6 $4pred_$1.gtf $4pred_$1_grouped.gtf

awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "transcript_id") print $(i+1) } }' $2/$1_iso_ground_truth.gtf | sed 's/\"//g' | sed 's/\;//g' > $4groundTruth_transidList_$1.txt

awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "transcript_id") print $(i-1)"\t"$(i+1) } }' $4pred_$1_grouped.gtf | sed 's/\"//g' | sed 's/\;//g' > $4locus_transidList_$1.txt

Rscript PredictionSplitter.R $4locus_transidList_$1.txt $4groundTruth_transidList_$1.txt $4$1_

$7 -l $4$1_locusTobeKept.txt -m blackG $4pred_$1_grouped.gtf $4transcriptsToBeRemoved_$1.gtf

#awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "transcript_id") print $(i+1) } }' $4transcriptsToBeRemoved_$1.gtf | sed 's/\"//g' | sed 's/\;//g' | grep TCONS > $4removedTransList_$1.txt
awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "transcript_id") print $(i+1) } }' $4transcriptsToBeRemoved_$1.gtf | sed 's/\"//g' | sed 's/\;//g' | grep -v ENST > $4removedTransList_$1.txt

uniq $4removedTransList_$1.txt > $4removedTransListUniq_$1.txt

$7 -l $4removedTransListUniq_$1.txt -m blackT $3 $4finalAssembly_$1.gtf


rm $4singleExon_pred_$1.gtf
rm $4pred_$1.gtf
rm $4pred_$1_grouped.gtf
rm $4groundTruth_transidList_$1.txt
rm $4locus_transidList_$1.txt
rm $4$1_locusTobeKept.txt
rm $4transcriptsToBeRemoved_$1.gtf
rm $4removedTransList_$1.txt
rm $4removedTransListUniq_$1.txt
