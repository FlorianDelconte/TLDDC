#!/bin/bash

#input: point cloud of log
#       point cloud of defects
#       calcul TP
#       calcul Precision
#       calcul Recal


#allPointFile=$1
detectedPointsFile=$1
gtPointsFile=$2

if [ ! -f "$detectedPointsFile" ] ; then
    echo "input detected index not found!"
    exit 0
fi

if [ ! -f "$gtPointsFile" ] ; then
    echo "input ground true index not found!"
    exit 0
fi

nbDT=`wc -l <$detectedPointsFile`
nbGT=`wc -l <$gtPointsFile`

#if (( nbDT*nbGT = 0 )) ; then
#    echo "zero input"
#    exit 0
#fi

#detectedPointsFile and gtPointsFile are sorted
#sort
sort -k 1b,1 $detectedPointsFile > /tmp/sorted
sort -k 1b,1 $gtPointsFile> /tmp/sortedGT
#True positive
#calcul join
join /tmp/sorted /tmp/sortedGT >tpIndex
TP=`wc -l <tpIndex`
FN=$(( nbGT-TP ))
FP=$(( nbDT-TP ))

echo $TP $FN $FP $nbDT $nbGT >>overall
#precision
Pre=`echo "scale=3;$TP/$nbDT"|bc`
#Recall
Rec=`echo "scale=3;$TP/$nbGT"|bc`

#Acc
#F measure
Fme=`echo "scale=3;2*($Pre*$Rec)/($Pre+$Rec)"|bc`

#echo "Precision= " $Pre 
#echo "Recall= " $Rec
#echo "F Measure= " $Fme

echo -n  "&" 0$Pre" & 0"$Rec "& " 0$Fme #"\\\\"

