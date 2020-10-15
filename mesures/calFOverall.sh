#!/bin/bash

#input: point cloud of log
#       point cloud of defects
#       calcul TP
#       calcul Precision
#       calcul Recal


#allPointFile=$1
iFile=$1
oFile=$2


if [ ! -f "$iFile" ] ; then
    echo "input file !"
    exit 0
fi

nbDTPatch=`awk '{sum+=$4} END {print sum}' $iFile`
#nbDTCyl=`awk 'NR%2==0{sum+=$4} END {print sum}' $iFile`
nbGT=`awk '{sum+=$5} END {print sum}' $iFile`
TPPatch=`awk '{sum+=$1} END {print sum}' $iFile`
#TPCyl=`awk 'NR%2==0{sum+=$1} END {print sum}' $iFile`


#echo $TP $FN $FP $nbDT $nbGT >>../overall
#precision
PrePatch=`echo "scale=3;$TPPatch/$nbDTPatch"|bc`
#Recall
RecPatch=`echo "scale=3;$TPPatch/$nbGT"|bc`

#Acc
#F measure
FmePatch=`echo "scale=3;2*($PrePatch*$RecPatch)/($PrePatch+$RecPatch)"|bc`

#echo "Precision= " $Pre
#echo "Recall= " $Rec
#echo "F Measure= " $Fme

echo "\hline" >> $oFile
echo -n  "Overall& 0"$PrePatch" & 0"$RecPatch "&  0"$FmePatch  #"\\\\"

#PreCyl=`echo "scale=3;$TPCyl/$nbDTCyl"|bc`
#Recall
#RecCyl=`echo "scale=3;$TPCyl/$nbGT"|bc`

#Acc
#F measure
#FmeCyl=`echo "scale=3;2*($PreCyl*$RecCyl)/($PreCyl+$RecCyl)"|bc`

#echo "Precision= " $Pre
#echo "Recall= " $Rec
#echo "F Measure= " $Fme
#echo -n  "&" 0$PreCyl" & 0"$RecCyl"& " 0$FmeCyl "\\\\" >> $oFile
