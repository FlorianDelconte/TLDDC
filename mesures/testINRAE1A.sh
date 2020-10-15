# !/bin/sh

#arguments for prediction
model="../models/testINRAE1a.hdf5"
tresh=127
#output file
rs=results.tex
overallVT=overall_van-tho
overallCyl=overall_Cylindrical
overallOur=overall_our
if [[ "$1" = "reset" ]]
then
    rm -f *-def-faces.id
fi
#results
resultDir=results
if [[ ! -d "$resultDir" ]]
then
    #mkdir "$resultDir"
    echo "Nothing"
fi

echo -n > $overallOur

for log in ../examples/INRAE1a/*.off
do
    logFile="$(basename $log)"
    logName=${logFile%.off}
    outputPrefix=$logName
    RM=$outputPrefix".pgm"
    groundtruthPoints="../examples/INRAE1a/"$logName"-groundtruth-points.id"
    outputUnrollPrefix=$logName
    defectUnrollPoints=$outputUnrollPrefix"-defect.id"




    if [[ ! -f $defectUnrollPoints ]]
    then
        # cylindrical fitting method
        echo "../build/segunroll -i $log -o $outputPrefix -n"
        ../build/segunroll -i $log -o $outputPrefix -n
        echo "../run/predict.py -i ${RM} ${model} ${tresh}"
        python3 ../run/predict.py ${RM} ${model} ${tresh}
        echo "../build/segToMesh -i $log -o $outputPrefix"
        ../build/segToMesh -i $log -o $outputPrefix
    fi

    echo -n $logName >>$rs
    sh ./calF.sh $defectUnrollPoints $groundtruthPoints $overallOur>> $rs
    echo "\\\\">> $rs

done

sh ./calFOverall.sh $overallOur $rs >> $rs
