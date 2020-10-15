#!/bin/bash
modelRepo="../models/"
K_number=5
tresh=127
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


groundtruth_sufix="-groundtruth-points.id"

#For each fold number
for number in {1..5}
do
  rs=results$number.tex
  echo -n > $rs
  overallOur=overall_our$number
  echo -n > $overallOur
done
logPath="../examples/INRAE1a/"
for i in {0..1}
do
  #loop on log
  for log in $logPath*.off
  do
      logFile="$(basename $log)"
      logName=${logFile%.off}
      outputPrefix=$logName
      RM=$outputPrefix".pgm"
      groundtruthPoints=$logPath$logName$groundtruth_sufix
      outputUnrollPrefix=$logName
      defectUnrollPoints=$outputUnrollPrefix"-defect.id"

      modelNumber=$(./searchModelNumber.sh $logName)
      modelName=$modelRepo"k"$modelNumber".hdf5"

      if [[ ! -f $defectUnrollPoints ]]
      then
          # cylindrical fitting method
          echo "../build/segunroll -i $log -o $outputPrefix -n"
          ../build/segunroll -i $log -o $outputPrefix -n
          echo "../run/predict.py -i ${RM} ${modelName} ${tresh}"
          python3 ../run/predict.py ${RM} ${modelName} ${tresh}
          echo "../build/segToMesh -i $log -o $outputPrefix"
          ../build/segToMesh -i $log -o $outputPrefix
      fi

      #write to latex file
      echo -n "K"$modelNumber" &">>results$modelNumber.tex
      echo -n $logName >>results$modelNumber.tex
      sh ./calF.sh $defectUnrollPoints $groundtruthPoints overall_our$modelNumber>> results$modelNumber.tex
      echo "\\\\">> results$modelNumber.tex

  done
  logPath="../examples/INRAE1b/"
done
#loop for overall
for i in {1..5}
do
  sh ./calFOverall.sh overall_our$i results$i.tex>> results$i.tex
done
