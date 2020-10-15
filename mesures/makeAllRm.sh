# !/bin/sh
#parameters

#output file
#rs=results.tex

if [[ "$1" = "reset" ]]
then
    rm -f *-def-faces.id
fi

for log in ../examples/INRAE1a/*.off
do
    logFile="$(basename $log)"
    logName=${logFile%.off}
    outputPrefix=$logName
    echo $logName


      if [[ ! -f $defectUnrollPoints ]]
      then
          echo "../build/segunroll -i $log -o  $logName -n --intensityPerCm 40 --grayscaleOrigin 0 "
          ../build/segunroll -i $log -o $logName -n --intensityPerCm 40 --grayscaleOrigin 0
      fi



done
