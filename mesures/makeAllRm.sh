# !/bin/sh
#parameters

#output file
#rs=results.tex

pathToMesh=$1

for log in $pathToMesh*.off
do
    logFile="$(basename $log)"
    logName=${logFile%.off}
    outputPrefix=$logName
    echo $logName


      if [[ ! -f $defectUnrollPoints ]]
      then
          echo "../build/segunroll -i $log -o  $logName -n "
          ../build/segunroll -i $log -o $logName -n
      fi



done
