# !/bin/sh
#parameters
binWidth=0.01
accRadius=200
voxelSize=5
trackStep=25

#output file
#rs=results.tex

if [[ "$1" = "reset" ]]
then
    rm -f *-def-faces.id
fi

for log in ../*off
do
    logFile=${log:3}
    logName=${logFile%.off}
    outputPrefix=$logName
    echo $logFile


      if [[ ! -f $defectUnrollPoints ]]
      then
          echo "../../build/segunroll -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
              --binWidth $binWidth --invertNormal --output  $outputUnrollPrefix"
          ../../build/segunroll -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
          --binWidth $binWidth --invertNormal --output  $logName
      fi



done
