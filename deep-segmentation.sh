#!/bin/bash
##################################################
#centrerline parameters FIXED
binWidth=0.01
accRadius=200
voxelSize=5
trackStep=20
#input path to the mesh
pathToLog=$1
#name of the mesh
temp="$(basename $pathToLog)"
logName="${temp%.*}"
#path to read relief map
inputReliefMap="output/RM_output/"$logName".png"
#output path to write segmentation
outputSegmentation="input/SEG_input/"
##################################################
#write reliefmap in ../output/RM_output
echo "build/./segunroll -i $pathToLog --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep --binWidth $binWidth --invertNormal --output  $logName"
build/./segunroll -i $pathToLog --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep  --binWidth $binWidth --invertNormal --output $logName
#make prediction + treshold
echo "python3 DeepLearningSegmentation/predictOne.py" $inputReliefMap $outputSegmentation
python3 DeepLearningSegmentation/predictOne.py $inputReliefMap $outputSegmentation
#write result on mesh
echo "build/./segunroll -i $pathToLog --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep --binWidth $binWidth --invertNormal --output  $logName -d true"
build/./segunroll -i $pathToLog --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep --binWidth $binWidth --invertNormal --output  $logName -d true
