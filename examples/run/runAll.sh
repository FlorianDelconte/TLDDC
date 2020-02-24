# !/bin/sh

#parameters
binWidth=0.01
accRadius=200
voxelSize=5
trackStep=25

#output file 
rs=results.tex

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

#empty result latex file
echo -n > $rs

for log in ../*off 
do
    logFile=${log:3}
    logName=${logFile%.off}
    outputPrefix=$logName
    groundtruth="../"$logName"-groundtruth.id"
    groundtruthPoints="../"$logName"-groundtruth-points.id"
    defectFaceIds=$outputPrefix"-def-faces.id"
    associateFile=$logName"-asso.off"

    defectPoints=$logName"-defect.id"


    # for cylindrical method
    outputCylPrefix=$logName"Cyl"

    defectCylPoints=$outputCylPrefix"-defect.id"

    if [[ ! -f $defectFaceIds ]]
    then
        #main command
        echo "../../build/segmentation -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep --binWidth $binWidth --invertNormal --output  $outputPrefix"
        ../../build/segmentation -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
        --binWidth $binWidth --invertNormal --output  $outputPrefix
    fi

    if [[ ! -f $defectCylPoints ]]
    then
        # cylindrical fitting method
        echo "../../build/segcyl -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
            --binWidth $binWidth --invertNormal --output  $outputCylPrefix"
        ../../build/segcyl -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
        --binWidth $binWidth --invertNormal --output  $outputCylPrefix
    fi
    #associate color of detected defects and ground truth
    #Yellow = defects ^ ground truth
    #Green = defects - ground truth
    #Red = ground truth - defects
    if [[ ! -f $associateFile ]]
    then
        echo "../../build/colorizeMesh -i $log -r $groundtruth -t $defectFaceIds -o $associateFile"  
        ../../build/colorizeMesh -i $log -r $groundtruth -t $defectFaceIds -o $associateFile
    fi

    #write to latex file
    echo -n $logName >>$rs
    calF.sh $defectPoints $groundtruthPoints >> $rs
    echo -n "&">>$rs
    calF.sh $defectCylPoints $groundtruthPoints >> $rs
    echo "\\\\">> $rs

done

# write overall to tex file
calFOverall.sh overall $rs
