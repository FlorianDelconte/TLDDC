# !/bin/sh
#parameters
binWidth=0.01
accRadius=200
voxelSize=5
trackStep=20


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

#empty result latex file
echo -n > $rs
#empty overal for vt method
echo -n > $overallVT
#empty overal for cyl method
echo -n > $overallCyl
#empty overall for our method
echo -n > $overallOur

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
    associateFileCyl=$outputCylPrefix"-asso.off"
    defectFaceIdsCyl=$outputCylPrefix"-def-faces.id"
    # for unroll method
    outputUnrollPrefix=$logName"Unroll"
    defectUnrollPoints=$outputUnrollPrefix"-defect.id"
    associateFileUnroll=$outputUnrollPrefix"-asso.off"
    defectFaceIdsUnroll=$outputUnrollPrefix"-def-faces.id"

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

    if [[ ! -f $defectUnrollPoints ]]
    then
        # cylindrical fitting method
        echo "../../build/segunroll -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
              --binWidth $binWidth --invertNormal --output  $outputUnrollPrefix"
        ../../build/segunroll -i $log --voxelSize $voxelSize --accRadius $accRadius --trackStep $trackStep \
        --binWidth $binWidth --invertNormal --output  $outputUnrollPrefix
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

    if [[ ! -f $associateFileCyl ]]
    then
        echo "../../build/colorizeMesh -i $log -r $groundtruth -t $defectFaceIdsCyl -o $associateFileCyl"
        ../../build/colorizeMesh -i $log -r $groundtruth -t $defectFaceIdsCyl -o $associateFileCyl
    fi

    if [[ ! -f $associateFileUnroll ]]
    then
        echo "../../build/colorizeMesh -i $log -r $groundtruth -t $defectFaceIdsUnroll -o $associateFileUnroll"
        ../../build/colorizeMesh -i $log -r $groundtruth -t $defectFaceIdsUnroll -o $associateFileUnroll
    fi

    #write to latex file
    echo -n $logName >>$rs
    sh ./calF.sh $defectPoints $groundtruthPoints $overallVT>> $rs
    echo -n "&">>$rs
    sh ./calF.sh $defectCylPoints $groundtruthPoints $overallCyl >> $rs
    echo -n "&">>$rs
    sh ./calF.sh $defectUnrollPoints $groundtruthPoints $overallOur>> $rs
    echo "\\\\">> $rs

done

# write overall to tex file
sh ./calFOverall.sh $overallVT $rs >> $rs
sh ./calFOverall.sh $overallCyl $rs >> $rs
sh ./calFOverall.sh $overallOur $rs >> $rs
