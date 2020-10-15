NameOfAssociatedFile="../models/KFoldAssociation"
exampleName=$1
modelNumberFind=-1
while IFS=' ' read -r col1 col2
do
  if [ $exampleName == $col1 ]
  then
    modelNumberFind=$col2
  fi

done < $NameOfAssociatedFile

echo $modelNumberFind
