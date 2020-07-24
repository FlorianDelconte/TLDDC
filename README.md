## TLDDC
<!--This code use a previous work (https://github.com/vanthonguyen/treelogdefectsegmentation). The previous work first compute a centerline from a tree mesh (by accumulation in 3D voxels) then an estimation of the surface radius is done by the studies of rectangular patch around vertices mesh. This allows the computation of delta distance for each vertices. Finally a treshold  (Rosin) is compute to distinguish wich vertices is a part of defect or not. To know more about the previous work, there is an online demo : http://ipol-geometry.loria.fr/~kerautre/ipol_demo/TDD_IPOLDemo/ and there is an article :  ICPR 2016 "Segmentation of defects on log surface from terrestrial Lidar data".

Here, we want to explore an other use of the centerline. The main idea, is using centerline to unroll tree mesh, then work on 2D image with an intensity able to represent defect shape. This image may be usefull to labeling task, segmentation task, and clasification task. -->
Repository of work submitted to ICPR 2020:

**Tree Defect Segmentation usingGeometric Features and CNN**

the run of the methods can be done without any installation with this online demonstration [here](http://kerautret.github.io/TLDDC)

This repository allow users to segment wood log surface defects on meshs. First, we compute a relief map based on height and circumference of the log. Then, segmentation is done on the relief map by our trained CNN model. Finally, segmentation result is send on mesh. Below, there is a pipeline of the project :  

![alt text](pipeline.png?raw=true "Pipeline")

## Dependencies
The program uses some C++ 11 feature, so we recommend the use of gcc 4.7 or later to compile. The program requires these libraries to be installed :
###### Boost 1.46 or later
###### DGtal version 0.9.4 To install DGtal see [DGtal installation] (https://github.com/DGtal-team/DGtal/blob/master/README.md)
###### PCL version 1.7 or later
###### Eigen3
###### GNU GSL
###### CMake 2.6 or later
###### OPENCV
Morover if you want to use our trained CNN model to make segmentation, you need theses dependencies to be installed :
###### tensorflow2.2
###### tensorflow-addons
## Install
```
mkdir build
cd build
cmake ..  -DDGtal_DIR=/path/to/DGtal
make
```
## Run the program
Generate a relief map :
```
./segunroll -i ../examples/INRAE1a/WildCherry2.off --voxelSize 5 --accRadius 200 --trackStep 20  --patchWidth 25 --patchHeight 100 --binWidth 0.01 --invertNormal false
```
Defect segmentation on surface log using our trained model :
```
cd runDemo
./deep-segmentation.sh ../examples/INRAE1a/Beech.off
```
