## TLDDC
<!--This code use a previous work (https://github.com/vanthonguyen/treelogdefectsegmentation). The previous work first compute a centerline from a tree mesh (by accumulation in 3D voxels) then an estimation of the surface radius is done by the studies of rectangular patch around vertices mesh. This allows the computation of delta distance for each vertices. Finally a treshold  (Rosin) is compute to distinguish wich vertices is a part of defect or not. To know more about the previous work, there is an online demo : http://ipol-geometry.loria.fr/~kerautre/ipol_demo/TDD_IPOLDemo/ and there is an article :  ICPR 2016 "Segmentation of defects on log surface from terrestrial Lidar data".

Here, we want to explore an other use of the centerline. The main idea, is using centerline to unroll tree mesh, then work on 2D image with an intensity able to represent defect shape. This image may be usefull to labeling task, segmentation task, and clasification task. -->
Repository of work submitted to RRPR 2020:

**Tree Defect Segmentation usingGeometric Features and CNN**

the run of the methods can be done without any installation with this online demonstration [here](http://kerautret.github.io/TLDDC)

This repository allow users to segment wood log surface defects on meshs. First, we compute a relief map based on height and circumference of the log. Then, segmentation is done on the relief map by our trained CNN model. Finally, segmentation result is send on mesh. Below, there is a pipeline of the project :  

![alt text](pipeline.png?raw=true "Pipeline")

## Dependencies
The program uses some C++ 14 feature, so we recommend the use of gcc 4.7 or later to compile. The program requires these libraries to be installed (we add instruction to install for ubuntu and debian users):
###### DGtal version 1.0.0 or later, see [DGtal installation] (https://github.com/DGtal-team/DGtal/blob/master/README.md)
###### Eigen3, using apt manager :
``` sudo apt install libeigen3-dev ```
###### GNU GSL, using apt manager :
``` sudo apt install libgsl-dev ```
###### PCL version 1.7 or later, using apt manager :
``` sudo apt install libpcl-dev ```

Morover if you want to use our trained CNN model to make segmentation, you need theses dependencies to be installed (virtual environnement is recommanded):
###### tensorflow2.2, see [TensorFlow installation] (https://www.tensorflow.org/install/pip)
###### tensorflow-addons, see [TensorFlow addons installation] (https://www.tensorflow.org/addons/overview)
###### openCV for python, see [OpenCV installation] (https://pypi.org/project/opencv-python/)
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
./segunroll -i *PathToMesh* -n
```
Generate all relief map from INRAE1a directory :
```
cd mesures/
./makeAllRm.sh
```
