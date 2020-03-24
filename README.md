#Introduction
Source code of the paper ICPR 2016 "Segmentation of defects on log surface from terrestrial Lidar data"
#Online demo
You can test the program [online]
(http://ipol-geometry.loria.fr/~kerautre/ipol_demo/TDD_IPOLDemo/)
#Installation
This guide is to install the program on linux operation system. This may work for unix-like operation system like MacOs or BSD, but we did not test. It could be compiled and run in windows operating system, however some additional software should be installed like mingw or MS C++ compiler.
## Dependencies
The program uses some C++ 11 feature, so we recommend the use of gcc 4.7 or later to compile. The program requires these libraries to be installed:
###### Boost 1.46 or later
###### DGtal version 0.9. To install DGtal see [DGtal installation] (https://github.com/DGtal-team/DGtal/blob/master/README.md)
###### PCL version 1.7 or later
###### Eigen3
###### GNU GSL 
###### CMake 2.6 or later
###### OPENCV 

## Install
```
mkdir build
cd build
cmake ..  -DDGtal_DIR=/path/to/DGtal
make
```
#User guide

=======
* Run the program
```
 ./segunroll -i ../examples/WildCherry2.off --voxelSize 5 --accRadius 100 --trackStep 20  --patchWidth 25 --patchHeight 100 --binWidth 0.01 --invertNormal false --output WC2

```
