#include <iostream>
#include <algorithm>
#include <fstream>
#include <utility>
#include <cmath>
#include <thread>
//debug
#include <stdlib.h>
#include <time.h>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"


#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Householder>

#include "DefectSegmentationUnroll.h"
#include "Regression.h"
#include "Statistic.h"
#include "IOHelper.h"
#include "MultiThreadHelper.h"
#include "UnrolledMap.h"

#include <opencv2/opencv.hpp>

using namespace DGtal;


void
DefectSegmentationUnroll::init(){
    allocate();
    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();
    convertToCcs();

    
    allocateExtra();
    
    computeEquations();
    //computeDistances();
    
}


void
DefectSegmentationUnroll::allocateExtra(){
  coefficients.resize(pointCloud.size());
}

void
DefectSegmentationUnroll::computeEquations(){
  trace.info()<<"Begin Compute Equation"<<std::endl;
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPcl(new pcl::PointCloud<pcl::PointXYZ>);
  for(int i = 0; i < pointCloud.size(); i++){
      Z3i::RealPoint p = pointCloud.at(i);
      cloudPcl->points.push_back(pcl::PointXYZ(p[0], p[1], p[2]));
  }

  pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
  kdtree.setInputCloud (cloudPcl);

  CylindricalPointOrder heightOrder;
  auto minMaxElem = std::minmax_element(myPoints.begin(), myPoints.end(), heightOrder);
  double minH= (*minMaxElem.first).height;
  double maxH = (*minMaxElem.second).height;

  int nbCores = getNumCores();
  std::vector<std::thread> ts;
  for(int i = 0; i < nbCores - 1; i++){
      ts.push_back(std::thread(&DefectSegmentationUnroll::computeEquationsMultiThread, this, i, nbCores, kdtree, minH, maxH));
  }
  computeEquationsMultiThread(nbCores - 1, nbCores, kdtree, minH, maxH);
  for(int i = 0; i < nbCores - 1; i++){
      ts[i].join();
  }
  trace.info()<<"finish eq"<<std::endl;
}

std::pair<double, double>
DefectSegmentationUnroll::computeEq(unsigned int idPoint,double searchRadius, double patchAngle,const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree){
  Z3i::RealPoint currentPoint = pointCloud.at(idPoint);
  CylindricalPoint mpCurrent = myPoints.at(idPoint);
  pcl::PointXYZ searchPoint(currentPoint[0], currentPoint[1], currentPoint[2]);

  std::vector<int> pointIdx;
  std::vector<double> radiusForEstimate;
  std::vector<double> lengthForEstimate;
  std::vector<float> pointRadiusSquaredDistance;
  std::pair<double, double> coefficient;

  if ( kdtree.radiusSearch (searchPoint, searchRadius, pointIdx, pointRadiusSquaredDistance) > 0 ){
    for (unsigned int idx = 0; idx < pointIdx.size (); ++idx){
      unsigned int foundedIndex = pointIdx.at(idx);
      Z3i::RealPoint found = pointCloud.at(foundedIndex);
      CylindricalPoint mpFound = myPoints.at(foundedIndex);
      double angleDiff = std::abs(mpFound.angle - mpCurrent.angle);
      if(angleDiff > patchAngle/2 && 2*M_PI - angleDiff > patchAngle / 2){
        continue;
      }
      radiusForEstimate.push_back(mpFound.radius);
      lengthForEstimate.push_back(mpFound.height);
    }
  }
  coefficient = Regression::PurgedByMedianlinearRegression(lengthForEstimate, radiusForEstimate);
  return coefficient;
}


void
DefectSegmentationUnroll::computeEquationsMultiThread(int threadId, int nbThread,const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, double minH, double maxH){
  std::pair<double, double> currentCoefficient;
  
  std::pair<double, double> currentCoefficient2;
  
  std::pair<double, double> currentCoefficient3;
  
  double patchAngle = arcLength / radii;
  double secondPatchAngle = patchAngle/2 + 1;
  double thirdPatchAngle = secondPatchAngle/2 + 1;
  double searchRadius = patchHeight / 2 + 1;
  double secondSearchRadius = searchRadius / 2 + 1;
  double thirdSearchRadius = secondSearchRadius / 2 + 1;
  for(unsigned int i = threadId; i < pointCloud.size();i+=nbThread){
   
    CylindricalPoint mpCurrent = myPoints.at(i);
    
    currentCoefficient = computeEq(i,searchRadius,patchAngle,kdtree);

    double estimateRadii = mpCurrent.height * currentCoefficient.first + currentCoefficient.second;


    double deltaDist= mpCurrent.radius - estimateRadii;
    //3 patch research
    //TODO : generalize
    if(deltaDist<0 ){
      currentCoefficient2 = computeEq(i,secondSearchRadius,secondPatchAngle,kdtree);
      double estimateRadii2 = mpCurrent.height * currentCoefficient2.first + currentCoefficient2.second;
      double deltaDist2= mpCurrent.radius - estimateRadii2;
      if(deltaDist2>deltaDist){
        currentCoefficient3 = computeEq(i,thirdSearchRadius,thirdPatchAngle,kdtree);
        double estimateRadii3 = mpCurrent.height * currentCoefficient3.first + currentCoefficient3.second;
        double deltaDist3= mpCurrent.radius - estimateRadii3;
        if(deltaDist3>deltaDist2){
          coefficients[i]=currentCoefficient3;
        }else{
          coefficients[i]=currentCoefficient2;
        }
      }else{
        coefficients[i]=currentCoefficient;
      }
    }else{
      coefficients[i]=currentCoefficient;
    }   
  }
}

void
DefectSegmentationUnroll::computeDistances(){
  //NOT USED, artefact from Van Tho code
}


std::vector<double>
DefectSegmentationUnroll::computeDeltaDistances(){
  std::vector<double> ddistance;
  ddistance.resize(myPoints.size());
  for(unsigned int i = 0; i < myPoints.size(); i++){
      std::pair<double, double> coeffs = coefficients[i];
      double estimateRadii = myPoints[i].height * coeffs.first + coeffs.second;
      if(coeffs.second == 0.0){
          ddistance[i] = 0;
      }else{
          ddistance[i] = myPoints[i].radius - estimateRadii;
      }
  }
  return ddistance;
}

std::vector<double>
DefectSegmentationUnroll::computeRadiusDistances(){
  std::vector<double> rdistance;
  rdistance.resize(myPoints.size());
  for(unsigned int i = 0; i < myPoints.size(); i++){  
    rdistance[i] = myPoints[i].radius;
  }
  return rdistance;
}

void
DefectSegmentationUnroll::getDefect(std::string outputFileName){
  //difference between reference distance and distance P
  //WARNING : this function need to be call after compute coefficients
  std::vector<double> DeltaDistances=computeDeltaDistances();
  //Radius of P
  std::vector<double> RadiusDistances=computeRadiusDistances();
  UnrolledMap unrolled_map(myPoints,DeltaDistances);
  cv::Mat normalizedMap=unrolled_map.computeNormalizedImageMultiScale();
  createVisuImage(outputFileName,normalizedMap);
}

//TODO : move to IOHelper
void 
DefectSegmentationUnroll::createVisuImage(std::string s,cv::Mat c){
    int grayscaleValue;
    double normalizedValue;
    int rows = c.rows;
    int cols = c.cols;
    trace.info()<<rows<<std::endl;
    trace.info()<<cols<<std::endl;
    cv::Mat grayscalemap(rows,cols,CV_8UC1,cv::Scalar(0));
    cv::Mat reliefPictures(rows, cols, CV_8UC3, cv::Scalar(110, 110, 110));

    for(unsigned int i = 0; i < rows; i++){
        for(unsigned int j = 0; j < cols; j++){
            normalizedValue=c.at<double>(i, j);
            grayscaleValue=((255/1)*(normalizedValue-1))+255;
            grayscalemap.at<uchar>(i, j) = grayscaleValue;
        }
    }
    cv::applyColorMap(grayscalemap, reliefPictures, cv::COLORMAP_JET);
    imwrite( "../unrollSurfaceOutput/"+s+".jpg", reliefPictures);
}
