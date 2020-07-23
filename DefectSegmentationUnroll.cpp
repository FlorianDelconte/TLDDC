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
#include "ImageAnalyser.h"

#include <opencv2/opencv.hpp>

#include <Python.h>
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

}


void
DefectSegmentationUnroll::allocateExtra(){
  coefficients.resize(pointCloud.size());
  ind_Patches.resize(pointCloud.size());
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
  std::vector<double> angleForEstimate;
  std::vector<double> lengthForEstimate;
  std::vector<unsigned int> indForEstimate;
  std::vector<float> pointRadiusSquaredDistance;


  if ( kdtree.radiusSearch (searchPoint, searchRadius, pointIdx, pointRadiusSquaredDistance) > 0 ){
    for (unsigned int idx = 0; idx < pointIdx.size (); ++idx){
      unsigned int foundedIndex = pointIdx.at(idx);

      CylindricalPoint mpFound = myPoints.at(foundedIndex);
      double angleDiff = std::abs(mpFound.angle - mpCurrent.angle);
      if(angleDiff > patchAngle/2 && 2*M_PI - angleDiff > patchAngle / 2){
        continue;
      }
      //fill the patches vector
      angleForEstimate.push_back(mpFound.angle);
      radiusForEstimate.push_back(mpFound.radius);
      lengthForEstimate.push_back(mpFound.height);
      indForEstimate.push_back(foundedIndex);
    }
  }
  struct coefs c= Regression::PurgedlinearRegression(lengthForEstimate, radiusForEstimate, angleForEstimate,indForEstimate);
  std::pair<double, double> coef;

  coef=c.coefficients;
  ind_Patches.at(idPoint)=c.ind_p;
  return coef;
}


void
DefectSegmentationUnroll::computeEquationsMultiThread(int threadId, int nbThread,const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, double minH, double maxH){
  std::pair<double, double> currentCoefficient;

  double patchAngle = arcLength / radii;

  double searchRadius = patchHeight / 2 + 1;

  for(unsigned int i = threadId; i < pointCloud.size();i+=nbThread){
    CylindricalPoint mpCurrent = myPoints.at(i);
    currentCoefficient = computeEq(i,searchRadius,patchAngle,kdtree);
    coefficients[i]=currentCoefficient;
  }
}

void
DefectSegmentationUnroll::computeDistances(){
  //NOT USED, artefact from Van Tho code
}



void
DefectSegmentationUnroll::computeDeltaDistances(){

  for(unsigned int i = 0; i < myPoints.size(); i++){
      std::pair<double, double> coeffs = coefficients[i];
      double estimateRadii = myPoints[i].height * coeffs.first + coeffs.second;
      if(coeffs.second == 0.0){
          distances[i] = 0;
      }else{
          distances[i] = myPoints[i].radius - estimateRadii;
      }

  }
}

void
DefectSegmentationUnroll::computeRadiusDistances(){

  for(unsigned int i = 0; i < myPoints.size(); i++){
    distances[i] = myPoints[i].radius;
  }

}


std::vector<unsigned int>
DefectSegmentationUnroll::getDefect(std::string outputFileName,std::string gtName, bool readSegmentation ){
  std::vector<unsigned int> concatPointInPixels;
  //compute vector of distances like distance=deltadistance
  computeDeltaDistances();
  //Compute vector of distance like distance=radius
  //computeRadiusDistances();
  //Compute vector of distancelike distance=deltadistance filtred by rosin
  //computeDeltaDistancesRosin();
  //Construct Unrolled_map with a vectorof distance (vector size = point cloud size)
  UnrolledMap unrolled_map(myPoints,distances);
  //compute the normalized image
  /*unrolled_map.computeNormalizedImage(1);
  unrolled_map.computeRGBImage();
  imwrite( "../../unrollSurfaceOutput/"+outputFileName+"RGB.png", unrolled_map.getImage());
  //unrolled_map.computeGRAYImage();*/
  //imwrite( "../../unrollSurfaceOutput/"+outputFileName+"GRAY.png", unrolled_map.getImage());
  //compute an rgb image from normalized image
  unrolled_map.computeNormalizedImageMultiScale();
  unrolled_map.computeGRAYImage();
  imwrite(outputFileName+".png", unrolled_map.getImage());
  std::cout << outputFileName+".png"<< std::endl;
  //std::cout << "size :"<<unrolled_map.getImage().size()<< std::endl;

  //uncomment to test drawing mesh
  if(readSegmentation){
    cv::Mat drawn;

    drawn = cv::imread(outputFileName+"SEG.png",0);
    if (drawn.empty())
    {
      std::cout << "!!! Failed imread(): drawinf input not found" << std::endl;
    }
    trace.info()<<outputFileName<<"SEG.png"<<std::endl;

    std::vector<unsigned int>  currentPointsInPixels;
    unsigned int IndP;
    int currentIntensity;
    for (int i=0; i < drawn.rows; ++i){
      for (int j=0; j < drawn.cols; ++j){
        currentIntensity=(int)drawn.at<uchar>(i, j);
        if(currentIntensity>0){
          currentPointsInPixels=unrolled_map.getPointsUnrolled_surface(i,j);

          for(std::vector<unsigned int>::iterator it = std::begin(currentPointsInPixels); it != std::end(currentPointsInPixels); ++it) {
            IndP=*it;
            concatPointInPixels.push_back(IndP);
          }
        }
      }
    }
    trace.info()<<concatPointInPixels.size()<<std::endl;
  }



  //uncomment to create groundTruth image
  /*std::vector<int> groundtrueIds;
  IOHelper::readIntsFromFile(gtName, groundtrueIds);
  std::cout << gtName<< std::endl;
  cv::Mat gtimg=unrolled_map.makeGroundTruthImage(groundtrueIds);
  imwrite( "../../unrollSurfaceOutput/"+outputFileName+"_GT.png",unrolled_map.makeGroundTruthImage(groundtrueIds));*/

  //Uncomment to analyse relief image
  //ImageAnalyser image_analyser(unrolled_map,ind_Patches,myPoints,coefficients);
  //image_analyser.analyse();

  return concatPointInPixels;
}
