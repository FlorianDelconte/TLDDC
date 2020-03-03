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
    computeDistances();

    //WARNING :Need to placed after computing eq
    computeMeasures();
    
}


void
DefectSegmentationUnroll::allocateExtra(){
  //allocate
  //pre allocate size of unrolled surface
  //unrolled_surface.resize(height_div);
  //for (int i = 0; i < height_div; ++i)
  //  unrolled_surface[i].resize(angle_div);
  //pre allocate size of coefficent of line patch
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
  double minHeight = (*minMaxElem.first).height;
  double maxHeight = (*minMaxElem.second).height;

  int nbCores = getNumCores();
  std::vector<std::thread> ts;
  for(int i = 0; i < nbCores - 1; i++){
      ts.push_back(std::thread(&DefectSegmentationUnroll::computeEquationsMultiThread, this, i, nbCores, kdtree, minHeight, maxHeight));
  }
  computeEquationsMultiThread(nbCores - 1, nbCores, kdtree, minHeight, maxHeight);
  for(int i = 0; i < nbCores - 1; i++){
      ts[i].join();
  }
  trace.info()<<"finish eq"<<std::endl;
}

void
DefectSegmentationUnroll::computeEquationsMultiThread(int threadId, int nbThread,const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, double minHeight, double maxHeight){
  double patchAngle = arcLength / radii;
  for(unsigned int i = threadId; i < pointCloud.size();i+=nbThread){
    Z3i::RealPoint currentPoint = pointCloud.at(i);
    CylindricalPoint mpCurrent = myPoints.at(i);
    pcl::PointXYZ searchPoint(currentPoint[0], currentPoint[1], currentPoint[2]);
    std::vector<int> pointIdx;
    std::vector<float> pointRadiusSquaredDistance;

    std::vector<double> radiusForEstimate;
    std::vector<double> lengthForEstimate;

    double searchRadius = patchHeight / 2 + 1;
    if(mpCurrent.height - minHeight < patchHeight /2){
        searchRadius += mpCurrent.height - minHeight;
    }else if(maxHeight - mpCurrent.height < patchHeight/2){
        searchRadius += maxHeight - mpCurrent.height;
    }
    if ( kdtree.radiusSearch (searchPoint, searchRadius, pointIdx, pointRadiusSquaredDistance) > 0 ){
        for (unsigned int idx = 0; idx < pointIdx.size (); ++idx){
            //index of dgtal and pcl is the same
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
  coefficients[i] = Regression::linearRegression(lengthForEstimate, radiusForEstimate);
  }
}

void
DefectSegmentationUnroll::computeDistances(){
  for(unsigned int i = 0; i < myPoints.size(); i++){
      std::pair<double, double> coeffs = coefficients[i];
      if(coeffs.second == 0.0){
          distances[i] = 0;
      }else{
          double estimateRadii = myPoints[i].height * coeffs.first + coeffs.second;
          distances[i] = myPoints[i].radius - estimateRadii;
          //if(distances[i]<0){
          //   distances[i]=0.0;
          //}
          //trace.info()<<"estimation du rayon:  "<<estimateRadii<<std::endl;
          //trace.info()<<"vrais rayon: "<<myPoints[i].radius<<std::endl;
          //trace.info()<<"difference: "<<distances[i]<<std::endl;
      }
  }

}



std::vector<unsigned int >
DefectSegmentationUnroll::getPointfromSegmentId(unsigned int id){
  std::vector<unsigned int> output;
  for(unsigned int i = 0; i < myPoints.size(); i++){

      if(myPoints.at(i).segmentId==id){

        output.push_back(i);
      }
  }
  return output;
}
/*******************************************************************************************************************************/
/**************************************************UNROLLED MAP*****************************************************************/
/*******************************************************************************************************************************/
void
DefectSegmentationUnroll::computeMeasures(){
  trace.info()<<"Compute measures..."<<std::endl;
  //compute mean circumference (approx by circle) => cricumference is two times Pi times radius
  double radius_average=0.;
  double circum;
  double height;
  CylindricalPoint mpCurrent;

  for(unsigned int i = 0; i < myPoints.size(); i++){
    mpCurrent=myPoints.at(i);
    radius_average+=mpCurrent.radius;
  }
  radius_average/=myPoints.size();
  //compute angle measures
  CylindricalPointOrderAngle angleOrder;
  auto minMaxElem = std::minmax_element(myPoints.begin(), myPoints.end(), angleOrder);
  minAngle = (*minMaxElem.first).angle;
  maxAngle = (*minMaxElem.second).angle;
  //compute radius measures
  CylindricalPointOrderRadius RadiusOrder;
  minMaxElem = std::minmax_element(myPoints.begin(), myPoints.end(), RadiusOrder);
  minRadius = (*minMaxElem.first).radius;
  maxRadius = (*minMaxElem.second).radius;
  circum= 2*M_PI*radius_average;
  //compute height measures
  CylindricalPointOrder heightOrder;
  minMaxElem = std::minmax_element(myPoints.begin(), myPoints.end(), heightOrder);
  minHeight = (*minMaxElem.first).height;
  maxHeight = (*minMaxElem.second).height;
  height=maxHeight-minHeight;
  //compute difference distance measures
  
  maxdiffDist=INT_MIN;
  mindiffDist=INT_MAX;
  for(unsigned int i = 0; i < myPoints.size(); i++){
    if(distances[i]>maxdiffDist){
      maxdiffDist=distances[i];
    }
    if(distances[i]<mindiffDist){
      mindiffDist=distances[i];
    }
    //trace.info()<<distances[i]<<std::endl;
  }
  
  //Discretisation of height and circum
  height_div=roundf(height);
  angle_div=roundf(circum);
  trace.info()<<"Dicrestisation : X = "<<angle_div<<" Y = "<<height_div<<std::endl;
  trace.info()<<"ANGLE : [ "<<minAngle<<" ,  "<<maxAngle<<"]"<<std::endl;
  trace.info()<<"RADIUS : [ "<<minRadius<<" ,  "<<maxRadius<<"]"<<std::endl;
  trace.info()<<"HEIGHT : [ "<<minHeight<<" ,  "<<maxHeight<<"]"<<std::endl;
  trace.info()<<"DISTDIFF : [ "<<mindiffDist<<" ,  "<<maxdiffDist<<"]"<<std::endl;
}

void
DefectSegmentationUnroll::unrollSurface() {
  //preallocate size of unrolledmap
  unrolled_surface.resize(height_div);
  for (int i = 0; i < height_div; ++i)
    unrolled_surface[i].resize(angle_div);

  trace.info()<<"Compute Unrolled map..."<<std::endl;
  CylindricalPoint mpCurrent;
  int posAngle, posHeight;
  for(unsigned int i = 0; i < myPoints.size(); i++){
    mpCurrent=myPoints.at(i);
    //change range [minAngle,maxAngle] to [0,angle_div-1]
    posAngle=roundf((((angle_div-1)/(maxAngle-minAngle))*(mpCurrent.angle-(maxAngle)))+(angle_div-1));
    //change range [minHeight,maxHeight] to [0,height_div-1]
    posHeight=roundf((((height_div-1)/(maxHeight-minHeight))*(mpCurrent.height-maxHeight))+(height_div-1));
    //add index point to the unrolled_surface
    unrolled_surface[posHeight][posAngle].push_back(i);
  }
  /*Uncomment to test mesh reconstruction from unrolled surface
  Mesh<Z3i::RealPoint> meshTest;
  for(unsigned int i = 0; i < height_div; i++){
    for(unsigned int j = 0; j < angle_div; j++){
      for(std::vector<unsigned int>::iterator it = std::begin(two_d[i][j]); it != std::end(two_d[i][j]); ++it) {
          meshTest.addVertex(pointCloud.at(*it));
      }
    }
  }
  IOHelper::export2OFF(meshTest, "segmentPoints.off");*/
}
bool
DefectSegmentationUnroll::detectCellsIn(unsigned int i, unsigned int j){
  //A cell is 'in' if we can't reach the top image or the bot image with empty cells
  bool CellsIn=true;
  
  std::vector<unsigned int > tempUp = unrolled_surface[i][j];
  std::vector<unsigned int > tempDown = unrolled_surface[i][j];
  //ind for reach the top image
  unsigned int iUp=i;
  //ind for reach the bot image
  unsigned int iDown=i;
  //decrement ind while cells is empty
  while(tempUp.empty() && (iUp>0)){
    iUp-=1;
    tempUp=unrolled_surface[iUp][j];
  }
  //increment ind while cells is empty
  while(tempDown.empty() && (iDown<height_div-1)){
    iDown+=1;
    tempDown=unrolled_surface[iDown][j];
  }
  //check reached top of bot
  if((iUp==0 && tempUp.empty()) || (iDown==height_div-1 && tempDown.empty())){
    CellsIn=false;
  }
  return CellsIn;
}

/*******************************************************************************************************************************/
/**************************************************NORMALIZED IMAGE*************************************************************/
/*******************************************************************************************************************************/
void
DefectSegmentationUnroll::computeNormalizedImage() {
  trace.info()<<"start compute normalized image with max resolution : X = "<<angle_div<< " Y = "<<height_div <<std::endl;
  normalizedMap=cv::Mat::zeros(height_div,angle_div,CV_32FC4);
  double moyenneRadius;
  double moyenneDiffDist;
  double normalizedRadius;
  double normalizeDiffDist;
  std::vector<unsigned int > temp;
  for(unsigned int i = 0; i < height_div; i++){
    for(unsigned int j = 0; j < angle_div; j++){
      temp=unrolled_surface[i][j];
      //moyenneRadius=getMaxRadius(unrolled_surface[i][j]);
      moyenneRadius=getMeansRadius(unrolled_surface[i][j]);
      moyenneDiffDist=getMeansDistDiff(temp);
      //moyenneRadius=getMaxRadius(unrolled_surface[i][j]);
      //
      if(!temp.empty()){
        normalizedRadius=normalizeRadius(moyenneRadius);
        normalizeDiffDist=normalizeDiffDistance(moyenneDiffDist);
        normalizedMap.at<double>(i, j) = normalizeDiffDist;  
      }
    }
  }
}

void
DefectSegmentationUnroll::computeNormalizedImage(int decreaseFactor) {
  trace.info()<<"start compute normalized image with decrease factor : 1 / "<<decreaseFactor<<std::endl;
  //resolution of relief image
  unsigned int resX=angle_div/decreaseFactor;
  unsigned int resY=height_div/decreaseFactor;
  trace.info()<<"new resolution : X = "<<resX<<" Y = "<<resY<<std::endl;
  //init normalized map with the new resolution
  normalizedMap=cv::Mat::zeros(resY+1,resX+1,CV_32FC4);
  int topLeftCornerHeight;
  int topLeftCornerTheta;
  double radius=0.;
  double normalizedRadius;
  std::vector<unsigned int > temp;
  int XNormMap,YNormMap;
  //loop on the top left corner of all new cells
  for(unsigned int i = 0; i < height_div; i+=decreaseFactor){
    for(unsigned int j = 0; j < angle_div; j+=decreaseFactor){
      //check if the cells is in the mesh -> to avoid some noise in the image
      if(detectCellsIn(i,j)){
        //[0;height_div] to [0;resY]
        YNormMap =i/decreaseFactor;
        //[0;angle_div] to [0;resX]
        XNormMap=j/decreaseFactor;
        //get all indice Points of i,j in a lower resolution
        temp=getIndPointsInLowerResolution(i,j,decreaseFactor);
        //decreaseFactor give in parameter need to be choosen for never have an empty vector in the lower resolution
        assert(!temp.empty());
        radius=getMaxRadius(temp);
        normalizedRadius=normalizeRadius(radius);
        normalizedMap.at<double>(YNormMap, XNormMap) = normalizedRadius;
      }
    }
  }
  
}

void
DefectSegmentationUnroll::computeNormalizedImageMultiScale() {
  trace.info()<<"start compute Multi-Scale normalized image..."<<std::endl;
  int decreaseHit;
  int maxDecreaseHit=0;
  int current_resX,current_resY;
  double moyenneRadius;
  double moyenneDiffDist;
  int topLeftCornerHeight,topLeftCornerTheta;

  normalizedMap=cv::Mat::zeros(height_div,angle_div,CV_32FC4);
  double normalizedRadius;
  double normalizedDistDiff;
  std::vector<unsigned int > temp;
  
  //loop on all cells of unrolled surface
  for(unsigned int i = 0; i < height_div; i++){
    for(unsigned int j = 0; j < angle_div; j++){
      if(detectCellsIn(i,j)){
        decreaseHit=2;
        moyenneRadius=0.;
        current_resX=angle_div;
        current_resY=height_div;
        //get the vector of point in cells i,j
        temp=unrolled_surface[i][j];
        //while this ector is empty find a little resolution where the corresponding i,j cells is not empty
        while(temp.empty()&& decreaseHit<32){
          //get all indice Points of i,j in a lower resolution
          temp=getIndPointsInLowerResolution(i,j,decreaseHit);
          //decrease resolution (1/decreaseHit)
          decreaseHit*=2;
          //keep the max decrease resolution
          if(decreaseHit>=maxDecreaseHit){
            maxDecreaseHit=decreaseHit;
          }
        }

        //if t is empty, thats means that even with a multi scale resolution( to 1/(2^5)) we can't find info
        if(!temp.empty()){
          //trace.info()<<temp.size()<<std::endl;
          //moyenneRadius=getMeansRadius(temp);
          moyenneDiffDist=getMeansDistDiff(temp);
          
          //trace.info()<<moyenneRadius<<std::endl;
          //normalizeDiffDistance(moyenneRadius);
          //moyenneRadius=getMaxRadius(t);
          //moyenneRadius=getMedianRadius(t);
          
          //normalizedRadius=normalizeDiffDistance(moyenneRadius);
          //normalizedMap.at<double>(i, j) = normalizedRadius;
        }
        //normalizedRadius=normalizeRadius(moyenneRadius);
        normalizedDistDiff=normalizeDiffDistance(moyenneDiffDist);
        //if(normalizedDistDiff>normalizedRadius){
        //  normalizedRadius=normalizedDistDiff;
        //}
        //trace.info()<<normalizedRadius<<std::endl;
        normalizedMap.at<double>(i, j) = normalizedDistDiff;
        //normalizedMap.at<double>(i, j) = normalizedRadius;

      }
    }
  }
  //computeNormalizedImage(maxDecreaseHit);
}
double 
DefectSegmentationUnroll::normalizeRadius(double value){
  //search min max of distance
  return ((1/(maxRadius-minRadius))*(value-maxRadius))+1;
}

double 
DefectSegmentationUnroll::normalizeDiffDistance(double value){
  return ((1/(maxdiffDist-mindiffDist))*(value-maxdiffDist))+1;
}

std::vector<unsigned int > 
DefectSegmentationUnroll::getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF){
  std::vector<unsigned int > outPutInd;
  int topLeftCornerHeight,topLeftCornerTheta;
  topLeftCornerHeight=(i/dF)*dF;
  topLeftCornerTheta=(j/dF)*dF;
  //loop on cells (of unrolledSurace) in region containing (i,j)
  for( int k = topLeftCornerHeight; k < (topLeftCornerHeight+dF); k++){
    for( int l = topLeftCornerTheta; l < (topLeftCornerTheta+dF); l++){
      //check if cells of region is in unrolled surface -> check no segmentation fault
      if((k < height_div)&&(l < angle_div)){
      //unlarge  vector size
      outPutInd.reserve(outPutInd.size() + unrolled_surface[k][l].size());
      //concat vector
      outPutInd.insert(outPutInd.end(), unrolled_surface[k][l].begin(),  unrolled_surface[k][l].end());
      }
    }
  }
  return outPutInd;
}
double
DefectSegmentationUnroll::getMeansDistDiff(std::vector<unsigned int > v) {
  double moyenneRadius=0.;
  if(v.size()>0){
    unsigned int IndP;
    for(std::vector<unsigned int>::iterator it = std::begin(v); it != std::end(v); ++it) {
      IndP=*it;
      moyenneRadius+=distances[IndP];
      //trace.info()<<moyenneRadius<<std::endl;
    }moyenneRadius
    /=v.size();
  }
  //if(moyenneRadius>10){
  //  trace.info()<<moyenneRadius<<std::endl;
  //}
  
  return moyenneRadius;
}
double
DefectSegmentationUnroll::getMeansRadius(std::vector<unsigned int > v) {
  double moyenneRadius=0.;
  if(!v.empty()){
    unsigned int IndP;
    for(std::vector<unsigned int>::iterator it = std::begin(v); it != std::end(v); ++it) {
      IndP=*it;
      CylindricalPoint mpCurrent = myPoints.at(IndP);
      moyenneRadius+=mpCurrent.radius;
    }
    moyenneRadius/=v.size();
  }
  //trace.info()<<moyenneRadius<<std::endl;
  return moyenneRadius;
}

double
DefectSegmentationUnroll::getMaxRadius(std::vector<unsigned int > v) {
  double maxRadius=INT_MIN;
  for(std::vector<unsigned int>::iterator it = std::begin(v); it != std::end(v); ++it) {
      CylindricalPoint mpCurrent = myPoints.at(*it);
      if(mpCurrent.radius >maxRadius){
        maxRadius=mpCurrent.radius;
      }
  }
  return maxRadius;
}
double
DefectSegmentationUnroll::getMedianRadius(std::vector<unsigned int > v) {
  std::vector<double> vectorOfRadius;
  unsigned int size = v.size();
  double median=0.0;
  if (size == 0)
  {
    median=0.;  // Undefined, really.
  }else{
    for(auto it = v.begin(); it != v.end(); it++) {
      double radiusCurrent = myPoints.at(*it).radius;
      vectorOfRadius.push_back(radiusCurrent);
    }

    std::sort(vectorOfRadius.begin(), vectorOfRadius.end());
    if (size % 2 == 0){
      median=(vectorOfRadius.at(size / 2 - 1) + vectorOfRadius.at(size / 2)) / 2;
    }else{
      median=vectorOfRadius.at(size / 2);
    }
  }
  return median;
}

/*******************************************************************************************************************************/
/**************************************************RGB IMAGE********************************************************************/
/*******************************************************************************************************************************/
void
DefectSegmentationUnroll::createVisuImage(std::string s) {
  
  int grayscaleValue;
  double normalizedValue;
  int rows = normalizedMap.rows;
  int cols = normalizedMap.cols;
  cv::Mat grayscalemap(rows,cols,CV_8UC1,cv::Scalar(0));
  cv::Mat reliefPictures(rows, cols, CV_8UC3, cv::Scalar(110, 110, 110));

  for(unsigned int i = 0; i < rows; i++){
    for(unsigned int j = 0; j < cols; j++){
      normalizedValue=normalizedMap.at<double>(i, j);
      grayscaleValue=((255/1)*(normalizedValue-1))+255;
      grayscalemap.at<uchar>(i, j) = grayscaleValue;
    }
  }
  cv::applyColorMap(grayscalemap, reliefPictures, cv::COLORMAP_JET);
  imwrite( "unrollSurfacePictures/"+s, reliefPictures);
}