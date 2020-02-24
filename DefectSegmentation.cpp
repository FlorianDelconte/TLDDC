#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>
#include <thread>
//debug
#include <stdlib.h>
#include <time.h>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Householder>

#include "DefectSegmentation.h"
#include "Regression.h"
#include "Statistic.h"
#include "IOHelper.h"
#include "MultiThreadHelper.h"

using namespace DGtal;


void
DefectSegmentation::init(){
    allocate();
    allocateExtra();

    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();

	convertToCcs();

    computeEquations();
    computeDistances();
    //debug
    //writeDebugInfo();
}

void
DefectSegmentation::allocateExtra(){
    //allocate
    coefficients.resize(pointCloud.size());
}


void
DefectSegmentation::computeEquations(){
    //w = patch width, wh = patch height
    //build kdtree using pcl
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
        ts.push_back(std::thread(&DefectSegmentation::computeEquationsMultiThread, this, i, nbCores, kdtree, minHeight, maxHeight));
    }
    computeEquationsMultiThread(nbCores - 1, nbCores, kdtree, minHeight, maxHeight);
    for(int i = 0; i < nbCores - 1; i++){
        ts[i].join();
    }
    trace.info()<<"finish eq"<<std::endl;
}
void
DefectSegmentation::computeEquationsMultiThread(int threadId, int nbThread,
        const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree,
        double minHeight, double maxHeight){

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
        //coefficients[i] = Regression::robustLinearOls(lengthForEstimate, radiusForEstimate);
        //coefficients2[i] = coefficients[i];// = Regression::robustLinearOls(lengthForEstimate, radiusForEstimate);
        coefficients[i] = Regression::linearRegression(lengthForEstimate, radiusForEstimate);
        //coefficients[i] = coefficients2[i];
    }
}


void DefectSegmentation::computeDistances(){

    for(unsigned int i = 0; i < myPoints.size(); i++){
        std::pair<double, double> coeffs = coefficients[i];
        if(coeffs.second == 0.0){
            distances[i] = 0;
        }else{
            double estimateRadii = myPoints[i].height * coeffs.first + coeffs.second;
            distances[i] = myPoints[i].radius - estimateRadii;
        }
    }
}

std::pair<double, double> DefectSegmentation::getCoeffs(unsigned int pointId){
    return coefficients.at(pointId);
}

std::pair<double, double> DefectSegmentation::getCoeffs2(unsigned int pointId){
    return coefficients.at(pointId);
}
