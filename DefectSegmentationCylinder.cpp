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

#include "Regression.h"
#include "Statistic.h"
#include "IOHelper.h"
#include "MultiThreadHelper.h"

#include "DefectSegmentationCylinder.h"

using namespace DGtal;


void
DefectSegmentationCylinder::init(){
    allocate();
    allocateExtra();
//    computeSegmentsAndRadii();

    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();
//    computeAngleOfPoints();
//    computeCells();
    computeEquations();
    computeDistances();
    //debug
    //writeDebugInfo();
}

void 
DefectSegmentationCylinder::allocateExtra(){
    //allocate
    coefficients.resize(nbSegment);
    extremities.resize(nbSegment);
}

void DefectSegmentationCylinder::computeEquations(){
    std::vector< std::vector<unsigned int> > cyls(fiber.size() - 1);
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint p = pointCloud.at(i);
        int segId = getSegment(p);
        cyls[segId].push_back(i);
    }
    for(unsigned int i = 0; i < cyls.size(); i++){
        std::vector<double> radiusOfSegment;
        std::vector<unsigned int> pointsOnSegment = cyls[i];
        for(unsigned int j = 0; j < pointsOnSegment.size(); j++){
            radiusOfSegment.push_back(getRadius(pointsOnSegment[j]));
        }
    }
    for(unsigned int i = 0; i < cyls.size(); i++){
        std::vector<unsigned int> pointsOnSegment = cyls[i];
        pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);
        pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
        for(unsigned int j = 0; j < pointsOnSegment.size(); j++){
            unsigned int pIndex = pointsOnSegment[j];
            Z3i::RealPoint p = pointCloud.at(pIndex);
            PointT pt(p[0], p[1], p[2]);
            cloud->points.push_back(pt);
        } pcl::NormalEstimation<PointT, pcl::Normal> ne;
        pcl::SACSegmentationFromNormals<PointT, pcl::Normal> seg; 

        pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT> ());
        pcl::PointCloud<pcl::Normal>::Ptr cloudNormals (new pcl::PointCloud<pcl::Normal>);
        pcl::ModelCoefficients::Ptr coeffs(new pcl::ModelCoefficients);
        pcl::PointIndices::Ptr inliers(new pcl::PointIndices);

        ne.setSearchMethod (tree);
        ne.setInputCloud (cloud);
        ne.setKSearch (30);
        ne.compute (*cloudNormals);

        seg.setOptimizeCoefficients (true);
        seg.setModelType (pcl::SACMODEL_CYLINDER);
        seg.setMethodType (pcl::SAC_RANSAC);
        seg.setNormalDistanceWeight (1);
        seg.setMaxIterations (10000);
        seg.setDistanceThreshold (1);
        seg.setRadiusLimits (20, 120);
        seg.setInputCloud (cloud);
        seg.setInputNormals (cloudNormals);
        seg.segment (*inliers, *coeffs);

        double x1 = coeffs->values[0];
        double y1 = coeffs->values[1];
        double z1 = coeffs->values[2];
        double x2 = x1 + coeffs->values[3];
        double y2 = y1 + coeffs->values[4];
        double z2 = z1 + coeffs->values[5];
        double ra = coeffs->values[6];

        Z3i::RealPoint aPointOnCenter(x1, y1, z1);
        Z3i::RealPoint aDirection(coeffs->values[3], coeffs->values[4], coeffs->values[5]);
        CylindricalCoefficients cylCoeffs;
        cylCoeffs.point = aPointOnCenter;
        cylCoeffs.direction = aDirection;
        cylCoeffs.radius = ra;

        
        coefficients[i] = cylCoeffs;
        Z3i::RealPoint center1(x1,y1,z1);
        Z3i::RealPoint center2(x2,y2,z2);
        std::pair<Z3i::RealPoint, Z3i::RealPoint> exs = getExtremityOfCylinder(pointCloud, pointsOnSegment, center1, center2);
        extremities[i] = exs;
        //std::vector<Z3i::RealPoint> centerLine;
        //centerLine.push_back(extremities.first);
        //centerLine.push_back(extremities.second);
        //Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, centerLine, ra, 0.1, DGtal::Color::Red);

        
        for(unsigned int j = 0; j < pointsOnSegment.size(); j++){
            unsigned int pIndex = pointsOnSegment[j];
            Z3i::RealPoint p = pointCloud.at(pIndex);
            Z3i::RealPoint aVect = aPointOnCenter - p;
            double radii = (aDirection.crossProduct(aVect)).norm()/aDirection.norm();
            
//double d = std::abs( (y2 - y1)*x0 - (x2 - x1)*y0 + x2*y1 -y2*x1 )/sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2-y1));
//std::cout<<radii <<"#"<<ra<<std::endl;
            distances[pIndex] = radii - ra;
        }
    }
}

void DefectSegmentationCylinder::computeDistances(){
}

std::vector<CylindricalCoefficients> DefectSegmentationCylinder::getCoeffs(){
    return coefficients;
}

std::vector<std::pair<Z3i::RealPoint, Z3i::RealPoint> > 
DefectSegmentationCylinder::getExtremityOfCylinders(){
    return extremities;
}

std::pair<Z3i::RealPoint, Z3i::RealPoint> 
DefectSegmentationCylinder::getExtremityOfCylinder(const std::vector<Z3i::RealPoint> &pointCloud, 
                       const std::vector<unsigned int> &pointsOnCylinder, 
                       Z3i::RealPoint center1, Z3i::RealPoint center2){
    Z3i::RealPoint vectDir = (center2 - center1).getNormalized();
    double minScalar = std::numeric_limits<double>::max();
    double maxScalar = -std::numeric_limits<double>::max();
    Z3i::RealPoint minPoint(0,0,0);
    Z3i::RealPoint maxPoint(0,0,0);
    for(unsigned int i = 0; i < pointsOnCylinder.size(); i++){
        unsigned int pIndex = pointsOnCylinder.at(i);
        Z3i::RealPoint p = pointCloud.at(pIndex);

        Z3i::RealPoint aVect = p - center1;
        double scalar = aVect.dot(vectDir);
        if(scalar < minScalar ){
            minScalar = scalar;
            minPoint = center1 + scalar*vectDir;
        }
        if(scalar > maxScalar){
            maxScalar = scalar;
            maxPoint = center1 + scalar*vectDir;
        }
    }
    std::pair<Z3i::RealPoint, Z3i::RealPoint> extre(minPoint, maxPoint);
    return  extre;
}
