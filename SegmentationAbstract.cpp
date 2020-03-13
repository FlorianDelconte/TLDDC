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

#include "SegmentationAbstract.h"
#include "Regression.h"
#include "Statistic.h"
#include "IOHelper.h"
#include "MultiThreadHelper.h"


using namespace DGtal;


bool
CylindricalPointOrder::operator() (CylindricalPoint p1, CylindricalPoint p2) {
        return p1.height < p2.height;
}

bool
CylindricalPointOrderRadius::operator() (CylindricalPoint p1, CylindricalPoint p2) {
        return p1.radius < p2.radius;
}
bool
CylindricalPointOrderAngle::operator() (CylindricalPoint p1, CylindricalPoint p2) {
        return p1.angle < p2.angle;
}

bool
CylindricalPointOrderArcLenght::operator() (CylindricalPoint p1, CylindricalPoint p2) {
        double arcLengtP1=p1.radius*p1.angle;
        double arcLengtP2=p2.radius*p2.angle;
        return arcLengtP1 < arcLengtP2;
}

void
SegmentationAbstract::init(){
    allocate();
//    computeSegmentsAndRadii();

    computeBeginOfSegment();
    computeVectorMarks();
    computePlaneNormals();
	convertToCcs();
//    computeAngleOfPoints();
//    computeCells();
    computeEquations();
    computeDistances();
    //debug
    //writeDebugInfo();
}

void
SegmentationAbstract::allocate(){
    nbSegment = fiber.size() - 1;
    //allocate
    myPoints.resize(pointCloud.size());
    distances.resize(pointCloud.size());
   //coefficients.resize(pointCloud.size());
    //debug
    //coefficients2.resize(pointCloud.size());

    beginOfSegment.resize(nbSegment);
    segmentLength.resize(nbSegment);
    vectMarks.resize(nbSegment);
    ns.resize(nbSegment);
}


double SegmentationAbstract::findThresholdRosin(){
    //build histogram
    double res = binWidth; //must be configurable?
//trace.error()<<"binWidth:"<<binWidth<<std::endl;
    double maxValue = *std::max_element(distances.begin(), distances.end());
    double minValue = *std::min_element(distances.begin(), distances.end());
    double range = maxValue - minValue;
    int nbInterval = range / res;

    std::vector<int> histogram(nbInterval, 0);
    for(unsigned int i = 0; i < myPoints.size(); i++){
        int index = (distances.at(i) - minValue)/res;
        histogram[index]++;
    }
    std::vector<int>::iterator maxFreq = std::max_element(histogram.begin(), histogram.end());

    int maxFreqValue = *maxFreq;
    int maxFreqIndex = std::distance(histogram.begin(), maxFreq);
    assert(maxFreqValue == histogram.at(maxFreqIndex));

    unsigned int lastIndex = histogram.size() - 1;
    int lastValue = histogram.at(lastIndex);


    for(unsigned int i = maxFreqIndex; i < histogram.size(); i++){
        if(histogram.at(i) == 0){
            lastIndex = i;
            lastValue = 0;
            break;
        }
    }

    double valueDiff = lastValue - maxFreqValue;
    double valueDiff2 = valueDiff *valueDiff;
    double indexDiff = lastIndex - maxFreqIndex;
    double indexDiff2 = indexDiff * indexDiff;
    double bestThresIndex = maxFreqIndex;
    double bestDist = 0;

    //line between maxFreq and last element of historgram
    double a = (lastValue - maxFreqValue)*1.0/(lastIndex - maxFreqIndex);
    double b = maxFreqValue - a * maxFreqIndex;

    //for (std::vector<int>::iterator it = maxFreq ; it != histogram.end(); ++it){
    for (unsigned int i = maxFreqIndex; i < lastIndex; i++){
        //trace.info()<<valueDiff *  i - indexDiff*histogram.at(i) + maxFreqIndex*lastValue - maxFreqValue * lastIndex<<std::endl;
        double dist = std::abs(valueDiff *  i - indexDiff*histogram.at(i) + maxFreqValue*lastIndex - maxFreqIndex * lastValue)/
            sqrt(valueDiff2 + indexDiff2 );
        if(dist > bestDist){
            bestDist = dist;
            bestThresIndex = i;
            //trace.info()<<"dist:"<<bestDist<< std::endl;
            //trace.info()<<"bestThresIndex:"<<bestThresIndex<< std::endl;
        }
    }

    /**
     * perpendicular line: -1/a + c passe through bestThresIndex
     * -1/ax1 + c = y1
     *  c = y1 + 1/ax1
     *  ax2 + b = y2
     * -1/ax2 + c = y2
     *  -1/ax2 + y1 + 1/ax1 = y2
     *  (a + 1/a)x2 = y1 + 1/ax1 -b
     */
    int bestVal = histogram.at(bestThresIndex);
    double x2 = (bestVal + bestThresIndex/a - b)/(a + 1/a);
    double y2 = b + a*x2;

    std::vector<std::pair<double, double>> forPlot;
    std::pair<double, double> maxPoint(maxFreqIndex * res + minValue, maxFreqValue);
    std::pair<double, double> lastPoint(lastIndex * res + minValue, lastValue);
    std::pair<double, double> bestPoint(bestThresIndex * res + minValue, bestVal);
    std::pair<double, double> projBestPoint(x2* res + minValue, y2);

    forPlot.push_back(maxPoint);
    forPlot.push_back(lastPoint);
    forPlot.push_back(bestPoint);
    forPlot.push_back(projBestPoint);
    IOHelper::export2Text(forPlot, "pointFile");

    //histogram
    std::vector<std::pair<double, double>> histForPlot;
    for(unsigned int i = 0; i< histogram.size(); i++){
        std::pair<double, double> aBin(i* res + minValue, histogram.at(i));
        histForPlot.push_back(aBin);
    }

    IOHelper::export2Text(histForPlot, "hist2d");
	trace.info()<<"threshold: "<< bestThresIndex*res + minValue<<std::endl;
    return bestThresIndex*res + minValue;
}

std::vector<unsigned int>
SegmentationAbstract::getDefect(double threshold){
    std::vector<unsigned int> defects;
    for(unsigned int i = 0; i < myPoints.size(); i++){
        if(distances[i] > threshold){
            defects.push_back(i);
        }
    }
    return defects;
}

std::vector<unsigned int>
SegmentationAbstract::getDefect(){
    double th = findThresholdRosin();
    return getDefect(th);
}

Z3i::RealPoint
SegmentationAbstract::getDirectionVector(const unsigned int &segmentId){
    Z3i::RealPoint vectDir = fiber[segmentId + 1] - fiber[segmentId];
    return vectDir/vectDir.norm();
}


unsigned int
SegmentationAbstract::getSegment(const Z3i::RealPoint &aPoint){
    double lastSign = 1;
    for (int i = 1; i < nbSegment; i++){
        Z3i::RealPoint aVect = aPoint - fiber.at(i);
        double sign = aVect.dot(ns[i]);
        if(sign*lastSign <= 0){
            return i - 1;
        }
        lastSign = sign;
    }
    return nbSegment - 1;
}

//aDirection is normalized!
Z3i::RealPoint
SegmentationAbstract::getRadialVector(const Z3i::RealPoint &aPoint, const Z3i::RealPoint &aDirection, const Z3i::RealPoint &p0){
    double dist = aDirection.dot(aPoint - p0);
    Z3i::RealPoint proj = p0 + dist*aDirection;
    return aPoint - proj;
}



void SegmentationAbstract::convertToCcs(){
    double sumRadii = 0.0;
    for(unsigned int i = 0; i < pointCloud.size(); i++){
        Z3i::RealPoint aPoint = pointCloud.at(i);
        unsigned int segmentId = getSegment(aPoint);

        myPoints[i].segmentId = segmentId;
        assert(segmentId < fiber.size() - 1);

        Z3i::RealPoint aDir = fiber[segmentId + 1] - fiber[segmentId];
        aDir = aDir / aDir.norm();
        Z3i::RealPoint vectRadial = getRadialVector(pointCloud[i], aDir, fiber[segmentId]);
        double ra = vectRadial.norm();
        sumRadii += ra;
        //radius of point
        myPoints[i].radius = ra;

        Z3i::RealPoint p0 = fiber[segmentId];
        Z3i::RealPoint vectDir = getDirectionVector(segmentId);
        Z3i::RealPoint vect = aPoint - p0;
        double dist = vect.dot(vectDir);
        //z
        myPoints[i].height = beginOfSegment[segmentId] + dist;


        double angle = acos(vectRadial.dot(vectMarks[segmentId])/
                vectMarks[segmentId].norm()/vectRadial.norm());
        //Z3i::RealPoint crossProduct = vectMarks[segmentId].crossProduct(vectRadial);
        Z3i::RealPoint u = vectMarks[segmentId].crossProduct(vectDir);
        if (u.dot(vectRadial) < 0)
        {
            angle = 2 * M_PI - angle;
        }
        //angle
        myPoints[i].angle = angle;

    }
    radii = sumRadii / pointCloud.size();
}


void
SegmentationAbstract::computeBeginOfSegment(){
    beginOfSegment[0] = 0.0;
    for (int i = 1; i < nbSegment; i++){
        Z3i::RealPoint vectDir = fiber.at(i) - fiber.at(i - 1);
        beginOfSegment[i] = beginOfSegment[i - 1] + vectDir.norm();
    }
}

void
SegmentationAbstract::computePlaneNormals(){
    for (unsigned int i = 0; i < fiber.size() - 1; i++){
        Z3i::RealPoint vectDir = fiber.at(i + 1) - fiber.at(i);
        if(i == 0){
            ns[i] = vectDir.getNormalized();
        }else{
            Z3i::RealPoint previousVectDir = fiber.at(i) - fiber.at(i - 1);
            ns[i] = (previousVectDir + vectDir).getNormalized();
        }
    }
}

void
SegmentationAbstract::computeVectorMarks(){
    Z3i::RealPoint lastVectMark(0, 1 , 0);
    //dummy
    for (unsigned int i = 0; i < fiber.size() - 1; i++){
        //direction of current line
        //vector w
        Z3i::RealPoint vectDir = fiber.at(i + 1) - fiber.at(i);
        //normal vector of centerline and lastVectMark
        Z3i::RealPoint vectNormal = vectDir.crossProduct(lastVectMark);
        //first point of current line
        Z3i::RealPoint p0 = fiber.at(i);
        //The plane (vectDir, lastVectMark)
        double d = -(p0(0)*vectNormal(0) + p0(1)*vectNormal(1) + p0(2)*vectNormal(2));

        Z3i::RealPoint somePoint;
        somePoint[1] = p0(1) + 1;
        Eigen::Matrix3d A;
        Eigen::Vector3d b;
        A << 0,1,0,  vectDir[0],vectDir[1],vectDir[2],  vectNormal(0),vectNormal(1),vectNormal(2);
        b << p0(1) + 100, p0(0)*vectDir[0] + p0(1)*vectDir[1] + p0(2)*vectDir[2], -d;
        Eigen::Vector3d aPoint = A.colPivHouseholderQr().solve(b);
        Z3i::RealPoint aP(aPoint(0), aPoint(1), aPoint(2));
        Z3i::RealPoint w = aP - p0;

        vectMarks[i] = w.getNormalized();

        lastVectMark = vectMarks[i];
    }
}


std::vector<double>
SegmentationAbstract::getDistances(){
    return distances;
}


int SegmentationAbstract::getNbSegment(){
    return nbSegment;
}

int SegmentationAbstract::getNbSector(){
    return nbSector;
}

double SegmentationAbstract::getRadius(unsigned int index){
    return myPoints.at(index).radius;
}
double SegmentationAbstract::getLength(unsigned int index){
    return myPoints.at(index).height;
}

CylindricalPoint
SegmentationAbstract::getPointInCylindric(unsigned int pId){
    return myPoints.at(pId);
}
//-----------------------------------------------------
//demo
//-----------------------------------------------------
std::vector<unsigned int> SegmentationAbstract::getPatch(unsigned int pointIndex){
    std::vector<unsigned int> pIds;
    double patchAngle = arcLength / radii;
    //w = patch width, wh = patch height
    //build kdtree using pcl
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




    //trace.info()<<minHeight<<std::endl;
    Z3i::RealPoint currentPoint = pointCloud.at(pointIndex);
    CylindricalPoint mpCurrent = myPoints.at(pointIndex);
    pcl::PointXYZ searchPoint(currentPoint[0], currentPoint[1], currentPoint[2]);
    std::vector<int> pointIdx;
    std::vector<float> pointRadiusSquaredDistance;

    //double searchRadius = maxHeight - mpCurrent.height < patchHeight / 2 ? patchHeight/2 + maxHeight - mpCurrent.height : patchHeight/2;
    //double rightLimit = mpCurrent.height - minHeight < patchHeight / 2 ? patchHeight + mpCurrent.height - minHeight : patchHeight/2;
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
            /*
               if((mpFound.height <= mpCurrent.height && mpCurrent.height - mpFound.height > leftLimit) ||
               (mpFound.height > mpCurrent.height  && mpFound.height - mpCurrent.height > rightLimit)){
               continue;
               }
               */
            pIds.push_back(foundedIndex);
        }
    }
    return pIds;
}
