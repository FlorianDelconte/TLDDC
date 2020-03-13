#ifndef UNROLL_MAP
#define UNROLL_MAP

#include <utility>
#include <iostream>
#include <vector>


#include <opencv2/opencv.hpp>

#include "CylindricalPoint.h"

class UnrolledMap{
  public:
    /**
    Constructor. TODO : prendre soit les points cylindrique soit la delta distance. 1 des 2. 
   **/
    UnrolledMap(std::vector<CylindricalPoint> CylindricalPoints, std::vector<double> DeltaDistance);
    /**
    return false if cells is consiedred as out of the mesh (function to skip noise in relief image)
    **/
    bool detectCellsIn(unsigned int i, unsigned int j);
    /**
    return the mean (on radius) of unrolled surface at position (i,j) . Specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double meansRadius(unsigned int i, unsigned int j,int dF);
    /**
    return the mean (dalta dist) of unrolled surface at position (i,j) . Specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double meansDeltaDist(unsigned int i, unsigned int j,int dF);
    /**
    return image (all values are in [0,1]) of unrolled map. dF is the decrease factor parameter : 1/dF.
    if df=1 return an image height_div*angle_div
     **/
    cv::Mat computeNormalizedImage(int dF);
    /**
    return a normalized image of unrolled map with the black pixel complete by a multi scale analyse of unrolled map
     **/
    cv::Mat computeNormalizedImageMultiScale();
  protected:
    //to compute normalization on radius
    double minRadius;
    double maxRadius;
    //to compute normalization on delta dist
    double minDeltaDist;
    double maxDeltaDist;
    /**
     Normalize radius betwen [0,1] with min and max radius
     **/
    double normalizeRadius(double value);
    /**
     Normalize radius betwen [0,1] with min and max radius
     **/
    double normalizeDeltaDist(double value);
    /**
    return the vector of ind at pos (i,j) in unrolled_surface with the decrease factor specify by df : 1/dF.
    **/
    std::vector<unsigned int > getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF);
    //unrolled surface representation : each cells contain some index points.
    std::vector<std::vector<std::vector<unsigned int>>> unrolled_surface;
    //Delta diff of all points
    std::vector<double> DDistance;
    //Cylindricales Points
    std::vector<CylindricalPoint> CPoints;
    //discretisation
    int height_div, angle_div;


};
#endif