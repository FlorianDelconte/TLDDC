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
    Copy constructor
    **/
    UnrolledMap();
    /**
    Constructor.
   **/
    UnrolledMap(std::vector<CylindricalPoint> CylindricalPoints, std::vector<double> DeltaDistance);
    /**
    Copy constructor
    **/
    UnrolledMap(const UnrolledMap &um);
    
    /**
    return false if cells is consiedred out of the mesh (function to skip noise in relief image)
    **/
    bool detectCellsIn(unsigned int i, unsigned int j);
    
    /**
    compute normalized image (all values are in [0,1]) of unrolled map. dF is the decrease factor parameter : 1/dF.
     **/
    void computeNormalizedImage(int dF);
    /**
    compute normalized image of unrolled map with the black pixel complete by a multi scale analyse of unrolled map
     **/
    void computeNormalizedImageMultiScale();
    /**
    compute rgb image from normalized image
    **/
    void computeRGBImage();
    /**
    return normalized image
    **/
    cv::Mat getNormalizedImage();
    /**
    return rgb image 
    **/
    cv::Mat getRgbImage();
    /**
    return the vector of ind at pos (i,j) in unrolled_surface with the decrease factor specify by df : 1/dF.
    **/
    std::vector<unsigned int > getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF);
    /**
    return the point in cylindrical coordinate
    **/
    CylindricalPoint getCPoint(unsigned int);
  protected:
    //to compute normalization on delta dist
    double minRelief;
    double maxRelief;
    /**
     Normalize radius betwen [0,1] with min and max radius
     **/
    double normalizeReliefRepresentation(double value);
    /**
    return the mean of relief representation. You can specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double meansReliefRepresentation(unsigned int i, unsigned int j,int dF);
    /**
    return the max of relief representation. You can specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double maxReliefRepresentation(unsigned int i, unsigned int j,int dF);
    /**
    return the median of relief representation. You can specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double medianReliefRepresentation(unsigned int i, unsigned int j,int dF);
    //unrolled surface representation : each cells contain some index points.
    std::vector<std::vector<std::vector<unsigned int>>> unrolled_surface;
    //Represenation of the relief, radius of deltadiff.
    std::vector<double> reliefRepresentation;
    //Cylindricales Points
    std::vector<CylindricalPoint> CPoints;
    //discretisation
    int height_div, angle_div;
    //normalizedimage
    cv::Mat normalizedImage;
    //rgb image
    cv::Mat rgbImage;


};
#endif