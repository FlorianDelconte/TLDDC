#ifndef IMAGE_ANALYSER
#define IMAGE_ANALYSER

#include <utility> 
#include <opencv2/opencv.hpp>

#include "UnrolledMap.h"
#include "CylindricalPoint.h"

class ImageAnalyser{
  public:
    /**
    Constructor.
   **/
    ImageAnalyser(UnrolledMap uM, std::vector<std::vector<unsigned int>> pR,std::vector<CylindricalPoint> cP, std::vector<std::pair<double, double> > coefs):unrolled_map(uM),ind_Patches(pR),CPoints(cP),ind_CoefsLines(coefs)
    {}
    /**
    display rgb image from unrolled map and allow user event
    **/
    void analyse();

    CylindricalPoint getShiftedPoint(std::pair<double, double>& coef,std::vector<unsigned int> pointsInPatch);
    /**
    menber function to setMouseCallBack
    **/
    static void onMouse( int event, int x, int y, int flags, void* param);

  protected:
    


    std::vector<std::vector<unsigned int>> ind_Patches;

    UnrolledMap unrolled_map;
    
    std::vector<CylindricalPoint> CPoints;

    std::vector<std::pair<double, double> > ind_CoefsLines;
};

#endif