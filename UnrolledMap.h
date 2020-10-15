#ifndef UNROLL_MAP
#define UNROLL_MAP

#include <utility>
#include <iostream>
#include <vector>

#include "CylindricalPoint.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/boards/Board2D.h"

using namespace DGtal;
  //Image of double to store normalized values
  typedef ImageContainerBySTLVector<Z2i::Domain, float> Image2dNormalized;
  //Image of Char to make a grayscale image
  typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2dGrayScale;
  //Image of Color to make a rgb image
  typedef ImageContainerBySTLVector<Z2i::Domain, Color > imageRGB;

class UnrolledMap{
  public:
    /**
    Copy constructor
    **/
    UnrolledMap();
    /**
    Constructor.
   **/
    UnrolledMap(std::vector<CylindricalPoint> CylindricalPoints, std::vector<double> DeltaDistance,int dF,int gs_ori,int gs_intensity):
      maxDecreaseFactor(dF),
      CPoints(CylindricalPoints),
      reliefRepresentation(DeltaDistance),
      reliefImage(Z2i::Domain()),
      reliefImageGrayScale(Z2i::Domain()),
      reliefImageRGB(Z2i::Domain()),
      maxIndTop(0),
      zero_level_intensity(gs_ori),
      intensity_per_cm(gs_intensity),
      minIndBot(0){};

    /**
    Copy constructor
    **/
    UnrolledMap(const UnrolledMap &um):
      maxDecreaseFactor(um.maxDecreaseFactor),
      CPoints(um.CPoints),
      reliefRepresentation(um.reliefRepresentation),
      reliefImage(um.reliefImage),
      reliefImageGrayScale(um.reliefImageGrayScale),
      reliefImageRGB(um.reliefImageRGB),
      maxIndTop(um.maxIndTop),
      minIndBot(um.minIndBot),
      height_div(um.height_div),
      angle_div(um.angle_div){};

    /**
    return false if cells is consiedred out of the mesh (function to skip noise in relief image)
    **/
    bool detectCellsIn(unsigned int i, unsigned int j);
    /**
    Compute discretisation
     **/
    void computeDicretisation();
    /**
    compute normalized image (all values are in [0,1]) of unrolled map. dF is the decrease factor parameter : 1/dF.
     **/
    void computeNormalizedImage(int dF);
    /**
    compute normalized image of unrolled map with the black pixel complete by a multi scale analyse of unrolled map
     **/
    void computeNormalizedImageMultiScale();
    /**
    write rgb image from normalized image
    **/
    void computeRGBImage();
    /**
    compute rgb image from normalized image
    **/
    void computeGRAYImage();

    /**
    return the vector of ind at pos (i,j) in unrolled_surface with the decrease factor specify by df : 1/dF.
    **/
    std::vector<unsigned int > getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF);
    /**
    return the vector of ind at pos (i,j) in unrolled_surface. Need unrolled_sruface to be build
    **/
    std::vector<unsigned int > getPointsUnrolled_surface(unsigned int i,unsigned int j);
    /**
    return ground truth image. NOT ANUMORE DISPONIBLE. TODO : IMPLEMENT IN DGTal
    **/
    //void makeGroundTruthImage(std::vector<int> gtId,std::string outputFileName);
    /**
    * Crop bot and top of normalizd image
    **/
    void cropTopBotImage();
    /*********
    **GETTERS**
    **********/
    /**
    return rgb image
    **/
    Image2dGrayScale getReliefImageGrayScale();
    /**
    return rgb image
    **/
    imageRGB getReliefImageRGB();
    /**
    return normalized image
    **/
    Image2dNormalized getNormalizedImage();
    /**
    return the point in cylindrical coordinate
    **/
    CylindricalPoint getCPoint(unsigned int);
    /**
    return discretisation vectors
    **/
    std::vector<std::vector<std::vector<unsigned int>>> getDiscretisation();
    /**
    return the nb row cropped
    **/
    int getRowCroppedBot();
    /**
    return the nb row cropped
    **/
    int getRowCroppedTop();
  protected:

    /**
     transforme reliefImage betwen [0,255] with min and max reliefREP
     **/
    Image2dGrayScale toGrayscaleImageMinMax();
    /**
     Normalize reliefImage betwen [0,255] with a fixed start and padding
     **/
    Image2dGrayScale toGrayscaleImageFixed(int intensityPerCm,double reliefValueforZero);

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
    //double maxReliefRepresentationv2(unsigned int i, unsigned int j,int dF);
    /**
    return the sum of relief representation. You can specify the resolution by dF : 1/dF.
    if df=1 return the mean at position (i,j) in unrolled surface.
    **/
    double sumReliefRepresentation(unsigned int i, unsigned int j,int dF);
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
    //the maximum decrease factor for multi resolution research (2^n) with n = decreaseFactor
    int maxDecreaseFactor;
    //index of lines to be cropped
    unsigned int maxIndTop,minIndBot;
    //normalizedimage with DGtal
    Image2dNormalized reliefImage;
    //grayscale image with DGTAL
    Image2dGrayScale reliefImageGrayScale;
    //rgb image with DGTAL
    imageRGB reliefImageRGB;
    //relief value for the 0 level intensity in grayscale map
    int zero_level_intensity;
    //number of intensity value to represent 1 cm of relief
    int intensity_per_cm;
};
#endif
