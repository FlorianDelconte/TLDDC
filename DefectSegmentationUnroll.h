#ifndef UNROLL_SURFACE
#define UNROLL_SURFACE

#include <utility>


#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>

#include "SegmentationAbstract.h"
#include "CylindricalPoint.h"

#include <opencv2/opencv.hpp>

using namespace DGtal;

class DefectSegmentationUnroll : public SegmentationAbstract {
  public:

    using SegmentationAbstract::SegmentationAbstract;

    void init() override;
    /**
    unroll the surface --> for the moment unroll onely one segment
    **/
    void unrollSurface();
    /**
    compute normalized image from unrolled surface
    **/
    void computeNormalizedImage();
    /**
    create visual multi-scaled image from unrolled surface
    **/
    void computeNormalizedImageMultiScale();
    /**
    create visual image from normalized image
    **/
    void createVisuImage(std::string s);
  protected:
     /**
    return false if cells is consiedred as out of the mesh (function to skip noise in relief image)
    **/
    bool detectCellsIn(unsigned int i, unsigned int j);
    /**
    Cylindrical point from segment of centerline
    **/
    std::vector<unsigned int >  getPointfromSegmentId(unsigned int id);
    /**
    return the mean (on radius) of vector
    **/
    double getMeansVector(std::vector<unsigned int > v);
    /**
    return the max (on radius) of vector
    **/
    double getMaxVector(std::vector<unsigned int > v);
    /**
    return the median (on radius) of vector
    **/
    double getMedianVector(std::vector<unsigned int > v);
    /**
    Method to find the discretisation i / j of angle and height
    **/
    void computeDiscretisation();

    /**
    transfrome a 2Darray of vector into an image
    **/
    void transformToImage( std::vector<unsigned int > **two_d);

    /** @TODO
    Don't need this method
     **/
    void allocateExtra() override;

    /** @TODO
    Don't need this method
     **/
    void computeEquations() override;

    /** @TODO
    Don't need this method
     **/
    void computeDistances() override;

    /**
    find gcd for two number
    **/
    int gcd(int a, int b) ;

    //Discretisation of height
    int height_div;
    //Discretisation of circumference
    int angle_div;
    //image of unrolled
    cv::Mat normalizedMap;
    //unrolled surface representation : each cells contain many points. Order by z and theta
    std::vector < std::vector < std::vector<unsigned int> > > unrolled_surface;
};

#endif
