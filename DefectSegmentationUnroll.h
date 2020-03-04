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
    compute normalized image from unrolled surface with max resolution
    **/
    void computeNormalizedImage();
    /**
    compute normalized image from unrolled surface with a decrease factor specif
    **/
    void computeNormalizedImage(int dF);
    /**
    create visual multi-scaled image from unrolled surface
    **/
    void computeNormalizedImageMultiScale();
    /**
    create visual image from normalized image
    **/
    void createVisuImage(std::string s);
    /**
    create GT image from vector of int (ind of point in groung truth file)
    **/
    void createGTImage(std::vector<int>,std::string s);
  protected:
   /**
    return normalized radius in [0;1] 
    **/
    double normalizeRadius(double value);
    /**
    return normalized difference of distance (patch distance - distance) in [0;1] value 
    **/
    double normalizeDiffDistance(double value);
    /**
    return a vector of Points
    **/
    std::vector<unsigned int > getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF);
     /**
    return false if cells is consiedred as out of the mesh (function to skip noise in relief image)
    **/
    bool detectCellsIn(unsigned int i, unsigned int j);
    /**
    Cylindrical point from segment of centerline
    **/
    std::vector<unsigned int >  getPointfromSegmentId(unsigned int id);
    /**
    return the mean of distance difference vector given in parameters
    **/
    double getMeansDistDiff(std::vector<unsigned int > v);
    /**
    return the max of distance difference vector given in parameters
    **/
    double getMaxDistDiff(std::vector<unsigned int > v);
    /**
    return the mean (on radius) of vector
    **/
    double getMeansRadius(std::vector<unsigned int > v);
    /**
    return the max (on radius) of vector
    **/
    double getMaxRadius(std::vector<unsigned int > v);
    /**
    return the median (on radius) of vector
    **/
    double getMedianRadius(std::vector<unsigned int > v);
    /**
    Method to find the discretisation i / j of angle and height
    **/
    void computeMeasures();
    /**
    transfrome a 2Darray of vector into an image
    **/
    void transformToImage( std::vector<unsigned int > **two_d);
    /**
     * Allocate space for unrollMapp
     **/
    void allocateExtra() override;
    /** @TODO
    Refactor with Van tho'code
     **/
    void computeEquations() override;
    /** Brief 
    Refactor with Van tho'code
    **/
    void computeEquationsMultiThread(int threadId, int nbThread, const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, double minHeight, double maxHeight);
    /** @TODO
    Refactor with Van tho'code
     **/
    void computeDistances() override;

    //minimal angle on cylindrical coordonate of all points
    double minAngle;
    //maximal angle on cylindrical coordonate of all points 
    double maxAngle;
    //minimal radius on cylindrical coordonate of all points
    double minRadius;
    //maximal radius on cylindrical coordonate of all points 
    double maxRadius;
    //mean of radius on cylindrical coordonate of all points
    double meanRadius;
    //minimal height on cylindrical coordonate of all points
    double minHeight;
    //maximal height on cylindrical coordonate of all points
    double maxHeight;
    //Discretisation of height
    int height_div;
    //maximal difference of distance
    double maxdiffDist;
    //minimal difference of distance
    double mindiffDist;
    //Discretisation of circumference
    int angle_div;
    //image of unrolled
    cv::Mat normalizedMap;
    //unrolled surface representation : each cells contain many points. Order by z and theta
    std::vector < std::vector < std::vector<unsigned int> > > unrolled_surface;
    //@TODO: refactoring with van Tho' code.
    //coefficients of regressed lines, one line for each windows a window = some bands consecutives
    std::vector<std::pair<double, double> > coefficients;
};

#endif
