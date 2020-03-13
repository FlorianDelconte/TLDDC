#ifndef UNROLL_SURFACE2
#define UNROLL_SURFACE2

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

class DefectSegmentationUnroll_copy : public SegmentationAbstract {
  public:

    using SegmentationAbstract::SegmentationAbstract;

    void init() override;
    void getDefect(std::string output);

  protected:
    /**
    Write a RGB image
    **/
    void createVisuImage(std::string s,cv::Mat c);
    /**
    Compute line coeficient for one id point
    **/
    std::pair<double, double> computeEq(unsigned int idPoint,double searchRadius, double patchAngle,const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree);
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
    //@TODO: refactoring with van Tho' code.
    //coefficients of regressed lines, one line for each windows a window = some bands consecutives
    std::vector<std::pair<double, double> > coefficients;
};

#endif
