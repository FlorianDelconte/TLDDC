#ifndef DEFECT_SEGMENTATION_CYLINDER_H
#define DEFECT_SEGMENTATION_CYLINDER_H

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

#include <pcl/filters/extract_indices.h>
#include <pcl/features/normal_3d.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>



#include "SegmentationAbstract.h"
#include "CylindricalPoint.h"

//#define BIN_SIZE 0.01
using namespace DGtal;


typedef pcl::PointXYZ PointT;

struct CylindricalCoefficients{
    Z3i::RealPoint point;
    Z3i::RealPoint direction;
    double radius;
};

class DefectSegmentationCylinder : public SegmentationAbstract {
    public:
        using SegmentationAbstract::SegmentationAbstract;

        void init() override;
        std::vector<CylindricalCoefficients> getCoeffs();

        std::vector<std::pair<Z3i::RealPoint, Z3i::RealPoint> > getExtremityOfCylinders();
    protected:
        std::vector<CylindricalCoefficients> coefficients;

                /** Brief 
         ** Allocate (resize) memory for array
         **/
        void allocateExtra() override;
        /** Brief 
         * Compute the equation of the relation between distance to center line and z
         * of patches associated to points
         */
        void computeEquations() override;

        /** Brief 
         * Call by computeEquations
         */
        void computeEquationsMultiThread(int threadId, int nbThread, const pcl::KdTreeFLANN<pcl::PointXYZ> &kdtree, 
                             double minHeight, double maxHeight);

        void computeDistances() override;

        std::pair<Z3i::RealPoint, Z3i::RealPoint> 
        getExtremityOfCylinder(const std::vector<Z3i::RealPoint> &pointCloud, 
                       const std::vector<unsigned int> &pointsOnCylinder, 
                       Z3i::RealPoint center1, Z3i::RealPoint center2);

        std::vector<std::pair<Z3i::RealPoint, Z3i::RealPoint> > extremities;
};
#endif
