#ifndef SURFACE_ANALYSE_H
#define SURFACE_ANALYSE_H

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

//#define BIN_SIZE 0.01
using namespace DGtal;




class DefectSegmentation : public SegmentationAbstract {
    public:
        //DefectSegmentation(std::vector<Z3i::RealPoint> &aCloud, std::vector<Z3i::RealPoint> &aFib, double arcLen, int wh):
        //               pointCloud(aCloud), fiber(aFib), arcLength(arcLen), patchHeight(wh){
        //}
        using SegmentationAbstract::SegmentationAbstract;

        void init() override;

        std::pair<double, double> getCoeffs(unsigned int index);
        std::pair<double, double> getCoeffs2(unsigned int index);

    protected:
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

        //coefficients of regressed lines, one line for each windows
        //a window = some bands consecutives
        std::vector<std::pair<double, double> > coefficients;
 };
#endif
