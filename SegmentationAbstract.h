#ifndef SURFACE_ANALYSE_ABSTRACT_H
#define SURFACE_ANALYSE_ABSTRACT_H

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


#include "CylindricalPoint.h"


using namespace DGtal;


typedef unsigned int uint;

struct CylindricalPointOrder {
    bool operator() (CylindricalPoint p1, CylindricalPoint p2);
};
struct CylindricalPointOrderRadius {
    bool operator() (CylindricalPoint p1, CylindricalPoint p2);
};

struct CylindricalPointOrderAngle {
    bool operator() (CylindricalPoint p1, CylindricalPoint p2);
};

struct CylindricalPointOrderArcLenght {
    bool operator() (CylindricalPoint p1, CylindricalPoint p2);
};


class SegmentationAbstract{
    public:
        SegmentationAbstract(std::vector<Z3i::RealPoint> &aCloud, std::vector<Z3i::RealPoint> &aFib, double arcLen, int wh, double bw):
                       pointCloud(aCloud), fiber(aFib), arcLength(arcLen), patchHeight(wh), binWidth(bw){
        }

        virtual void init() = 0;
        std::vector<unsigned int> getDefect();

        std::vector<unsigned int> getDefect(double threshold);

        std::vector<double> getDistances();
        std::vector<std::vector<unsigned int> > getCells();

        double getRadius(unsigned int index);
        double getLength(unsigned int index);

        //for demo
        std::vector<unsigned int> getPatch(unsigned int pointId);
        CylindricalPoint getPointInCylindric(unsigned int pId);

        double findThresholdRosin();


        int getNbSegment();
        int getNbSector();

        /** Brief
         * get the segment that aPoint belongs to
         */
        unsigned int getSegment(const Z3i::RealPoint &aPoint);

    protected:
        /** Brief
         ** Allocate (resize) memory for array
         **/
        void allocate();

        virtual void allocateExtra() = 0;
        /** Brief
         * Compute the equation of the relation between distance to center line and z
         * of patches associated to points
         */
        virtual void computeEquations() = 0;
        void computeAngleOfPoints();
        void computeSegmentsAndRadii();
        void computeBeginOfSegment();
        void computeSegments();
		    void convertToCcs();

        //should be change to computeLocalCoordinate vectors
        void computeVectorMarks();
        virtual void computeDistances() = 0;
        void computePlaneNormals();

        Z3i::RealPoint
        getRadialVector(const Z3i::RealPoint &aPoint, const Z3i::RealPoint &aDirection, const Z3i::RealPoint &p0);

        /** Brief
         * get the sector that aPoint belongs to
         */
        unsigned int
        getSector(const Z3i::RealPoint &aPoint);

        Z3i::RealPoint
        getDirectionVector(const unsigned int &segmentId);

        void writeDebugInfo();
        /////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<Z3i::RealPoint> &pointCloud;
        std::vector<Z3i::RealPoint> &fiber;

        //band width
        double arcLength;
        //number of band used to estimate equation
        int patchHeight;

        //the number of angular sector
        int nbSector;

        //the number of segment
        int nbSegment;
        //store segment id of each point
        //@TOTO:change to char
        std::vector<CylindricalPoint> myPoints;
        //store distance of each point to reference surface
        //@TOTO:change to char

        std::vector<double> beginOfSegment;
        std::vector<double> segmentLength;

        //normal vector of the plane passing Oi and divide the log
        std::vector<Z3i::RealPoint> ns;
        //local vector Oy
        std::vector<Z3i::RealPoint> vectMarks;

        //coefficients of regressed lines, one line for each windows
        //a window = some bands consecutives
        //
        //difference between
        std::vector<double> distances;
        //bin width used by Rosin method
        double binWidth;

        double radii;
};
#endif
