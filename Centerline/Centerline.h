#ifndef CENTERLINE_H
#define CENTERLINE_H 

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/writers/VolWriter.h"
#include <DGtal/base/BasicFunctors.h>
#include <DGtal/shapes/implicit/ImplicitBall.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/shapes/EuclideanShapesDecorator.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>

#include "DGtal/shapes/Mesh.h"



///////////////////////////////////////////////////////////////////////////////
// class Centerline
/**
 * Description of class 'Centerline' <p>
 * 
 * @brief Class to compute centerline of log using volumic accumulation from normal vector field
 * and Splines to smoothing
 *
 *
 */

using namespace DGtal;

// types of image containers:
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<DGtal::Z3i::Domain, double> Image3DDouble;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

typedef ConstImageAdapter<Image3D, Z2i::Domain, functors::Point2DEmbedderIn3D<Z3i::Domain>,
                                 unsigned int, functors::Identity >  ImageAdapterExtractor;


class Centerline{

// ----------------------- Standard methods ------------------------------
public:
    Centerline( const Mesh<Z3i::RealPoint> &aMesh, const double aRadius, double step, bool iNormal):  
        mesh(aMesh),
        accRadius(aRadius), 
        trackStep(step),
        invertNormal(iNormal),
        accImage(Z3i::Domain()), dirImage(Z3i::Domain()){
            std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> boudingBox = mesh.getBoundingBox();
            Z3i::RealPoint ptLow = boudingBox.first;
            Z3i::RealPoint ptUp = boudingBox.second;
            domain = Z3i::Domain(Z3i::Point((int) ptLow[0], (int) ptLow[1], (int) ptLow[2]),
                    Z3i::Point((int) ptUp[0], (int) ptUp[1], (int) ptUp[2]));
            accImage = Image3D(domain);
            dirImage = ImageVector(domain);
        }

    std::vector<Z3i::RealPoint> compute();

//protected functions
protected:

    // Optimize fiber according sections vertex
    std::vector<Z3i::RealPoint>
    optimizeElasticForces(std::vector<Z3i::RealPoint> aFiberRaw, double epsilon);

    // Used to detect end of tracking
    bool
    isFurtherInside(const Z3i::RealPoint &aPoint,
                const Z3i::RealPoint &aPreviousPoint, double aDistance );

    // Track patch center in one direction 
    std::vector<Z3i::RealPoint>
    trackPatchCenter(const Z3i::Point &aStartingPoint, bool firstDirection);

    /**
     * Main method for tracking skeleton points from image. 
     * Tracking the centerline from the maximum density point by two directions
     * @return raw centerline
     */
    std::vector<Z3i::RealPoint>
    trackCenterline(const Z3i::Point &aStartingPoint);

    /**
     * Main method for computing tube accumulation values.  
     * It also computes the vector image representing the main
     * axis directionfrom vectorial products from the intersection map (to
     * help the tracking process).  
     * @return the maximum accumulation point.
     */
    Z3i::Point
    accumulate(double epsilonArea);

    // protected attributes:
protected:
    Mesh<Z3i::RealPoint> mesh;
    double accRadius; // the maximal radius of accumulation
    Z3i::Domain domain; // the domain of the mesh

    Image3D accImage;  //store accumulation value
    ImageVector dirImage; //store direction
    double trackStep;	//the distance between 2 steps using by tracking algo
    bool invertNormal;

};

#endif //end CENTERLINE_H

