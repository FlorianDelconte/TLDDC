#ifndef CENTERLINE_HELPER_H
#define CENTERLINE_HELPER_H
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <numeric>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"
#include "DGtal/shapes/GaussDigitizer.h"


//#include "DGtal/geometry/curves/AlphaThickSegmentComputer.h"

class CenterlineHelper{


public:

    typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned int> Image3D;
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> Image3DChar;
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  DGtal::Z3i::RealPoint> ImageVector;
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  std::vector<DGtal::Z3i::RealPoint> > ImagePointAssociation;


    template<typename TPoint>
    static std::vector<TPoint>
    getSmoothCenterlineBSplines(const DGtal::Z3i::Domain &aDomain, const std::vector<TPoint> &fib){
        double minAngle = M_PI / 4*3;

        std::vector<TPoint> fibOut;
        //check if first point is good
        TPoint v1 = fib.at(0) - fib.at(1);
        TPoint v2 = fib.at(2) - fib.at(1);
        double cosA = v1.dot(v2)/v1.norm()/v2.norm();
        double angle = acos(cosA);
        unsigned start = 1;

        if(angle > minAngle){
            fibOut.push_back(fib.at(0));
        }else{
            fibOut.push_back(fib[1]);
            start = 2;
        }

        for (unsigned int i = start; i < fib.size() - 2; i++){
            TPoint v1 = fibOut.at(fibOut.size() -1) - fib.at(i);
            TPoint v2 = fib.at(i+1) - fib.at(i);
            double cosA = v1.dot(v2)/v1.norm()/v2.norm();
            double angle = acos(cosA);
            if(angle > minAngle){
            fibOut.push_back(fib.at(i));
            }
        }
        fibOut.push_back(fib.at(fib.size() - 2));

        //splines after remove point with detours
        //double meanLength = std::accumulate(fibOut.begin(), fibOut.end(), 0.0)/fibOut.size();
        std::vector<double> ds(fibOut.size(), 0);
        for (unsigned int i = 0; i < fibOut.size() - 1; i++){
            ds[i] = (fibOut.at(i+1) - fibOut.at(i)).norm();
        }

        double meanDs = std::accumulate(ds.begin(), ds.end(),0.0)/ds.size();
        std::vector<TPoint> goodFib;
        for (unsigned int i = 1; i < fibOut.size(); i++){
            if(ds.at(i) > meanDs / 2){
                goodFib.push_back(fibOut.at(i));
            }
        }

        goodFib.push_back(fibOut.at(fibOut.size() - 1));

        unsigned int nbGoodPoint = goodFib.size();


        unsigned int nbPointForSample = 8;
        if(nbPointForSample > nbGoodPoint){
            nbPointForSample = nbGoodPoint;
        }

        std::vector<TPoint> sample(nbPointForSample);
        sample[0] = goodFib.at(0);
        sample[nbPointForSample - 1] = goodFib.at(goodFib.size() -1);
        double step = goodFib.size()*1.0 / (nbPointForSample - 1);
        //trace.info()<<"step:"<<step<<std::endl;
        unsigned int index = 1;
        while(index < nbPointForSample - 1){
            unsigned int i = index * step;
            sample[index++] = goodFib.at(i);
        }

        //Image3D anImage(domain);
        double* T = (double*)malloc(nbPointForSample*sizeof(double));
        double* X = (double*)malloc(nbPointForSample*sizeof(double));
        double* Y = (double*)malloc(nbPointForSample*sizeof(double));
        double* Z = (double*)malloc(nbPointForSample*sizeof(double));
        //splines after remove point with detours
        for(unsigned int i = 0; i< nbPointForSample; i++){
            T[i] = i;
            X[i] = sample[i][0];
            Y[i] = sample[i][1];
            Z[i] = sample[i][2];
        }

        gsl_interp_accel *accX = gsl_interp_accel_alloc ();
        gsl_spline *splineX = gsl_spline_alloc (gsl_interp_cspline, nbPointForSample);

        gsl_interp_accel *accY = gsl_interp_accel_alloc ();
        gsl_spline *splineY = gsl_spline_alloc (gsl_interp_cspline, nbPointForSample);

        gsl_interp_accel *accZ = gsl_interp_accel_alloc ();
        gsl_spline *splineZ = gsl_spline_alloc (gsl_interp_cspline, nbPointForSample);

        gsl_spline_init (splineX, T, X, nbPointForSample);
        gsl_spline_init (splineY, T, Y, nbPointForSample);
        gsl_spline_init (splineZ, T, Z, nbPointForSample);

        std::vector<TPoint> smoothFib;
        for(double t = 0; t <= nbPointForSample - 1; t += 0.025){
            double x = gsl_spline_eval(splineX, t, accX);
            double y = gsl_spline_eval(splineY, t, accY);
            double z = gsl_spline_eval(splineZ, t, accZ);
            TPoint sPoint(x, y, z);
            smoothFib.push_back(sPoint);
        }

        gsl_spline_free (splineX);
        gsl_interp_accel_free (accX);

        gsl_spline_free (splineY);
        gsl_interp_accel_free (accY);

        gsl_spline_free (splineZ);
        gsl_interp_accel_free (accZ);

        free(T);
        free(X);
        free(Y);
        free(Z);

        TPoint firstVect = smoothFib[0] - smoothFib[1];
        TPoint firstPoint = smoothFib[0];

        TPoint lastVect = smoothFib[smoothFib.size()-1] - smoothFib[smoothFib.size()-2];
        TPoint lastPoint = smoothFib[smoothFib.size()-1];

        std::vector<TPoint> frontVect;

        TPoint nextFront = firstPoint + firstVect;

        while(aDomain.isInside(nextFront)){
            frontVect.push_back(nextFront);
            nextFront = nextFront + firstVect;
        }

        std::vector<TPoint> backVect;
        TPoint nextBack = lastPoint + lastVect;
        while(aDomain.isInside(nextBack)){
            backVect.push_back(nextBack);
            nextBack = nextBack + lastVect;
        }

        std::reverse(frontVect.begin(), frontVect.end());

        frontVect.insert(frontVect.end(), smoothFib.begin(), smoothFib.end());
        frontVect.insert(frontVect.end(), backVect.begin(), backVect.end());
        smoothFib = frontVect;

        return smoothFib;

    }

    /**
     * Get all the mesh faces as indice associated to a fiber point.
     *  => Method by projecting in the given direction.
     *
     *  @param aMesh the source mesh
     *  @param aFiberPt the fiber for which we want the sections
     *  @param aDirection the main axis direction vector
     *  @param aSectionSize the size of the section according the main axis
     *  @param radius: use to eliminate faces outside the mesh.
     **/

    template<typename TPoint>
    static std::vector<unsigned int>
    getSectionFacesFromDirection(const DGtal::Mesh<TPoint> &aMesh, const TPoint &aFiberPt,
                                 const TPoint &aDirection, double aSectionSize, double radius){

        typedef typename DGtal::Mesh<TPoint>::MeshFace Face;
        std::vector<unsigned int>  vectResult;

        for(unsigned int j = 0; j<aMesh.nbFaces(); j++){

            Face aFace = aMesh.getFace(j);
            TPoint center;
            center[0]=0; center[1]=0; center[2]=0;
            for (unsigned int i = 0; i< aFace.size(); i++){
                TPoint aPoint  = aMesh.getVertex(aFace.at(i));
                center+=aPoint;
            }
            center /= (double) aFace.size();
            // Projection on the directional vector (aDirection)
            TPoint vectCenter = center - aFiberPt;
            double fact = std::abs(aDirection.dot(vectCenter));
            double radDist  = std::abs(vectCenter.crossProduct(aDirection).norm()/aDirection.norm());
            if (fact <= aSectionSize && radDist < radius ){
                vectResult.push_back(j);
            }
        }
        return vectResult;
    }




    template<typename TImage>
    static DGtal::Z2i::Point
    getMaxCoords(const TImage &anImage){
        DGtal::Z2i::Point ptMax = *(anImage.domain().begin());
        typename TImage::Value valMax = anImage(ptMax);
        for( typename TImage::Domain::ConstIterator it = anImage.domain().begin(); it!= anImage.domain().end(); it++){
            typename TImage::Value val = anImage(*it);
            if(val>valMax){
                valMax = val;
                ptMax = *it;
            }
        }
        return ptMax;
    }





    template<typename TImage>
    static void
    getBallOrientedSurfaceSet(TImage anImage, std::vector<DGtal::Z3i::Point> &aVectPoint,
                const DGtal::Z3i::RealPoint &aPoint, const DGtal::Z3i::RealPoint aPreviousPoint,
                double aRadius, bool filterOrientation, int threshold=200 ){

        typedef DGtal::ImplicitBall<DGtal::Z3i::Space> EuclideanBall;
        typedef DGtal::deprecated::EuclideanShapesMinus< EuclideanBall, EuclideanBall > Minus;
        typedef DGtal::GaussDigitizer<DGtal::Z3i::Space, Minus> DigitalShape;
        EuclideanBall aBallp (DGtal::Z3i::Point(aPoint[0], aPoint[1], aPoint[2]), aRadius+1);
        EuclideanBall aBallm (DGtal::Z3i::Point(aPoint[0], aPoint[1], aPoint[2]), aRadius-1);
        Minus aBallSurface ( aBallp, aBallm );

        DigitalShape  gaussDig;
        gaussDig.attach (aBallSurface);
        gaussDig.init(DGtal::Z3i::Point(aPoint[0]-(int)aRadius-1, aPoint[1]-(int)aRadius-1, aPoint[2]-(int)aRadius-1),
                DGtal::Z3i::Point(aPoint[0]+(int)aRadius+1, aPoint[1]+(int)aRadius+1, aPoint[2]+(int)aRadius+1), 1);
        DGtal::Z3i::Domain dom = gaussDig.getDomain();
        if(aPoint==aPreviousPoint){
            DGtal::trace.info() << "Point identique...." <<std::endl;
            return;
        }
        DGtal::Z3i::RealPoint dirRefNormalized = (aPoint - aPreviousPoint)/(aPoint - aPreviousPoint).norm();
        for( DGtal::Z3i::Domain::ConstIterator it = dom.begin(); it!=dom.end(); it++){
            if( *it != aPoint){
                DGtal::Z3i::RealPoint dirCurrentNormalized = (*it-aPoint)/(*it-aPoint).norm();
                if(gaussDig(*it) && anImage.domain().isInside(*it) &&
                        (!filterOrientation || dirRefNormalized.dot(dirCurrentNormalized)>=0.5)
                        && anImage(*it)>= threshold){
                    aVectPoint.push_back(*it);
                }
            }
        }
    }

};
#endif
