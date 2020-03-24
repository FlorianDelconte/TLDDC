#ifndef REGRESSION_H
#define REGRESSION_H

#include <utility>
#include <vector>
#include <cassert>
#include <iostream>

#include <eigen3/Eigen/Dense>

#include "Statistic.h"
#include <gsl/gsl_multifit.h>


#define MAX_SD 1000
#define MIN_NB_POINTS 10
#define EPS 1

struct coefs
{
    std::pair<double, double> coefficients;
    std::vector<unsigned int> ind_p;
};

class Regression
{
public:
    Regression(){}

    static std::pair<double, double> robustLinearOls(const std::vector<double> &xs, const std::vector<double> &ys){
        const size_t p = 2; /* linear fit */
        gsl_matrix *X, *cov;
        gsl_vector *x, *y, *c;
        size_t n = ys.size();
        X = gsl_matrix_alloc (n, p);
        x = gsl_vector_alloc (n);
        y = gsl_vector_alloc (n);

        c = gsl_vector_alloc (p);
        cov = gsl_matrix_alloc (p, p);

        for (size_t i = 0; i < n; i++){
            gsl_vector_set (x, i, xs[i]);
            gsl_vector_set (y, i, ys[i]);
            gsl_matrix_set (X, i, 0, 1.0);
            gsl_matrix_set (X, i, 1, xs[i]);
        }

        dofit(gsl_multifit_robust_ols, X, y, c, cov);
        #define C(i) (gsl_vector_get(c,(i)))
        std::pair<double, double> coeffs;
        coeffs.first = C(1);
        coeffs.second = C(0);

        gsl_matrix_free (X);
        gsl_vector_free (x);
        gsl_vector_free (y);

        gsl_vector_free (c);
        gsl_matrix_free (cov);

        return coeffs;
    }
    static std::pair<double, double> robustLinear(const std::vector<double> &xs, const std::vector<double> &ys){
        const size_t p = 2; /* linear fit */
        gsl_matrix *X, *cov;
        gsl_vector *x, *y, *c;
        size_t n = ys.size();
        X = gsl_matrix_alloc (n, p);
        x = gsl_vector_alloc (n);
        y = gsl_vector_alloc (n);

        c = gsl_vector_alloc (p);
        cov = gsl_matrix_alloc (p, p);

        for (size_t i = 0; i < n; i++){
            gsl_vector_set (x, i, xs[i]);
            gsl_vector_set (y, i, ys[i]);
            gsl_matrix_set (X, i, 0, 1.0);
            gsl_matrix_set (X, i, 1, xs[i]);
        }

        dofit(gsl_multifit_robust_fair, X, y, c, cov);
        #define C(i) (gsl_vector_get(c,(i)))
        std::pair<double, double> coeffs;
        coeffs.first = C(1);
        coeffs.second = C(0);

        gsl_matrix_free (X);
        gsl_vector_free (x);
        gsl_vector_free (y);

        gsl_vector_free (c);
        gsl_matrix_free (cov);

        return coeffs;
    }

    /*static std::pair<double, double> linearRegression(const std::vector<double> &xs, const std::vector<double> &ys){
        std::pair<double, double> coefficients;
        size_t N = ys.size();
        assert(xs.size() == N);
        //too little infos
        if( N < MIN_NB_POINTS ){
            return coefficients;
        }
        double mean = Statistic::getMean(ys);
        double sd = Statistic::standardDeviation(ys, mean);
        //double mediane = Statistic::standardDeviation(ys);
///std::cout<<"bf:"<< xs.size()<< "  "<< mean<< "  "<< sd << std::endl;
        std::vector<double> inXs;
        std::vector<double> inYs;
        //std::cout<<"sd:"<<sd<<std::endl;
        double th = mean + 2*sd;
        double th2 = mean - 2*sd;
        for(size_t i = 0; i < N; i++)
        {
            //if(data[i][1] < upThreshold)
            if(ys[i] < th)
            {
                inXs.push_back(xs[i]);
                inYs.push_back(ys[i]);
            }
        }
//std::cout<<"at:" << inXs.size()<< "  "<< mean<< "  "<< sd << std::endl;
        if(sd < MAX_SD){
            return rmse(inXs, inYs);
            //robustLinearOls(xs, ys);
        }

//std::cout<<"ransac:"<<std::endl;
        return ransac(inXs, inYs, 2, 100);

    }*/
    static struct coefs PurgedByMedianlinearRegression(const std::vector<double> &xs, const std::vector<double> &ys, const std::vector<unsigned int> &indP){
        struct coefs c;
        //std::pair<double, double> coefficients;
        unsigned int  N = ys.size();
        assert(xs.size() == N);
        //too little infos
        if( N < MIN_NB_POINTS ){
            return c;
        }
        double mean = Statistic::getMean(ys);
        double sd = Statistic::standardDeviation(ys, mean);
        std::vector<double> inXs;
        std::vector<double> inYs;
        std::vector<unsigned int> indPselect;
        //filtre sur la moyenne
        double th = mean + 2*sd;
        double th2 = mean - 2*sd;
        for(unsigned int  i = 0; i < N; i++)
        {
            if(ys[i] < th)
            {
                indPselect.push_back(indP.at(i));
                inXs.push_back(xs[i]);
                inYs.push_back(ys[i]);
            }
        }
        //filtre by medianne
        std::vector<double> inXs2;
        std::vector<double> inYs2;
        std::vector<unsigned int> indPselect2;
        unsigned int T=inYs.size();
        double mediane=Statistic::getMedian(inYs);
        for(unsigned int i = 0; i < T; i++)
        {
            if(inYs[i] < mediane   )
            {
                indPselect2.push_back(indPselect[i]);
                inXs2.push_back(inXs[i]);
                inYs2.push_back(inYs[i]);
            }
        }
        c.ind_p=indPselect2;
        if(sd < MAX_SD){
            c.coefficients=rmse(inXs2, inYs2);
        }else{
            c.coefficients=ransac(inXs2, inYs2, 2, 100);
        }
        return c;
    }

    static std::pair<double, double> ransac(const std::vector<double> xs, std::vector<double> ys, double epsilon, int minNbIter){
        //use pcl instead of mine?
        std::pair<double, double> coefficients;
        //struct coefs c;
    srand (time(NULL));
    std::pair<double, double> coeffrmse = rmse(xs, ys);

    int numberOfIteration = std::numeric_limits<int>::max();
    double p = 0.9999;
    int n = xs.size();

    int nbSamples = 0;
    int bestNbInliers = 0;
    while (nbSamples < numberOfIteration)
    {
        //random select 2 points
        int p1 = n*(1.0*rand() / RAND_MAX);
        int p2 = n*(1.0*rand() / RAND_MAX);
        if(p1 == p2){
            continue;
        }
        double y1 = ys.at(p1);
        double y2 = ys.at(p2);

        double x1 = xs.at(p1);
        double x2 = xs.at(p2);

        double a = (y1 - y2)/(x1 - x2);
        double b = y1 - a*x1;

        if( std::abs(a - coeffrmse.first) > 0.5 ){
            continue;
        }
        
        int nbInliers = 0;
        //calculate the number of inliers
        for (unsigned int i = 0; i < xs.size(); i++)
        {
            double x = xs.at(i);
            double y = ys.at(i);
            if(std::abs(a * x + b - y) < epsilon){
                nbInliers++;
                //inliers.push_back(i);
            }
        }

        if(nbInliers > bestNbInliers)
        {
            bestNbInliers = nbInliers;
            coefficients.first = a;
            coefficients.second = b;
        }

        double w = 1.0*bestNbInliers/xs.size();
        numberOfIteration = log10(1 - p)/log10(1 - pow(w, 2));
        if(numberOfIteration < minNbIter) numberOfIteration = minNbIter;
        nbSamples++;
        //qDebug()<<"Nb Iter: "<< nbSamples;
        //qDebug()<<"Nb Iterxxx: "<< numberOfIteration;
        //qDebug()<<"Max Nb Inliers /n : "<< bestNbInliers << "/"<<n;

    }

    //get inliers
    std::vector<double> inXs;
    std::vector<double> inYs;
    for (unsigned int i = 0; i < xs.size(); i++)
    {
            double x = xs.at(i);
            double y = ys.at(i);
            if(std::abs(coefficients.first * x + coefficients.second - y) < epsilon){
                inXs.push_back(x);
                inYs.push_back(y);
            }
    }

    return rmse(inXs, inYs);

    }
    static std::pair<double, double>  rmse(const std::vector<double> &xs, const std::vector<double> &ys){
        //struct coefs c;
        std::pair<double, double> coefficients;
        size_t N = xs.size();

        //least square fitting with inliers
        Eigen::MatrixX2d A(N,2);
        Eigen::VectorXd R(N);

        for (size_t i = 0; i < N; i++){
            double x = xs.at(i);
            double y = ys.at(i);
            A(i, 0) = x;
            A(i, 1) = 1;
            R[i] = y;
        }

        Eigen::Vector2d P = A.colPivHouseholderQr().solve(R);
        coefficients.first = P[0];
        coefficients.second = P[1];
        return coefficients;
    }
private:
    static int
    dofit(const gsl_multifit_robust_type *T,
          const gsl_matrix *X, const gsl_vector *y,
          gsl_vector *c, gsl_matrix *cov){
        int s;
        gsl_multifit_robust_workspace * work 
            = gsl_multifit_robust_alloc (T, X->size1, X->size2);
        //gsl_multifit_robust_maxiter(100000, work);
        //gsl_multifit_robust_tune(4, work);
        s = gsl_multifit_robust (X, y, c, cov, work);
        gsl_multifit_robust_free (work);
        return s;
    }

};
#endif // REGRESSION_H
