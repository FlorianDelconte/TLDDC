#include "Statistic.h"
#include <stdlib.h>
#include <math.h>

Statistic::Statistic(){
}

double Statistic::standardDeviation(const std::vector<double> &v, const double mean){
    double s = 0.0;
    for (size_t i = 0; i < v.size(); i++)
    {
        s += (v[i] - mean) * (v[i] - mean);
    }
    return sqrt(s/(v.size() - 1));
}

double Statistic::getMean(const std::vector<double> &v){
    //calculate mean
    double sum = 0.0;
    for (unsigned int i = 0; i < v.size(); i++)
    {
        sum += v[i];
    }

    return sum / v.size();
}

//double Statistic::getMedian(std::vector<double> v){
//    return 0.0;
//}
