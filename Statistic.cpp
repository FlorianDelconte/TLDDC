#include "Statistic.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>

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

double Statistic::getMedian(std::vector<double> v){
    double median=0.0;
    unsigned int size = v.size();
    if (size == 0){
        median=0.;  // Undefined, really.
    }else{
        std::sort(v.begin(), v.end());
        if (size % 2 == 0){
            median=(v.at(size / 2 - 1) + v.at(size / 2)) / 2;
        }else{
            median=v.at(size / 2);
        }
    }
    return median;
}
