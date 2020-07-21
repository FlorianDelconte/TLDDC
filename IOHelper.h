#ifndef IO_HELPER_H
#define IO_HELPER_H

#include <iostream>
#include <utility>

//#include <opencv2/opencv.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/shapes/Mesh.h"



using namespace DGtal;

class IOHelper{
public:
    IOHelper(){
    }
    template<typename T, typename = typename
            std::enable_if<std::is_arithmetic<T>::value>::type>
    static void export2Text(const std::vector<T> &v, const std::string &filename){
        std::ofstream outStream;
        outStream.open(filename.c_str(), std::ofstream::out);
        for(unsigned int i = 0; i < v.size();i++){
            outStream << v.at(i)<<std::endl;
        }
        outStream.close();
    }

    static void export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename);
    static void export2Text(const std::vector<std::pair<double, double> > &v, const std::string &filename);
    static void export2Text(const std::vector<DGtal::Z3i::RealPoint> &pointCloud,
            const std::vector<unsigned int> &indices, const std::string &filename);
    static void readDistanceFromFile(const std::string &fileName, std::vector<double> &vectDistances);
    //not generic!!!!
    static void export2OFF(const Mesh<Z3i::RealPoint> &mesh, std::string fileName);
    static void export2OBJ(const Mesh<Z3i::RealPoint> &mesh, std::string fileName);
    static void readIntsFromFile(const std::string &fileName, std::vector<int> &rs);
    //not generic!!!!
    //static void createVisuImage(std::string name,cv::Mat image);
};

#endif
