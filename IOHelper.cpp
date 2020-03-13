#include <iostream>
#include <fstream>

//#include <opencv2/opencv.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"


#include "IOHelper.h"


using namespace DGtal;

void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        Z3i::RealPoint p = v[i];
        outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
    }
    outStream.close();
}

void IOHelper::export2Text(const std::vector<DGtal::Z3i::RealPoint> &pointCloud, 
        const std::vector<unsigned int> &indices, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < indices.size();i++){
        Z3i::RealPoint p = pointCloud.at(indices.at(i));
        outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
    }
    outStream.close();
}
void IOHelper::export2Text(const std::vector<std::pair<double, double> > &v, const std::string &filename){
    std::ofstream outStream;
    outStream.open(filename.c_str(), std::ofstream::out);
    for(unsigned int i = 0; i < v.size();i++){
        std::pair<double, double> pa = v.at(i);
        outStream << pa.first << " "<< pa.second <<std::endl;
    }
    outStream.close();
}

void IOHelper::readDistanceFromFile(const std::string &fileName, std::vector<double> &vectDistances){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
      if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
          vectDistances.push_back(std::stod(str));
      }
      getline(infile, str);
    }
  }


void IOHelper::export2OFF(const Mesh<Z3i::RealPoint> &mesh, std::string fileName){
    std::ofstream offMesh (fileName.c_str());
    DGtal::MeshWriter<Z3i::RealPoint>::export2OFF(offMesh, mesh);
    offMesh.close();
}


void IOHelper::readIntsFromFile(const std::string &fileName, std::vector<int> &rs){
    std::ifstream infile;
    infile.open (fileName.c_str(), std::ifstream::in);
    std::string str;
    getline(infile, str );
    while ( infile.good() ){
        if ( ( str != "" ) && ( str[ 0 ] != '#' ) ){
            rs.push_back(std::stoi(str));
        }
        getline(infile, str);
    }
}

/*void IOHelper::createVisuImage(std::string s,cv::Mat c){
    int grayscaleValue;
    double normalizedValue;
    int rows = c.rows;
    int cols = c.cols;
    cv::Mat grayscalemap(rows,cols,CV_8UC1,cv::Scalar(0));
    cv::Mat reliefPictures(rows, cols, CV_8UC3, cv::Scalar(110, 110, 110));

    for(unsigned int i = 0; i < rows; i++){
        for(unsigned int j = 0; j < cols; j++){
            normalizedValue=c.at<double>(i, j);
            grayscaleValue=((255/1)*(normalizedValue-1))+255;
            grayscalemap.at<uchar>(i, j) = grayscaleValue;
        }
    }
    cv::applyColorMap(grayscalemap, reliefPictures, cv::COLORMAP_JET);
    imwrite( "unrollSurfacePictures/"+s, reliefPictures);
}*/
