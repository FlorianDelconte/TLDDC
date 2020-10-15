#include <iostream>
#include <fstream>

#include <stdlib.h>

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

void IOHelper::writeDiscretisationToFile(const std::vector<std::vector<std::vector<unsigned int>>> &discretisation,const int &rowcroppedBot,const int &rowcroppedTop, const std::string &fileName){
  trace.info()<<"Writting discretisation ..."<<std::endl;
  std::ofstream outStream;
  outStream.open(fileName.c_str(), std::ofstream::out);
  //first line is for dimension X/Y of discretisation
  outStream<<discretisation.size()<<" "<<discretisation[0].size()<<std::endl;
  //Second line is for nb line cropped
  outStream<<rowcroppedBot<<" "<<rowcroppedTop<<std::endl;
  //then the discretisation
  int rows=discretisation.size();
  int cols=discretisation[0].size();

  for (int i=0; i < cols; ++i){
    for (int j=0; j < rows; ++j){
      std::vector<unsigned int> v =discretisation.at(j).at(i) ;
      if(!v.empty()){
        for(auto it = v.begin(); it != v.end(); ++it) {
          outStream << *it;
          outStream << " ";
        }
        outStream<<std::endl;
      }else{
        outStream<<"-1 "<<std::endl;
      }
    }
  }
  outStream.close();
}
void IOHelper::readDiscretisationFromFile(const std::string &fileName, std::vector<std::vector<std::vector<unsigned int>>> &discretisation, int &rowcroppedBot,int &rowcroppedTop){
  std::ifstream infile;
  infile.open(fileName.c_str(), std::ifstream::in);
  std::string currentLine;
  std::string delimiter=" ";
  std::string subLine;

  getline(infile, currentLine);
  //Here we want to read the dimension of the image. In discretisation.txt, dimension are delimited by a single space : dimY dimX. Dimension is located at the first line of the file
  int cols,rows;

  subLine=currentLine.substr(0, currentLine.find(delimiter));
  if(subLine.empty()){
    trace.info()<<"Problem in discretisation.txt. First line shoulde be : dimX dimY"<<std::endl;
  }else{
    rows=std::stoi(subLine);
  }
  subLine=currentLine.substr(currentLine.find(delimiter)+1,currentLine.size() );
  if(subLine.empty()){
    trace.info()<<"Problem in discretisation.txt. First line shoulde be : dimX dimY"<<std::endl;
  }else{
    cols=std::stoi(subLine);
  }
  //trace.info()<<"cols "<<cols<<std::endl;
  //trace.info()<<"rows "<<rows<<std::endl;
  //Here we want to read the number of rowcropped (second line of the file discretisation.txt)
  getline(infile, currentLine);
  if(currentLine.empty()){
    trace.info()<<"Problem in discretisation.txt. Second line shoulde be : numberRowCroppedFromBot numberRowCroppedFromTop"<<std::endl;
  }else{
    rowcroppedBot=std::stoi(currentLine);
  }
  subLine=currentLine.substr(currentLine.find(delimiter)+1,currentLine.size() );
  if(subLine.empty()){
    trace.info()<<"Problem in discretisation.txt. Second line shoulde be : numberRowCroppedFromBot numberRowCroppedFromTop"<<std::endl;
  }else{
    rowcroppedTop=std::stoi(subLine);
  }
  //Here we want to fill input discretisation with value in discretisation.txt starting at the trhird line
  //init dicretisation size with row and cols read before
  discretisation.resize(rows);
  for (int i = 0; i < rows; ++i){
      discretisation[i].resize(cols);
  }
  //loop on other line in the discretisation.txt
  size_t pos = 0;
  //actual position to add in discretisation
  int i=0,j=0;//i for cols and j for rows
  //there are exactly the same number of line in discretisation.txt than the number of cels un discretisation
  for (int i=0; i < cols; ++i){
    for (int j=0; j < rows; ++j){
      getline(infile, currentLine);
      std::stringstream linestream(currentLine);
      std::string value;
      while(getline(linestream,value,' ')){
        if(std::stod(value)!=-1){
          discretisation[j][i].push_back( std::stod(value));
        }
      }
    }
  }


}

void IOHelper::export2OFF(const Mesh<Z3i::RealPoint> &mesh, std::string fileName){
    std::ofstream offMesh (fileName.c_str());
    DGtal::MeshWriter<Z3i::RealPoint>::export2OFF(offMesh, mesh);
    offMesh.close();
}

void IOHelper::export2OBJ(const Mesh<Z3i::RealPoint> &mesh, std::string fileName){
    std::ofstream offMesh (fileName.c_str());
    DGtal::MeshWriter<Z3i::RealPoint>::export2OBJ(offMesh, mesh);
    offMesh.close();
}

void IOHelper::readIntsFromFile(const std::string &fileName, std::vector<int> &rs){

    std::ifstream infile;
      //std::cout<<fileName<<std::endl;
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
