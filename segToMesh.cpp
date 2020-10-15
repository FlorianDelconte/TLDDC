
#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "IOHelper.h"

#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace DGtal;
namespace po = boost::program_options;
typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

int
main(int argc,char **argv)
{
  //params
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input mesh.")
    ("segmentationMapName,o", po::value<std::string>()->default_value("output"), "prefix of the segmentation map");
  bool parseOK=true;
  po::variables_map vm;
  try{
      po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
      trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
      parseOK=false;
  }
  //the name if segmentation map
  std::string segmentationMap = vm["segmentationMapName"].as<std::string>();
  //vector of id defects
  std::vector<unsigned int> idOfDefect;
  //the discretisation map
  std::vector<std::vector<std::vector<unsigned int>>> discretisation;
  //number of row to jump
  int rowCroppedBot=0;
  int rowCroppedTop=0;
  //read from discretisation.txt

  IOHelper::readDiscretisationFromFile("discretisation.txt",discretisation,rowCroppedBot,rowCroppedTop);

  //the segmentation image
  typedef ImageContainerBySTLVector < Z2i::Domain, unsigned char> Image;
  Image image2D = GenericReader<Image>::import(segmentationMap+"SEGTRESH.pgm" );

  //cols = width || rows = height
  int cols=image2D.domain().upperBound()[0];
  int rows=image2D.domain().upperBound()[1];

  //loop on segmentation map and fill a vector of defects indices
  int currentIntensity;
  std::vector<unsigned int>  currentPointsInPixels;
  unsigned int IndP;
  for (int i=0; i < cols; ++i){
    for (int j=0; j < rows; ++j){
      currentIntensity=image2D(Z2i::Point(i,j));
      if(currentIntensity>0){
        currentPointsInPixels=discretisation.at(j+rowCroppedBot).at(i);
        for(std::vector<unsigned int>::iterator it = std::begin(currentPointsInPixels); it != std::end(currentPointsInPixels); ++it) {
          IndP=*it;
          idOfDefect.push_back(IndP);
        }
      }
    }
  }
  //Read mesh file
  DGtal::Mesh<Z3i::RealPoint> mesh(true);
  std::string inputMeshName = vm["input"].as<std::string>();
  MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, mesh, false);
  //all the point in mesh
  std::vector<Z3i::RealPoint> pointCloud(mesh.nbVertex());
  //vector of boolean to match
  std::vector<bool> defectFlags(pointCloud.size(), false);
  for(unsigned int i = 0; i< idOfDefect.size(); i++){
      defectFlags[idOfDefect.at(i)] = true;
  }
  std::vector<unsigned int> facesToDelete;
  //color defect mesh
  for (unsigned int i = 0; i < mesh.nbFaces(); i++){
      Face aFace = mesh.getFace(i);
      unsigned int c = 0;
      for (unsigned int k = 0; k < aFace.size(); k++){
          //trace.info()<<aFace.at(k)<<std::endl;
          if(defectFlags.at(aFace.at(k))){
              c++;
          }
      }
      if(c >=  aFace.size()){
          mesh.setFaceColor(i, DGtal::Color::Green);
          facesToDelete.push_back(i);
      }
  }
  std::string defectFile = segmentationMap + "defect.off";
  std::string defectIdFile = segmentationMap + "-defect.id";
  IOHelper::export2OFF(mesh,defectFile);
  IOHelper::export2Text(idOfDefect,defectIdFile);

  return 0;
}
