#include <iostream>
#include <fstream>
#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include "DefectSegmentationUnroll.h"
#include "IOHelper.h"
#include "Centerline/Centerline.h"
#include "Centerline/CenterlineHelper.h"

//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>


using namespace DGtal;
namespace po = boost::program_options;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;



int
main(int argc,char **argv)
{
  /**************************/
  /**Gestion des paramètres**/
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "input mesh.")
        ("accRadius,r", po::value<double>()->default_value(200), "accumulation radius.")
        ("trackStep,s", po::value<double>()->default_value(20), "tracking step.")
        ("invertNormal,n", "invert normal to apply accumulation.")
        ("binWidth,b", po::value<double>()->default_value(0.01), "bin width used to compute threshold")
        ("patchWidth,a", po::value<double>()->default_value(25), "Arc length/ width of patch")
        ("patchHeight,e", po::value<int>()->default_value(100), "Height of patch")
        ("voxelSize", po::value<int>()->default_value(5), "Voxel size")
        ("decreaseFactor,d", po::value<int>()->default_value(4), "Max decrease factor for multi resolution search")
        ("grayscaleOrigin", po::value<int>()->default_value(-5), "relief value for 0 level in grayscale intensity")
        ("intensityPerCm", po::value<int>()->default_value(10), "number of grayscale intensity to represente 1cm of relief")
        ("output,o", po::value<std::string>()->default_value("output"), "output prefix: output-defect.off, output-def-faces-ids, ...");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }

    po::notify(vm);
    if(vm.count("help") || argc<=1 || !parseOK || !vm.count("input") || !vm.count("accRadius") || !vm.count("trackStep")){
        if(!vm.count("input")){
            trace.error()<<"the input mesh is required!"<<std::endl;
        }else if( !vm.count("accRadius") ){
            trace.error()<<"the accRadius is required!"<<std::endl;
        }else if( !vm.count("trackStep") ){
            trace.error()<<"the trackStep is required!"<<std::endl;
        }
        trace.info()<< "Segmentation log defects" <<std::endl << "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }

    int voxelSize = vm["voxelSize"].as<int>();
    assert(voxelSize > 0);

    DGtal::Mesh<Z3i::RealPoint> scaledMesh(true);

    double accRadius = vm["accRadius"].as<double>() / voxelSize;
    double trackStep = vm["trackStep"].as<double>() / voxelSize;
    bool invertNormal = vm.count("invertNormal");
    //trace.info()<< "invertNormal :::::::::" <<invertNormal<<std::endl;
    double binWidth = vm["binWidth"].as<double>();
    int maxDecreaseFactor=vm["decreaseFactor"].as<int>();
    int gs_origin=vm["grayscaleOrigin"].as<int>();
    int intensity_cm=vm["intensityPerCm"].as<int>();
    DGtal::Mesh<Z3i::RealPoint> oriMesh(true);
    std::string inputMeshName = vm["input"].as<std::string>();

    /**************************/


    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, scaledMesh, false);
    MeshReader<Z3i::RealPoint>::importOFFFile(inputMeshName, oriMesh, false);
    std::vector<Z3i::RealPoint> pointCloud(scaledMesh.nbVertex());

    //adjust by voxelSize
    for ( int i = 0; i < scaledMesh.nbVertex(); i++ ){
        Z3i::RealPoint &p = scaledMesh.getVertex(i);
        p /= voxelSize;
    }
    trace.info()<<"Cloud size : "<< pointCloud.size()<< std::endl;

    Centerline acc(scaledMesh, accRadius, trackStep, invertNormal);

    std::copy(oriMesh.vertexBegin(), oriMesh.vertexEnd(), pointCloud.begin());


    //@TODO:check input mesh and fiber here
    std::vector<Z3i::RealPoint> fiber = acc.compute();
    //unscale fiber for more accuracy Splines
    for(unsigned int i = 0; i < fiber.size(); i++){
        fiber[i] = fiber[i]*voxelSize;
    }

    std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> boudingBox = oriMesh.getBoundingBox();
    Z3i::RealPoint ptLow = boudingBox.first;
    Z3i::RealPoint ptUp = boudingBox.second;
    Z3i::Domain domain = Z3i::Domain(Z3i::Point((int) ptLow[0], (int) ptLow[1], (int) ptLow[2]),
            Z3i::Point((int) ptUp[0], (int) ptUp[1], (int) ptUp[2]));

    std::vector<Z3i::RealPoint> centerline = CenterlineHelper::getSmoothCenterlineBSplines(domain, fiber);
    trace.info()<<"centerline size : "<< fiber.size()<< std::endl;
    trace.info()<<"centerline smoothed size : "<< centerline.size()<< std::endl;
    // Uncomment to test interpolated centerline
    //write centerline
    Mesh<Z3i::RealPoint> transMesh = oriMesh;
    for(unsigned int i =0; i< transMesh.nbFaces(); i++){
        transMesh.setFaceColor(i, DGtal::Color(120, 120 ,120, 180));
    }
    Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, fiber, 1, 0.1, DGtal::Color::Blue);
    Mesh<Z3i::RealPoint>::createTubularMesh(transMesh, centerline, 1, 0.1, DGtal::Color::Red);

    IOHelper::export2OFF(transMesh, "centerline.off");


    double patchWidth = vm["patchWidth"].as<double>();
    int patchHeight = vm["patchHeight"].as<int>();
    std::string outputPrefix = vm["output"].as<std::string>();
    size_t lastindex = inputMeshName.find_last_of(".");

    std::string GtFileName = inputMeshName.substr(0, lastindex)+"-groundtruth-points.id";




    DefectSegmentationUnroll sa(pointCloud,centerline,patchWidth,patchHeight,binWidth);
    sa.init();
    sa.makeRM(outputPrefix,GtFileName, maxDecreaseFactor,gs_origin,intensity_cm);



    return 0;


}
