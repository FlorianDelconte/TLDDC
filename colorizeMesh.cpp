#include <iostream>
#include <fstream>

#include <ctime>

#ifndef Q_MOC_RUN
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



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/shapes/Mesh.h"

#endif
#include "IOHelper.h"

using namespace DGtal;
namespace po = boost::program_options;




typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;
typedef ImageContainerBySTLVector<Z3i::Domain,  std::vector<Z3i::RealPoint> > ImagePointAssociation;



typedef DGtal::ConstImageAdapter<Image3D, Z2i::Domain, DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,
        unsigned int,  DGtal::functors::Identity >  ImageAdapterExtractor;



int
main(int argc,char **argv)

{
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("mesh1,i", po::value<std::string>(), "input file name of mesh vertex given as OFF format.")
        ("groundtrue,r", po::value<std::string>(), "face ids")
        ("result,t", po::value<std::string>(), "face ids")
        ("output,o", po::value<std::string>(), "output file name of mesh vertex given as OFF format.");

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }catch(const std::exception& ex){
        trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
        parseOK=false;
    }
    po::notify(vm);
    if(vm.count("help")||argc<=1|| !parseOK)
    {
        trace.info()<< "Options: "<<std::endl
            << general_opt << "\n";
        return 0;
    }
    std::string meshName = vm["mesh1"].as<std::string>();
    std::string groundtrueIdsFile = vm["groundtrue"].as<std::string>();
    std::string resultIdsFile = vm["result"].as<std::string>();

    DGtal::Mesh<Z3i::RealPoint> mesh1(true);
    DGtal::Mesh<Z3i::RealPoint> mesh2(true);
    MeshReader<Z3i::RealPoint>::importOFFFile(meshName, mesh1, false);
    MeshReader<Z3i::RealPoint>::importOFFFile(meshName, mesh2, false);

    std::vector<int> groundtrueIds;
    IOHelper::readIntsFromFile(groundtrueIdsFile, groundtrueIds);
    std::vector<int> resultIds;
    IOHelper::readIntsFromFile(resultIdsFile, resultIds);
    trace.info()<< "Import mesh 1 ok, new size: "<< mesh1.nbVertex() <<  std::endl;
    std::vector<bool> groundtrueFlags(mesh1.nbFaces(),false);
    std::vector<bool> resultFlags(mesh1.nbFaces(),false);

    for(unsigned int i = 0; i < groundtrueIds.size(); i++){
        unsigned int index = groundtrueIds.at(i);
        groundtrueFlags[index] = true;
    }
    for(unsigned int i = 0; i < resultIds.size(); i++){
        unsigned int index = resultIds.at(i);
        resultFlags[index] = true;
    }
    //DGtal::Mesh<Z3i::RealPoint> outMesh = mesh1;
    for(unsigned int i =0; i< mesh1.nbFaces(); i++){
        DGtal::Color color = DGtal::Color(150, 150, 150, 255);
        if(groundtrueFlags.at(i) && resultFlags.at(i)){
            color = DGtal::Color::Yellow;
        }else if(groundtrueFlags.at(i)){
            color = DGtal::Color::Red;
        }else if(resultFlags.at(i)){
            color = DGtal::Color::Green;
        }
        //mesh1.setFaceColor(i, DGtal::Color(200, 200 ,200, 255));
        mesh1.setFaceColor(i, color);
    }
    trace.info()<< "Concat mesh 2 into mesh 1 ok, new size: "<< mesh1.nbVertex() <<  std::endl;
    std::string outName = vm["output"].as<std::string>();
    IOHelper::export2OFF(mesh1, outName);

    //std::ofstream outMesh(outName.c_str());
    //std::vector<unsigned int> fids;
    //for(unsigned int i = 0; i < groundtrueFlags.size(); i++){
    //    if(!groundtrueFlags[i]){
    //        fids.push_back(i);
    //    }
    //}
    //mesh1.removeFaces(fids);
    //IOHelper::export2OFF(mesh1, "test"+outName);
    //export yellow part
}
