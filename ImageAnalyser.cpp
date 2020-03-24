#include "ImageAnalyser.h"
#include "UnrolledMap.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/opencv.hpp>

#include <boost/tuple/tuple.hpp>
//#include <boost/array.hpp>
//#include <boost/range/adaptor/transformed.hpp>
//#include <boost/range/irange.hpp>
//#include <boost/bind.hpp>
//namespace plt = matplotlibcpp;
//using namespace cv;


void 
ImageAnalyser::onMouse( int event, int x, int y, int flags, void* param)
{
    if( event != cv::EVENT_LBUTTONDOWN )
        return;
    cv::Point seed = cv::Point(x,y);
    ImageAnalyser *p_ianalyser=(ImageAnalyser*) param;
    std::vector<unsigned int > indInCells=p_ianalyser->unrolled_map.getIndPointsInLowerResolution(seed.y,seed.x,1);
    std::cout<<"nb points in cells : "<<indInCells.size()<<std::endl;

    if(!indInCells.empty()){
        //delete old .dat and png
        system("rm ../clickedPatchOutput/*.dat");
        system("rm ../clickedPatchOutput/*.png");
        unsigned int IndP;
        std::ofstream dataFile;
        CylindricalPoint mpFound;
        unsigned int nbPatches=indInCells.size();
        std::string datafilename;
        //create new .dat files
        for (unsigned int idPatches = 0; idPatches < nbPatches; ++idPatches){
            datafilename="../clickedPatchOutput/data"+std::to_string(idPatches)+".dat";
            dataFile.open(datafilename);
            dataFile<<"#X\tY\n";
            for(std::vector<unsigned int>::iterator it = std::begin(p_ianalyser->ind_Patches[indInCells[idPatches]]); it != std::end(p_ianalyser->ind_Patches[indInCells[idPatches]]); ++it) {
                IndP=*it;
                mpFound = p_ianalyser->CPoints.at(IndP);
                dataFile << std::fixed << std::setprecision(8) <<mpFound.height<<"\t"<<mpFound.radius <<"\n";
                
            }
            dataFile.close();
        }
        //create png with gnuplot
        FILE *fp;
        fp = popen("cd ../clickedPatchOutput/; gnuplot","w");
        fprintf(fp,"set terminal png\n");
        fprintf(fp,"set xlabel \"x\"\n");
        fprintf(fp,"set ylabel \"y\"\n");
        for (unsigned int idPatches = 0; idPatches < nbPatches; ++idPatches){
            fprintf(fp,"set output \"patch%d.png\"\n",idPatches);
            fprintf(fp, "plot \"data%d.dat\" using 1:2 title \"patch%d\"\n",idPatches,idPatches);
        }
        pclose(fp);
    }
    
    
}

void
ImageAnalyser::analyse(){
    cv::Mat image=unrolled_map.getRgbImage();
    //std::cout<<ind_Patches.at(100).size()<<std::endl;
    cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );// Create a window for display.
    cv::setMouseCallback("Display window", &ImageAnalyser::onMouse,this);
    cv::imshow( "Display window", image );                   // Show our image inside it.
    while (char(cv::waitKey(1) != 'q')){
        int c;
        c = cv::waitKey(1);
        if ((char)c == 'u'){
            cv::pyrUp(image, image, cv::Size(image.cols * 2, image.rows * 2));
            std::cout << " ** Zoom In: Image x 2"<<std::endl;
        }else if ((char)c == 'd'){       
            cv::pyrDown(image, image, cv::Size(image.cols / 2, image.rows / 2));
            std::cout<<" ** Zoom Out: Image / 2"<<std::endl;
        }
        cv::imshow("Display window", image);
    }
}


