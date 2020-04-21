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
    std::cout<<"x:"<<seed.x<<"y: "<<seed.y<<std::endl;
    ImageAnalyser *p_ianalyser=(ImageAnalyser*) param;
    std::vector<unsigned int > indInCells=p_ianalyser->unrolled_map.getIndPointsInLowerResolution(seed.y,seed.x,1);
    std::cout<<"nb points in cells "<<indInCells.size()<<std::endl;
    
    if(!indInCells.empty()){
        //delete old .dat and png
        system("rm ../clickedPatchOutput/*.dat");
        system("rm ../clickedPatchOutput/*.png");
        unsigned int IndP,indLine;
        std::ofstream dataFile;
        CylindricalPoint mpFound;
        unsigned int nbPatches=indInCells.size();
        std::string datafilename;
        double a,b;
        //fill  patches.dat files
        for (unsigned int idPatches = 0; idPatches < nbPatches; ++idPatches){
           
            datafilename="../clickedPatchOutput/patche"+std::to_string(idPatches)+".dat";
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
        fprintf(fp,"set xlabel \"z (mm)\"\n");
        fprintf(fp,"set ylabel \"distance to the centerline (mm)\"\n");
        CylindricalPoint shiftedPoint;
        for (unsigned int idPatches = 0; idPatches < nbPatches; ++idPatches){
            IndP=indInCells[idPatches];
            mpFound = p_ianalyser->CPoints.at(IndP);
            //shiftedPoint=p_ianalyser->getShiftedPoint(p_ianalyser->ind_CoefsLines.at(IndP),p_ianalyser->ind_Patches.at(IndP));
            //std::cout<<"h:"<<mpFound.height<<"r: "<<mpFound.radius<<std::endl;
            fprintf(fp,"set output \"patch%d.png\"\n",idPatches);

            fprintf(fp,"a=%f\n", p_ianalyser->ind_CoefsLines.at(IndP).first);
            fprintf(fp,"b=%f\n", p_ianalyser->ind_CoefsLines.at(IndP).second);
            fprintf(fp,"f(x) = a*x + b\n");

            fprintf(fp, "plot \"patche%d.dat\" using 1:2 title \"points in the patch\",  ",idPatches,idPatches);
            fprintf(fp, "\"<echo '%f %f'\" pt 7 ps 2 title \"processing point\",",mpFound.height,mpFound.radius);
            //fprintf(fp, "\"<echo '%f %f'\" pt 7 ps 3 title \"farest poit\",",shiftedPoint.height,shiftedPoint.radius);
            fprintf(fp, "f(x) title 'fitted line' with lines linestyle 1 lc 4\n");
        }
        pclose(fp);
    }
    
    
}
CylindricalPoint
ImageAnalyser::getShiftedPoint(std::pair<double, double>& coef,std::vector<unsigned int> pointsInPatch){
    unsigned int shiftedPoint;
    CylindricalPoint mp;
    CylindricalPoint maxP;
    double A=coef.first;
    double B=-1;
    double C=coef.second;
    double distanceLine;
    //loop on points in current
    //std::cout<<" a : "<<A<<std::endl;
    //std::cout<<" b : "<<B<<std::endl;
    //std::cout<<" c : "<<C<<std::endl;
    
    double maxdistance=INT_MIN;
    for(std::vector<unsigned int>::iterator it = std::begin(pointsInPatch); it != std::end(pointsInPatch); ++it) {
        mp = CPoints.at(*it);
        //yp < f(xp)
        if(mp.radius - (coef.first*mp.height + coef.second) < 0 ){
            //std::cout<<" point at "<<mp.height<<" ; " << mp.radius<<"is below the fitted line"<<std::endl;
            distanceLine=abs((A*mp.height)+(B*mp.radius)+C)/sqrt((A*A)+(B*B));
            if(distanceLine>maxdistance){
                maxdistance=distanceLine;
                maxP=mp;
            }
        }
    }
    //std::cout<<" max distance is  : "<<maxdistance<<std::endl;
    //std::cout<<" attaint par le point : ["<<maxP.height<<" ; "<<maxP.radius<<" ] "<<std::endl;
    coef.second-=maxdistance;
    return maxP;
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


