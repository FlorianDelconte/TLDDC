#include "UnrolledMap.h"
#include <iostream>
#include <algorithm>
#include "CylindricalPoint.h"
#include "SegmentationAbstract.h"

//#include <opencv2/opencv.hpp>
#include <chrono>
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ArrayImageAdapter.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
using namespace DGtal;

void
UnrolledMap::computeDicretisation(){
    CylindricalPointOrder heightOrder;
    auto minMaxHeight = std::minmax_element(CPoints.begin(), CPoints.end(), heightOrder);
    double minHeight = (*minMaxHeight.first).height;
    double maxHeight = (*minMaxHeight.second).height;
    height_div=roundf(maxHeight-minHeight);
    //compute angle discretisation
    double meanRadius=0.;
    CylindricalPoint mpCurrent;
    for(unsigned int i = 0; i < CPoints.size(); i++){
        mpCurrent=CPoints.at(i);
        meanRadius+=mpCurrent.radius;
    }
    meanRadius/=CPoints.size();
    angle_div=roundf(2*M_PI*meanRadius);
    //compute min and max angle
    CylindricalPointOrderAngle angleOrder;
    auto minMaxAngle = std::minmax_element(CPoints.begin(), CPoints.end(), angleOrder);
    double minAngle = (*minMaxAngle.first).angle;
    double maxAngle = (*minMaxAngle.second).angle;

    //preallocate size for unrolled map
    unrolled_surface.resize(height_div);
    for (int i = 0; i < height_div; ++i){
        unrolled_surface[i].resize(angle_div);
    }
    //fill unrolled map
    trace.info()<<"Compute discretisation..."<<std::endl;
    int posAngle, posHeight;
    for(unsigned int i = 0; i < CPoints.size(); i++){
        mpCurrent=CPoints.at(i);
        //change range [minAngle,maxAngle] to [0,angle_div-1]
        posAngle=roundf((((angle_div-1)/(maxAngle-minAngle))*(mpCurrent.angle-(maxAngle)))+(angle_div-1));
        //change range [minHeight,maxHeight] to [0,height_div-1]
        posHeight=roundf((((height_div-1)/(maxHeight-minHeight))*(mpCurrent.height-maxHeight))+(height_div-1));
        //add index point to the unrolled_surface
        unrolled_surface[posHeight][posAngle].push_back(i);
    }
      trace.info()<<"Discretisation : [ "<<height_div<<" ; "<<angle_div<<" ]"<<std::endl;
}

/*TODO : change openCV*/
/*void
UnrolledMap::makeGroundTruthImage(std::vector<int> gtId,std::string outputFileName){
    trace.info()<<"create GT file..."<<std::endl;

    //cv::Mat gtImage(height_div,angle_div,CV_8UC1,cv::Scalar(0));
    //std::cout<<gtImage.size()<<std::endl;
    //Image2dGrayScale GT_image=Image2dGrayScale(  Z2i::Domain(Z2i::Point(0,0), Z2i::Point(angle_div,height_div)));
    //min and max angle
    CylindricalPointOrderAngle angleOrder;
    auto minMaxAngle = std::minmax_element(CPoints.begin(), CPoints.end(), angleOrder);
    double minAngle = (*minMaxAngle.first).angle;
    double maxAngle = (*minMaxAngle.second).angle;

    //max and min height
    CylindricalPointOrder heightOrder;
    auto minMaxHeight = std::minmax_element(CPoints.begin(), CPoints.end(), heightOrder);
    double minHeight = (*minMaxHeight.first).height;
    double maxHeight = (*minMaxHeight.second).height;

    //make gtImage
    int posAngle, posHeight, id;
    CylindricalPoint mpCurrent;
    for(unsigned int h = 0; h < gtId.size(); h++){
        id=gtId.at(h);
        mpCurrent=CPoints.at(id);
        //change range [minAngle,maxAngle] to [0,angle_div-1]
        posAngle=roundf((((angle_div-1)/(maxAngle-minAngle))*(mpCurrent.angle-(maxAngle)))+(angle_div-1));
        //change range [minHeight,maxHeight] to [0,height_div-1]
        posHeight=roundf((((height_div-1)/(maxHeight-minHeight))*(mpCurrent.height-maxHeight))+(height_div-1));
        //GT_image.setValue( Z2i::Point(posAngle,posHeight),255);
        //gtImage.at<uchar>(posHeight, posAngle) = 255;
    }

    //cv::Rect myROI(0, maxIndTop, angle_div, (minIndBot-maxIndTop)+1);
    //gtImage = gtImage(myROI);
    //std::cout<<gtImage.size()<<std::endl;
    //std::cout<<gtImage.size()<<std::endl;
    //crop top and bot
    Z2i::Domain subDomain(Z2i::Point(0,maxIndTop), Z2i::Point(angle_div-1,minIndBot));
    auto subGt_image = makeArrayImageAdapterFromImage( GT_image, subDomain );
    Image2dGrayScale croppedImage=Image2dGrayScale(subDomain);
    std::copy( subGt_image.begin(), subGt_image.end(), croppedImage.begin() );*/
    //Aply closing and opening operation
    /*int close_size = 1;
    int open_size = 1;
    cv::Mat elementClose = cv::getStructuringElement( cv::MORPH_RECT,cv::Size( 2*close_size + 1, 2*close_size+1 ),cv::Point( close_size, close_size ) );
    cv::Mat elementOpen= cv::getStructuringElement( cv::MORPH_RECT,cv::Size( 2*open_size + 1, 2*open_size+1 ),cv::Point( open_size, open_size ) );
    cv::morphologyEx( gtImage, gtImage, cv::MORPH_CLOSE , elementClose );
    cv::morphologyEx( gtImage, gtImage, cv::MORPH_OPEN, elementOpen );
    cv::flip(gtImage,gtImage,0);
    cv::imwrite( outputFileName+"_GT.pgm",gtImage);

}*/


bool
UnrolledMap::detectCellsIn(unsigned int i, unsigned int j){
    //A cell is 'in' if we can't reach the top image or the bot image with empty cellsgrayscalemap.at<uchar>(i, j) = grayscaleValue;
    bool CellsIn=true;

    std::vector<unsigned int > tempUp = unrolled_surface[i][j];
    std::vector<unsigned int > tempDown = unrolled_surface[i][j];

    //ind for reach the top image
    unsigned int iUp=i;
    //ind for reach the bot image
    unsigned int iDown=i;
    //decrement ind while cells is empty
    while(tempUp.empty() && (iUp>0)){
        iUp-=1;
        tempUp=unrolled_surface[iUp][j];
    }
    //increment ind while cells is empty
    while(tempDown.empty() && (iDown<height_div-1)){
        iDown+=1;
        tempDown=unrolled_surface[iDown][j];
    }
    //check reached top of bot
    if((iUp==0 && tempUp.empty()) || (iDown==height_div-1 && tempDown.empty())){
        CellsIn=false;
    }
    return CellsIn;
}

std::vector<unsigned int >
UnrolledMap::getIndPointsInLowerResolution(unsigned int i,unsigned int j,int dF){
    int pad=pow(2,dF);
    std::vector<unsigned int > outPutInd;
    int topLeftCornerHeight,topLeftCornerTheta;
    topLeftCornerHeight=(i/pad)*pad;
    topLeftCornerTheta=(j/pad)*pad;
    //loop on cells (of unrolledSurace) in region containing (i,j)
    for( int k = topLeftCornerHeight; k < (topLeftCornerHeight+pad); k++){
        for( int l = topLeftCornerTheta; l < (topLeftCornerTheta+pad); l++){
            //check if cells of region is in unrolled surface -> no segmentation fault
            if((k < height_div)&&(l < angle_div)){
                //unlarge  vector size
                outPutInd.reserve(outPutInd.size() + unrolled_surface[k][l].size());
                //concat vector
                outPutInd.insert(outPutInd.end(), unrolled_surface[k][l].begin(),  unrolled_surface[k][l].end());
            }
        }
    }
    return outPutInd;
}



std::vector<unsigned int>
UnrolledMap::getPointsUnrolled_surface(unsigned int i,unsigned int j){

    return unrolled_surface[i+maxIndTop][j];
}

double
UnrolledMap::sumReliefRepresentation(unsigned int i, unsigned int j,int dF){
    assert(dF!=0);
    double sumRelief=0.;
    std::vector<unsigned int> v = getIndPointsInLowerResolution(i,j,dF);
    if(!v.empty()){
        unsigned int IndP;
        for(std::vector<unsigned int>::iterator it = std::begin(v); it != std::end(v); ++it) {
            IndP=*it;
            sumRelief+=reliefRepresentation.at(IndP);
        }
    }else{
        sumRelief=-1;
    }
    return sumRelief;
}



double
UnrolledMap::meansReliefRepresentation(unsigned int i, unsigned int j,int dF){
    double moyenneRelief=0.;
    std::vector<unsigned int> v = getIndPointsInLowerResolution(i,j,dF);
    if(!v.empty()){
        unsigned int IndP;
        for(std::vector<unsigned int>::iterator it = std::begin(v); it != std::end(v); ++it) {
            IndP=*it;
            moyenneRelief+=reliefRepresentation.at(IndP);
        }
        moyenneRelief/=v.size();
    }else{
        moyenneRelief=-1;
    }
    return moyenneRelief;
}



double
UnrolledMap::maxReliefRepresentation(unsigned int i, unsigned int j,int dF){
    double maxReliefResolution=INT_MIN;
    double currentRelief;
    std::vector<unsigned int> v = getIndPointsInLowerResolution(i,j,dF);
    if(!v.empty()){
        unsigned int IndP;
        for(std::vector<unsigned int>::iterator it = std::begin(v); it != std::end(v); ++it) {
            IndP=*it;
            currentRelief=reliefRepresentation.at(IndP);
            if(currentRelief>maxReliefResolution){
                maxReliefResolution=currentRelief;
            }
        }
    }else{
        maxReliefResolution=-1;
    }
    return maxReliefResolution;
}

double
UnrolledMap::medianReliefRepresentation(unsigned int i, unsigned int j,int dF) {
    double medianRelief;
    std::vector<unsigned int> v = getIndPointsInLowerResolution(i,j,dF);
    unsigned int size = v.size();
    if(size!=0){
        std::sort(v.begin(), v.end());

        if (size % 2 == 0){
            double m1=reliefRepresentation.at(v.at(size / 2 - 1));
            double m2=reliefRepresentation.at(v.at(size / 2));
            medianRelief=(m1+m2)/2;
        }else{
            medianRelief=reliefRepresentation.at(v.at(size / 2));
        }
    }else{
        medianRelief=-1;
    }
    return medianRelief;
}

Image2dGrayScale
UnrolledMap::toGrayscaleImageMinMax(){
    int grayscaleValue;
    double reliefValue;
    /*DGtal*/
    Image2dGrayScale grayScaleReliefImage=Image2dGrayScale(reliefImage.domain());
    //SEARCH MIN MAX IN RELIEFIMAGE
    double minDG, maxDG;
    minDG=*min_element(reliefImage.range().begin(), reliefImage.range().end());
    maxDG=*max_element(reliefImage.range().begin(), reliefImage.range().end());
    trace.info()<<"min relief image dgtal :"<<minDG<<std::endl;
    trace.info()<<"max relief image dgtal:"<<maxDG<<std::endl;
    GrayscaleColorMap<float> grayShade(minDG,maxDG);

    for ( auto point : reliefImage.domain() ){
      reliefValue=reliefImage(point);
      grayscaleValue=((reliefValue-minDG)/(maxDG-minDG))*255;
      grayScaleReliefImage.setValue(point,grayscaleValue);
    }

    return grayScaleReliefImage;
}

Image2dGrayScale
UnrolledMap::toGrayscaleImageFixed(int intensityPerCm, double reliefValueforZero){
    float pad=1./intensityPerCm;
    double minDG=reliefValueforZero;
    double maxDG=(255*pad)+reliefValueforZero;
    int grayscaleValue;
    double reliefValue;
    //CREATE GRAYSCALE IMAGE
    Image2dGrayScale grayScaleReliefImage=Image2dGrayScale(reliefImage.domain());
    std::cout<<"GRAY : "<<reliefImage.domain()<<std::endl;
    //FILL THE GRAYSCALLE IMAGE
    for ( auto point : reliefImage.domain() ){
        reliefValue=reliefImage(point);
        //trace.info()<<point<<std::endl;
        grayscaleValue=((reliefValue-minDG)/(maxDG-minDG))*255;
        if(grayscaleValue<0){
            grayscaleValue=0;
        }
        if(grayscaleValue>255){
            grayscaleValue=255;
        }
        grayScaleReliefImage.setValue(point,grayscaleValue);
    }
    return grayScaleReliefImage;
}


void
UnrolledMap::computeNormalizedImage(int dF) {
    trace.info()<<"Compute normalized image with decrease factor : 1 / "<<pow(2,dF)<<" ..."<<std::endl;
    //resolution of relief imagev
    int pad=pow(2,dF);
    unsigned int resX=(angle_div-1)/pad;
    unsigned int resY=(height_div-1)/pad;
    trace.info()<<"resolution : Y = "<<resY<<" X = "<<resX<<std::endl;
    //resize reliefImageDGtal
    Z2i::Domain domain(Z2i::Point(0,0), Z2i::Point(resX,resY));
    reliefImage = Image2dNormalized(domain);
    //init normalized map with the new resolution
    double relief=0.;
    int XNormMap,YNormMap;
    //loop on the top left corner of all new cells
    for(unsigned int i = 0; i < height_div-1; i+=pad){
        for(unsigned int j = 0; j < angle_div-1; j+=pad){
            //check if the cells is in the mesh -> to avoid some noise in the image
            if(detectCellsIn(i,j)){
                YNormMap =i/pad;
                XNormMap=j/pad;
                //get radius in dF resolution
                relief=maxReliefRepresentation(i,j,pad);
                //Radius = -1 when cells is empty
                if(relief!=-1){
                  //reliefImage.at<float>(YNormMap, XNormMap) = relief;
                  reliefImage.setValue(Z2i::Point(XNormMap,YNormMap),relief);
                }
            }

        }
    }
    //can't crop top and bot here
}


void
UnrolledMap::computeNormalizedImageMultiScale(){
    trace.info()<<"Compute normalized image in multi scale ..."<<std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    //resize reliefImageDGtal
    Z2i::Domain domain(Z2i::Point(0,0), Z2i::Point(angle_div,height_div));
    reliefImage = Image2dNormalized(domain);
    int dF;
    int maxDecreaseHit=0;
    double relief;
    double normalizedRelief;
    //counter
    int counterdf0=0;
    int counterdf1=0;
    int counterdf2=0;
    int counterdf3=0;
    int counterdf4=0;
    int counter_out=0;
    int counter_un=0;
    //loop on all cells of unrolled surface
    for(unsigned int i = 0; i < height_div; i++){
        for(unsigned int j = 0; j < angle_div; j++){
            if(detectCellsIn(i,j)){
                relief=maxReliefRepresentation(i,j,0);
                dF=1;


                //while this vector is empty then find a little resolution where the corresponding i,j cells is not empty
                while(relief==-1 && dF<maxDecreaseFactor){
                    relief=maxReliefRepresentation(i,j,dF);
                    //std::cout<<dF<<std::endl;

                    //keep the max decrease resolution
                    if(dF>=maxDecreaseHit){
                        maxDecreaseHit=dF;
                    }
                    //decrease resolution (1/decreaseHit)
                    dF+=1;
                }
                reliefImage.setValue(Z2i::Point(j,i),relief);
            }else{
              counter_out+=1;
            }
        }
    }
    /*auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int nbPixels=(angle_div)*(height_div);
    trace.info()<<"decrease paramater down to : "<<maxDecreaseHit<<std::endl;
    trace.info()<<"RATIO (for : "<<nbPixels<<" nb picels)"<<std::endl;
    trace.info()<<"at df=0 completion : "<<double(counterdf0)/nbPixels<< "nb pixel search in df=0 "<<counterdf0<<std::endl;
    trace.info()<<"at df=1 completion : "<<double(counterdf1)/nbPixels<< "nb pixel search in df=1 "<<counterdf1<<std::endl;
    trace.info()<<"at df=2 completion : "<<double(counterdf2)/nbPixels<< "nb pixel search in df=2 "<<counterdf2<<std::endl;
    trace.info()<<"at df=3 completion : "<<double(counterdf3)/nbPixels<< "nb pixel search in df=3 "<<counterdf3<<std::endl;
    trace.info()<<"at df=4 completion : "<<double(counterdf4)/nbPixels<< "nb pixel search in df=4 "<<counterdf4<<std::endl;

    trace.info()<<counter_out<< " pixels not traited so : "<<double(counter_out)/nbPixels<<" ratio  "<<std::endl;
    trace.info()<<"at beginin they are  : "<<nbPixels-counter_un<<" pixels set"<<std::endl;
    trace.info()<<"end multi scale research ... Duration : "<<duration.count()<<" microseconds"<<std::endl;*/
    cropTopBotImage();
}

void
UnrolledMap::cropTopBotImage(){
    trace.info()<<"start crop top and bot ..."<<std::endl;
    int colsDG=reliefImage.domain().upperBound()[0];
    int rowsDG=reliefImage.domain().upperBound()[1];

    double currentPixelValue;
    unsigned int x,y;
    //Search the max ind point (different zeros) from top normalized image
    unsigned int maxIT_DGtal=0;
    for(x = 0; x < colsDG; x++){
        currentPixelValue=reliefImage(Z2i::Point(x,0));
        y=0;
        while (currentPixelValue==0) {
            y=y+1;
            currentPixelValue=reliefImage(Z2i::Point(x,y));
        }
        if(y>maxIT_DGtal){
            maxIT_DGtal=y;
        }
    }
    unsigned int minIB_DGtal=height_div-1;
    for(x = 0; x < colsDG; x++){
      currentPixelValue=reliefImage(Z2i::Point(x,height_div-1));
      y=height_div-1;
      while (currentPixelValue==0) {
          y=y-1;
          currentPixelValue=reliefImage(Z2i::Point(x,y));
      }
      if(y<minIB_DGtal){
          minIB_DGtal=y;
      }
    }
    typedef ConstImageAdapter<Image2dNormalized, Z2i::Domain, functors::Identity, Image2dNormalized::Value, functors::Identity > ConstImageAdapterForSubImage;
    functors::Identity df;


    //Crop an image
    Z2i::Domain subDomain(Z2i::Point(0,maxIT_DGtal), Z2i::Point(angle_div-1,minIB_DGtal-1));
    auto subImageSTL = makeArrayImageAdapterFromImage( reliefImage, subDomain );
    Image2dNormalized croppedImage=Image2dNormalized(subDomain);
    std::copy( subImageSTL.begin(), subImageSTL.end(), croppedImage.begin() );
    //update attribut maxIndTop and minIndBot
    maxIndTop=maxIT_DGtal;
    minIndBot=minIB_DGtal;
    trace.info()<< "max IndTop :  "<<maxIndTop<< "min indBot : "<<minIndBot<< std::endl;
    reliefImage=croppedImage;

}

void
UnrolledMap::computeRGBImage(){
    double reliefValue;
    Color reliefColor;
    //first : compute in grayScale
    Image2dGrayScale reliefGray = toGrayscaleImageFixed(intensity_per_cm,zero_level_intensity);
    //second : convert grayscale to rgb
    imageRGB imagergb=imageRGB(reliefGray.domain());

    float min=*min_element(reliefGray.range().begin(), reliefGray.range().end());
    float max=*max_element(reliefGray.range().begin(), reliefGray.range().end());
    GradientColorMap<float,CMAP_COPPER> gradient( min, max);
    for ( auto point : reliefGray.domain() ){
        reliefValue=reliefGray(point);
        //trace.info()<<point<<std::endl;
        reliefColor=gradient(reliefValue);

        imagergb.setValue(point,reliefColor);
    }

    reliefImageRGB=imagergb;
}
void
UnrolledMap::computeGRAYImage(){

    if(intensity_per_cm==-1){
      trace.info()<<"make grayscale image with [min ;max]"<< std::endl;
      reliefImageGrayScale=toGrayscaleImageMinMax();
    }else{
      trace.info()<<"make grayscale with zero = "<<zero_level_intensity<< " and pad = "<<intensity_per_cm<< std::endl;
      reliefImageGrayScale=toGrayscaleImageFixed(intensity_per_cm,zero_level_intensity);
    }

}
/**GETTERS**/
Image2dNormalized
UnrolledMap::getNormalizedImage(){
    return reliefImage;
}

Image2dGrayScale
UnrolledMap::getReliefImageGrayScale(){
    return reliefImageGrayScale;
}


imageRGB
UnrolledMap::getReliefImageRGB(){
    return reliefImageRGB;
}

CylindricalPoint
UnrolledMap::getCPoint(unsigned int i){
    return CPoints.at(i);
}
std::vector<std::vector<std::vector<unsigned int>>>
UnrolledMap::getDiscretisation(){
  return unrolled_surface;
}
int
UnrolledMap::getRowCroppedBot(){
  return maxIndTop;
}
int
UnrolledMap::getRowCroppedTop(){
  return minIndBot;
}
