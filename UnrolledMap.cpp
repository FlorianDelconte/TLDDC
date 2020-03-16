#include "UnrolledMap.h"
#include <iostream>
#include <algorithm>
#include "CylindricalPoint.h"
#include "SegmentationAbstract.h"
#include <opencv2/opencv.hpp>


UnrolledMap::UnrolledMap(std::vector<CylindricalPoint> CylindricalPoints,std::vector<double> rRepresenation){
    trace.info()<<"Construct unrolled map..."<<std::endl;
    //init attribut
    CPoints=CylindricalPoints;
    reliefRepresentation=rRepresenation;
    //compute min and max on DDistance for normalisation
    maxRelief=INT_MIN;
    minRelief=INT_MAX;
    for(unsigned int i = 0; i < reliefRepresentation.size(); i++){
        if(reliefRepresentation[i]>maxRelief){
            maxRelief=reliefRepresentation[i];
        }
        if(reliefRepresentation[i]<minRelief){
            minRelief=reliefRepresentation[i];
        }
    }
    //compute Height discretisation
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
    trace.info()<<"Discretisation : [ "<<height_div<<" ; "<<angle_div<<" ]"<<std::endl;
    //preallocate size for unrolled map
    unrolled_surface.resize(height_div);
    for (int i = 0; i < height_div; ++i){
        unrolled_surface[i].resize(angle_div);
    }
    //fill unrolled map
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
}

bool 
UnrolledMap::detectCellsIn(unsigned int i, unsigned int j){
    //A cell is 'in' if we can't reach the top image or the bot image with empty cells
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
    std::vector<unsigned int > outPutInd;
    int topLeftCornerHeight,topLeftCornerTheta;
    topLeftCornerHeight=(i/dF)*dF;
    topLeftCornerTheta=(j/dF)*dF;
    //loop on cells (of unrolledSurace) in region containing (i,j)
    for( int k = topLeftCornerHeight; k < (topLeftCornerHeight+dF); k++){
        for( int l = topLeftCornerTheta; l < (topLeftCornerTheta+dF); l++){
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

double
UnrolledMap::meansReliefRepresentation(unsigned int i, unsigned int j,int dF){
    assert(dF!=0);
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
    assert(dF!=0);
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
    assert(dF!=0);
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

double 
UnrolledMap::normalizeReliefRepresentation(double value){
  return ((1/(maxRelief-minRelief))*(value-maxRelief))+1;
}

cv::Mat
UnrolledMap::computeNormalizedImage(int dF) {
    trace.info()<<"start compute normalized image with decrease factor : 1 / "<<dF<<" ..."<<std::endl;
    //resolution of relief image
    unsigned int resX=angle_div/dF;
    unsigned int resY=height_div/dF;
    trace.info()<<"resolution : Y = "<<resY<<" X = "<<resX<<std::endl;
    //init normalized map with the new resolution
    cv::Mat normalizedMap=cv::Mat::zeros(resY+1,resX+1,CV_32FC4);
    double relief=0.;
    double normalizedRelief;
    int XNormMap,YNormMap;
    //loop on the top left corner of all new cells
    for(unsigned int i = 0; i < height_div; i+=dF){
        for(unsigned int j = 0; j < angle_div; j+=dF){
            //check if the cells is in the mesh -> to avoid some noise in the image
            if(detectCellsIn(i,j)){
                //[0;height_div] to [0;resY]
                YNormMap =i/dF;
                //[0;angle_div] to [0;resX]
                XNormMap=j/dF;
                //get radius in dF resolution
                relief=medianReliefRepresentation(i,j,dF);
                //Radius = -1 when cells is empty
                if(relief!=-1){
                    normalizedRelief=normalizeReliefRepresentation(relief);
                    normalizedMap.at<double>(YNormMap, XNormMap) = normalizedRelief;
                }
            }
        }
    }
    return normalizedMap;
}


cv::Mat
UnrolledMap::computeNormalizedImageMultiScale(){
    trace.info()<<"start compute normalized image in multi scale ..."<<std::endl;
    int dF;
    int maxDecreaseHit=0;
    double relief;
    cv::Mat normalizedMap=cv::Mat::zeros(height_div,angle_div,CV_32FC4);
    double normalizedRelief;
    //loop on all cells of unrolled surface
    for(unsigned int i = 0; i < height_div; i++){
        for(unsigned int j = 0; j < angle_div; j++){
            if(detectCellsIn(i,j)){
                dF=2;
                relief=-1;
                //while this ector is empty find a little resolution where the corresponding i,j cells is not empty
                while(relief==-1 && dF<32){
                    relief=medianReliefRepresentation(i,j,dF);
                    //decrease resolution (1/decreaseHit)
                    dF*=2;
                    //keep the max decrease resolution -> not used for the moment
                    if(dF>=maxDecreaseHit){
                        maxDecreaseHit=dF;
                    }
                }
                normalizedRelief=normalizeReliefRepresentation(relief);
                normalizedMap.at<double>(i, j) = normalizedRelief;   
            }
        }
    }
    return normalizedMap;
}