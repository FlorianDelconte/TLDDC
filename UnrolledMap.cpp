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
    //init reliefImage
    reliefImage = cv::Mat::zeros(height_div,angle_div,CV_32F);
  
}

UnrolledMap::UnrolledMap(const UnrolledMap &um){
    minRelief = um.minRelief;
    maxRelief = um.maxRelief;
    unrolled_surface=um.unrolled_surface;
    reliefRepresentation=um.reliefRepresentation;
    CPoints=um.CPoints;
    height_div=um.height_div;
    angle_div=um.angle_div;
    reliefImage=um.reliefImage;
    image=um.image;
}
cv::Mat
UnrolledMap::makeGroundTruthImage(std::vector<int> gtId){
    trace.info()<<"create GT file..."<<std::endl;
    cv::Mat gtImage(height_div,angle_div,CV_8UC1,cv::Scalar(0));
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
        //trace.info()<<"max"<<height_div<<" "<<angle_div<<std::endl;
        //trace.info()<<"coord in maillage "<<mpCurrent.height<<" "<<mpCurrent.angle<<std::endl;
        //change range [minAngle,maxAngle] to [0,angle_div-1]
        posAngle=roundf((((angle_div-1)/(maxAngle-minAngle))*(mpCurrent.angle-(maxAngle)))+(angle_div-1));
        //change range [minHeight,maxHeight] to [0,height_div-1]
        posHeight=roundf((((height_div-1)/(maxHeight-minHeight))*(mpCurrent.height-maxHeight))+(height_div-1));
        //trace.info()<<"coord in image "<<posHeight<<" "<<posAngle<<std::endl;
        //add index point to the unrolled_surface
        gtImage.at<uchar>(posHeight, posAngle) = 255;
    }

    //crop bot and top
    cv::Rect myROI(0, maxIndTop, angle_div, minIndBot-maxIndTop);
    gtImage = gtImage(myROI);
    //Aply closing and opening operation
    int close_size = 2;
    int open_size = 2;
    cv::Mat elementClose = cv::getStructuringElement( cv::MORPH_RECT,cv::Size( 2*close_size + 1, 2*close_size+1 ),cv::Point( close_size, close_size ) );
    cv::Mat elementOpen= cv::getStructuringElement( cv::MORPH_RECT,cv::Size( 2*open_size + 1, 2*open_size+1 ),cv::Point( open_size, open_size ) );
    cv::morphologyEx( gtImage, gtImage, cv::MORPH_CLOSE , elementClose );
    cv::morphologyEx( gtImage, gtImage, cv::MORPH_OPEN, elementOpen );
    return gtImage;
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

cv::Mat
UnrolledMap::toGrayscaleImageMinMax(){
    int grayscaleValue;
    double reliefValue;
    int rows = reliefImage.rows;
    int cols = reliefImage.cols;
    cv::Mat grayscalemap(rows,cols,CV_8UC1,cv::Scalar(0));
    double min, max;
    cv::minMaxLoc(reliefImage, &min, &max);
    trace.info()<<"min relief image:"<<min<<std::endl;
    trace.info()<<"max relief image:"<<max<<std::endl;
    for(unsigned int i = 0; i < rows; i++){
        for(unsigned int j = 0; j < cols; j++){
            reliefValue=reliefImage.at<float>(i, j);
            grayscaleValue=((reliefValue-min)/(max-min))*255;
            grayscalemap.at<uchar>(i, j) = grayscaleValue;
        }
    }
    double mingm, maxgm;
    cv::minMaxLoc(grayscalemap, &mingm, &maxgm);
    trace.info()<<"min gray map:"<<mingm<<std::endl;
    trace.info()<<"max gray map:"<<maxgm<<std::endl;
    return grayscalemap;
}
double
UnrolledMap::ownNormalizeReliefRepresentation(double value){
    int minValue=0;
    int maxValue=25;
    if(value>=maxValue){
        trace.info()<<"true in min"<<std::endl;
        value=maxValue;
    }
    if(value<=minValue){
         trace.info()<<"true in max"<<std::endl;
        value=minValue;
    }
    return (value-minValue)/(maxValue-minValue);
}


void
UnrolledMap::computeNormalizedImage(int dF) {
    trace.info()<<"start compute normalized image with decrease factor : 1 / "<<dF<<" ..."<<std::endl;
    //resolution of relief image
    unsigned int resX=angle_div/dF;
    unsigned int resY=height_div/dF;
    trace.info()<<"resolution : Y = "<<resY<<" X = "<<resX<<std::endl;
    //init normalized map with the new resolution
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
                relief=meansReliefRepresentation(i,j,dF);
                //Radius = -1 when cells is empty
                if(relief!=-1){
                    reliefImage.at<float>(YNormMap, XNormMap) = normalizedRelief;       
                }
            }

        }
    }
    cropTopBotImage();
    //normalizeImage();
}


void
UnrolledMap::computeNormalizedImageMultiScale(){
    trace.info()<<"start compute normalized image in multi scale ..."<<std::endl;
    int dF;
    int maxDecreaseHit=0;
    double relief;
    double normalizedRelief;
    //loop on all cells of unrolled surface
    for(unsigned int i = 0; i < height_div; i++){
        for(unsigned int j = 0; j < angle_div; j++){
            if(detectCellsIn(i,j)){
                dF=2;
                relief=-1;
                //while this ector is empty find a little resolution where the corresponding i,j cells is not empty
                while(relief==-1 && dF<32){
                    relief=maxReliefRepresentation(i,j,dF);
                    //decrease resolution (1/decreaseHit)
                    dF*=2;
                    //keep the max decrease resolution -> not used for the moment
                    if(dF>=maxDecreaseHit){
                        maxDecreaseHit=dF;
                    }
                }
                reliefImage.at<float>(i, j) = relief;
            }
        }
    }
    cropTopBotImage();
    //normalizeImage();
   
    
}

//@TODO: Est ce qu'il ne vaut pas mieux crop la carte de relief ??
void
UnrolledMap::cropTopBotImage(){
    int rows = reliefImage.rows;
    int cols = reliefImage.cols;
    double currentPixelValue;
    unsigned int i,j;
    //Search the max ind point (different zeros) from top normalized image
    unsigned int maxIT=0;
    for(j = 0; j < cols; j++){
        //trace.info()<<"------"<<std::endl;
        currentPixelValue=reliefImage.at<float>(0, j);
        i=0;
        while (currentPixelValue==0) {
            i=i+1;
            currentPixelValue=reliefImage.at<float>(i, j);
            //trace.info()<<currentPixelValue<<std::endl;
        }
        if(i>maxIT){
            maxIT=i;
        }
    }
    //@HERE
    maxIndTop=maxIT;
    trace.info()<<maxIndTop<<std::endl;
    //Search the min ind point (different zeros) from bot normalized image
    unsigned int minIB=height_div-1;
    trace.info()<<height_div<<std::endl;
    for(j = 0; j < cols; j++){
        //trace.info()<<"------"<<std::endl;
        currentPixelValue=reliefImage.at<float>(height_div-1, j);
        i=height_div-1;
        while (currentPixelValue==0) {
            i=i-1;
            currentPixelValue=reliefImage.at<float>(i, j);
            //trace.info()<<currentPixelValue<<std::endl;
        }
        if(i<minIB){
            minIB=i;
        }
    }
    minIndBot=minIB;
    //@HERE
    trace.info()<<minIndBot<<std::endl;
    //crop top and bot
    cv::Rect myROI(0, maxIndTop, cols, minIndBot-maxIndTop);
    reliefImage = reliefImage(myROI);
    trace.info()<<minIndBot<<std::endl;
}

void
UnrolledMap::computeRGBImage(){
    cv::applyColorMap(toGrayscaleImageMinMax(), image, cv::COLORMAP_JET);
}
void
UnrolledMap::computeGRAYImage(){
    image=toGrayscaleImageMinMax();
}

/**GETTERS**/
cv::Mat
UnrolledMap::getNormalizedImage(){
    return reliefImage;
}

cv::Mat
UnrolledMap::getImage(){
    return image;
}

CylindricalPoint
UnrolledMap::getCPoint(unsigned int i){
    return CPoints.at(i);
}
