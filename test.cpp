 /**************************************************************************************************************************/ 
void
DefectSegmentationUnroll::test() {

  int decreaseHit=4;
  int current_resX,current_resY;
  int theta =159;
  int hauteur=200;
  if(unrolled_surface[hauteur][theta].empty()){
    trace.info()<<"cette case est vide "<<std::endl;
  }else{
      trace.info()<<"cette case est pas vide "<<std::endl;
  }
//  for(unsigned int i = 0; i < height_div; i++){
//    for(unsigned int j = 0; j < angle_div; j++){
      current_resX=angle_div;
      current_resY=height_div;
      trace.info()<<"res de base X : "<<current_resX<<std::endl;
      trace.info()<<"res de base Y : "<<current_resY<<std::endl;
      std::vector<unsigned int > accumulate_vector;
//      t=getPointsVectorFromOtherResolution(i,j,current_resX,current_resY);
      //virtual decrease of unrolled surface
//      if(t.empty()){*/
        //current_resX*=1./decreaseHit;
        //current_resY*=1./decreaseHit;
        int counter=0;
        //loop to find i and j of left-top corner of the square containing P
        for(unsigned int i = 0; i < height_div; i+=decreaseHit){
            for(unsigned int j = 0; j < angle_div; j+=decreaseHit){
              //trace.info()<< "i :"<<i<<std::endl;
              //trace.info()<< "j :"<<j<<std::endl;
              //check if P is in square
              if((hauteur>=i && hauteur<i+decreaseHit)&&(theta>=j&&theta<j+decreaseHit)){
                //accumulate vector of all cells in the square
                /*for(unsigned int k = i; k < i+decreaseHit; k++){
                  for(unsigned int h = j; h < j+decreaseHit; h++){
                    if((k<height_div)&&(h<angle_div)){
                      accumulate_vector.reserve(accumulate_vector.size() + unrolled_surface[k][h].size());
                      accumulate_vector.insert(accumulate_vector.end(), unrolled_surface[k][h].begin(),  unrolled_surface[k][h].end());
                      counter+=1;
                    }

                  }
                }*/
                trace.info()<<i<<" "<<j<<std::endl;
              }
            }
        }
        int newhauteur;
        int newtheta;
        //trace.info()<<"res 1/2 X : "<<current_resX<<std::endl;
        //trace.info()<<"res 1/2 Y : "<<current_resY<<std::endl;
        //newhauteur=roundf(((double((current_resY-1))/(height_div-1))*(hauteur-(height_div-1)))+(current_resY-1));
        //newtheta=roundf(((double((current_resX-1))/(angle_div-1))*(theta-(angle_div-1)))+(current_resX-1));
        newhauteur=(hauteur/decreaseHit)*decreaseHit;
        newtheta=(theta/decreaseHit)*decreaseHit;
        //trace.info()<<(hauteur/decreaseHit)<<std::endl;
        //trace.info()<<(theta/decreaseHit)<<std::endl;
        trace.info()<<"nex hauteur :"<<newhauteur<<" new theta : "<<newtheta<<std::endl;

        //trace.info()<<"nb cells considered "<<counter<<std::endl;
        //for(std::vector<unsigned int>::iterator it = std::begin(accumulate_vector); it != std::end(accumulate_vector); ++it) {
            //meshTest.addVertex(pointCloud.at(*it));
        //    trace.info()<<*it<<std::endl;
        //}
/*      }

    }
  }*/
}

 /*************************************************************************************************************************/ 
 * *********** multiple patch 
 */
double estimateRadii = mpCurrent.height * currentCoefficient.first + currentCoefficient.second;
    double deltaDist= mpCurrent.radius - estimateRadii;
    //3 patch research
    //TODO : generalize
    if(deltaDist<0 ){
      currentCoefficient2 = computeEq(i,secondSearchRadius,secondPatchAngle,kdtree);
      double estimateRadii2 = mpCurrent.height * currentCoefficient2.first + currentCoefficient2.second;
      double deltaDist2= mpCurrent.radius - estimateRadii2;
      if(deltaDist2>deltaDist){
        currentCoefficient3 = computeEq(i,thirdSearchRadius,thirdPatchAngle,kdtree);
        double estimateRadii3 = mpCurrent.height * currentCoefficient3.first + currentCoefficient3.second;
        double deltaDist3= mpCurrent.radius - estimateRadii3;
        if(deltaDist3>deltaDist2){
          //currentCoefficient4 = computeEq(i,fouthSearchRadius,fouthPatchAngle,kdtree);
          //double estimateRadii4 = mpCurrent.height * currentCoefficient4.first + currentCoefficient4.second;
          //double deltaDist4= mpCurrent.radius - estimateRadii4;
          //if(deltaDist4>deltaDist3){
          //  coefficients[i]=currentCoefficient4;
          //}else{
          coefficients[i]=currentCoefficient3;
          //}
        }else{
          coefficients[i]=currentCoefficient2;
        }
      }else{
        coefficients[i]=currentCoefficient;
      }
    }else{
      coefficients[i]=currentCoefficient;
    }

 /*************************************************************************************************************************/   
/**
search the point with min or max height in function of angle
**/
//unsigned int pointSearchByAngle(std::vector<double>& xs, std::vector<double>& ys, std::vector<unsigned int>& indP,CylindricalPoint mpCurrent, std::pair<double, double> coef);
unsigned int
DefectSegmentationUnroll::pointSearchByAngle(std::vector<double>& xs, std::vector<double>& ys, std::vector<unsigned int>& indP,CylindricalPoint mpCurrent,std::pair<double, double> coef){
  double a=coef.first;
  double b=coef.second;
  double angle=std::atan(b/a);
  //trace.info()<<"angle : "<<angle<<std::endl;
  unsigned int indFound;
  if(angle>0){
    double min=INT_MAX;
    for (unsigned int i = 0; i < ys.size(); ++i){
      if(xs[i]<min){
        min=xs[i];
        indFound=indP[i];
      }
    }
  }else{
    double max=INT_MIN;
    for (unsigned int i = 0; i < ys.size(); ++i){
      if(xs[i]>max){
        max=xs[i];
        indFound=indP[i];
      }
    }
  }
  //trace.info()<<"ind found : "<<indFound<<std::endl;
  return indFound;
}



 /*************************************************************************************************************************/ 


/**
Compute coef of the shifted line, line is shifted on y axis
**/
void tresholdPatchbyHeight(std::vector<double>& xs, std::vector<double>& ys, std::vector<unsigned int>& indP,CylindricalPoint mpCurrent,std::pair<double, double> coef);
   
void 
DefectSegmentationUnroll::tresholdPatchbyHeight(std::vector<double>& xs, std::vector<double>& ys, std::vector<unsigned int>& indP,CylindricalPoint mpCurrent,std::pair<double, double> coef){
  std::vector<double> new_xs;
  std::vector<double> new_ys;
  std::vector<unsigned int> new_indP;
  double a=coef.first;
  double b=coef.second;
  double angle=std::atan(b/a);
  /*if(angle>0){
    for (unsigned int i = 0; i < ys.size(); ++i){
      if(xs[i]<mpCurrent.radius){

      }
    }
  }else{

  }*/
  
  xs=new_xs;
  ys=new_ys;
  indP=new_indP;
}


 /*************************************************************************************************************************/ 
/**
 Compute coef of the shifted line, line is shifted on y axis
**/
void tresholdPatchbyRadius(std::vector<double>& xs, std::vector<double>& ys, std::vector<unsigned int>& indP,CylindricalPoint mpCurrent);
 
 void 
DefectSegmentationUnroll::tresholdPatchbyRadius(std::vector<double>& xs, std::vector<double>& ys, std::vector<unsigned int>& indP,CylindricalPoint mpCurrent){
  std::vector<double> new_xs;
  std::vector<double> new_ys;
   std::vector<unsigned int> new_indP;

  for (unsigned int i = 0; i < ys.size(); ++i){
    if(ys[i]<mpCurrent.radius){
      new_xs.push_back(xs[i]);
      new_ys.push_back(ys[i]);
      new_indP.push_back(indP[i]);
    }
  }
  xs=new_xs;
  ys=new_ys;
  indP=new_indP;
}

/*************************************************************************************************************************/ 
/********************CODE POUR GROSSIR UN PATCH DUN COTE**************************************************/ 
/*************************************************************************************************************************/ 
//compute delta distance 
  double deltadist;
  double estimateRadii = mpCurrent.height * c.coefficients.first + c.coefficients.second;
  if(c.coefficients.second == 0.0){
      deltadist = 0;
  }else{
      deltadist = mpCurrent.radius - estimateRadii;
  }
  //If negative delta distance
  if(deltadist<0){
    //search right or left poiint
    unsigned int id=pointSearchByAngle(lengthForEstimate,radiusForEstimate,indForEstimate,mpCurrent, c.coefficients);
    CylindricalPoint mpCyl = myPoints.at(id);
    Z3i::RealPoint mpCart = pointCloud.at(id);
    std::vector<int> pointIdx2;
    //compute distance on the hieght between current point and founded point
    double radiusTosearch=abs(mpCyl.height-mpCurrent.height);
    pcl::PointXYZ mpPCL(mpCart[0], mpCart[1], mpCart[2]);
    if ( kdtree.radiusSearch (mpPCL, searchRadius, pointIdx2, pointRadiusSquaredDistance) > 0 ){
      for (unsigned int idx = 0; idx < pointIdx2.size (); ++idx){
        unsigned int foundedIndex = pointIdx2.at(idx);
        CylindricalPoint mpFound = myPoints.at(foundedIndex);
        double angleDiff = std::abs(mpFound.angle - mpCyl.angle);
        //filtre en patch rectangulaire
        if(angleDiff > patchAngle/2 && 2*M_PI - angleDiff > patchAngle / 2){
          continue;
        }
        //fill the patches vector
        radiusForEstimate.push_back(mpFound.radius);
        lengthForEstimate.push_back(mpFound.height);
        indForEstimate.push_back(foundedIndex);
      }
    }
    unique( radiusForEstimate.begin(), radiusForEstimate.end() );
    unique( lengthForEstimate.begin(), lengthForEstimate.end() );
    unique( indForEstimate.begin(), indForEstimate.end() );
    c= Regression::PurgedlinearRegression(lengthForEstimate, radiusForEstimate, indForEstimate);
  }
//c= Regression::PurgedlinearRegression(lengthForEstimate, radiusForEstimate, indForEstimate);

/*************************************************************************************************************************/ 
/********************CODE POUR CALCULER LE TREHOLD DES DELTADIST**************************************************/ 
/*************************************************************************************************************************/ 

/**
compute rosin rehold on vectorgiven in paramtyers
**/
double findDistanceThresholdRosin();

double 
DefectSegmentationUnroll::findDistanceThresholdRosin(){
  //build histogram
  double res = binWidth; //must be configurable?
  //trace.info()<<"binWidth:"<<binWidth<<std::endl;
  double maxValue = *std::max_element(distances.begin(), distances.end());
  double minValue = *std::min_element(distances.begin(), distances.end());
  double range = maxValue - minValue;
  trace.info()<<"range de l'histo "<< range <<std::endl;
  int nbInterval = range / res;
  trace.info()<<"taille de l'histo : "<< nbInterval <<std::endl;
  std::vector<int> histogram(nbInterval, 0);
  for(unsigned int i = 0; i < myPoints.size(); i++){
      int index = (distances[i] - minValue)/res;
      histogram[index]+=1;
  }
  
  std::vector<int>::iterator maxFreq = std::max_element(histogram.begin(), histogram.end());

  int maxFreqValue = *maxFreq;
  int maxFreqIndex = std::distance(histogram.begin(), maxFreq);
  assert(maxFreqValue == histogram.at(maxFreqIndex));

  unsigned int lastIndex = histogram.size() - 1;
  int lastValue = histogram.at(lastIndex);


  for(unsigned int i = maxFreqIndex; i < histogram.size(); i++){
      if(histogram.at(i) == 0){
          lastIndex = i;
          lastValue = 0;
          break;
      }
  }
  double valueDiff = lastValue - maxFreqValue;
  double valueDiff2 = valueDiff *valueDiff;
  double indexDiff = lastIndex - maxFreqIndex;
  double indexDiff2 = indexDiff * indexDiff;
  double bestThresIndex = maxFreqIndex;
  double bestDist = 0;
  //line between maxFreq and last element of historgram
  double a = (lastValue - maxFreqValue)*1.0/(lastIndex - maxFreqIndex);
  double b = maxFreqValue - a * maxFreqIndex;
  //for (std::vector<int>::iterator it = maxFreq ; it != histogram.end(); ++it){
  for (unsigned int i = maxFreqIndex; i < lastIndex; i++){
      //trace.info()<<valueDiff *  i - indexDiff*histogram.at(i) + maxFreqIndex*lastValue - maxFreqValue * lastIndex<<std::endl;
      double dist = std::abs(valueDiff *  i - indexDiff*histogram.at(i) + maxFreqValue*lastIndex - maxFreqIndex * lastValue)/
          sqrt(valueDiff2 + indexDiff2 );
      if(dist > bestDist){
          bestDist = dist;
          bestThresIndex = i;
      }
  }
  int bestVal = histogram.at(bestThresIndex);
  double x2 = (bestVal + bestThresIndex/a - b)/(a + 1/a);
  double y2 = b + a*x2;

  std::vector<std::pair<double, double>> forPlot;
  std::pair<double, double> maxPoint(maxFreqIndex * res + minValue, maxFreqValue);
  std::pair<double, double> lastPoint(lastIndex * res + minValue, lastValue);
  std::pair<double, double> bestPoint(bestThresIndex * res + minValue, bestVal);
  std::pair<double, double> projBestPoint(x2* res + minValue, y2);

  forPlot.push_back(maxPoint);
  forPlot.push_back(lastPoint);
  forPlot.push_back(bestPoint);
  forPlot.push_back(projBestPoint);
  //IOHelper::export2Text(forPlot, "pointFile");

  //histogram
  std::vector<std::pair<double, double>> histForPlot;
  for(unsigned int i = 0; i< histogram.size(); i++){
      std::pair<double, double> aBin(i* res + minValue, histogram.at(i));
      histForPlot.push_back(aBin);
  }

  IOHelper::export2Text(histForPlot, "hist2d");
	trace.info()<<"threshold: "<< bestThresIndex*res + minValue<<std::endl;
  
  return bestThresIndex*res + minValue;

}

/*************************************************************************************************************************/ 
/********************CODE POUR CALCULER LES DELTADIST SEUILLER PAR ROSIN**************************************************/ 
/*************************************************************************************************************************/

/**
filtre delta distance by rosin
**/
void computeDeltaDistancesRosin();
   

void
DefectSegmentationUnroll::computeDeltaDistancesRosin(){
  //compute delta distance forall points
  computeDeltaDistances();
  //compute treshold by rosin
  double th = findDistanceThresholdRosin();
  //loop on delta distance
  for(unsigned int i = 0; i < distances.size(); i++){
      if(distances[i] <= th){
        distances[i]=0.0;
      }else{
        distances[i]=1.0;
      }
  }

}