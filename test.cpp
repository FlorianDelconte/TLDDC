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
