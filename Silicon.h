#ifndef Silicon_h
#define Silicon_h
#include "TObject.h"
#include <vector>
#include <iostream>

#define BRANCH Silicon
#define DETECTOR SiliconDetector

class DETECTOR{
public:
Int_t id;                   // Detector number
Float_t time30;               // Signal time 30 %
Float_t time80;               // Signal time 80 %
Float_t ampl;               // Signal amplitude


bool ok(){ return 
    time30 > 0;
//     && time30<3840.&& 
//     ampl > 0.  && ampl < 3840. 
//     && time80 > 0.&& time80<3840.;
}

Float_t deltaT(){return time30 - time80;}

};




class BRANCH: public TObject{

public:
    Float_t tplsi;
    std::vector<DETECTOR> det;

Int_t mult()const{ // Multiplicity
    return det.size();
}

void Reset(){
    tplsi = 0;
    det.clear();
}

void Add(DETECTOR & p){
    det.push_back(p);
}

BRANCH(){}

~BRANCH(){}

ClassDef(BRANCH, 2)
};

#undef BRANCH
#undef DETECTOR

#endif
