#ifndef Plastic_h
#define Plastic_h
#include "TObject.h"
#include <vector>
#include <iostream>

#define BRANCH Plastic
#define DETECTOR PlasticDetector

class DETECTOR{
public:
Int_t id;           // Detector number
Float_t time;       // Signal time in ns
Float_t ampl;       // Signal amplitude
Float_t pdst;       // Signal pedestal
Short_t form[400];  // Full waveform

Int_t kratta()const{             // id of the kratta module
    return id/4;
}

bool ok()const{
    return time > 0;
}

bool intime()const{
    return time > 620 && time < 720;
}

};

class BRANCH: public TObject{

public:
std::vector<DETECTOR> det;
    
Int_t mult()const{ // Multiplicity
    return det.size();
}

void Reset(){
    det.clear();
}

void Add(DETECTOR & p){
    det.push_back(p);
}

BRANCH(){}

~BRANCH(){}

ClassDef(BRANCH, 1)
};

#undef BRANCH
#undef DETECTOR

#endif
