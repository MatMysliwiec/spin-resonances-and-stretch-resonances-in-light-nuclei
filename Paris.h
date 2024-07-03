#ifndef Paris_h
#define Paris_h
#include "TObject.h"
#include <vector>
#include <iostream>

#define WAVEFORMS

#define BRANCH Paris
#define DETECTOR ParisDetector

class Wave{
public:
  Float_t time;     // Signal time in ns
  Float_t ampl;     // Signal amplitude
  Float_t pdst;     // Signal pedestal
  Short_t form[400];  // Full waveform
  
};


class DETECTOR{
public:
Int_t id;       // Detector number
Float_t time;   // Signal time in ns
Int_t   qshort; // Q integral in short gate
Int_t   qlong;  // Q integral in long gate
#ifdef WAVEFORMS
Wave wave;  // Full waveform
#endif

bool ok(){
    return time > 0;
}

Int_t x(){      // x position
        return id/3;
}

Int_t y(){      // y position
        return id%3;
}
};



class BRANCH: public TObject{
public:
    static const int Nmod = 2;
    Double_t tref[Nmod];
    std::vector<DETECTOR> det;
   
Int_t mult()const{  // Multiplicity
    return det.size();
}

void Reset(){
    det.clear();
    memset(tref, 0, sizeof(tref));
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
