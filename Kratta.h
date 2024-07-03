#ifndef Kratta_h
#define Kratta_h
#include "TObject.h"
#include <vector>
#include <iostream>

#define WAVEFORMS

#define BRANCH Kratta
#define DETECTOR KrattaDetector

class PD0{
public:
  Float_t time;     // Signal time in ns
  Float_t ampl;     // Signal amplitude
  Float_t pdst;     // Signal pedestal
#ifdef WAVEFORMS
  Short_t form[512];  // Full waveform
#endif

// Build-in functions

bool ok()const{
    return time > 0;
}

bool intime()const{
    return (time > 1040 && time < 1100); // trig = 2
//    return (time > 600 && time < 900); // trig = 3
//    return (time > 1200 && time < 1250); // trig = 4
}

};

class PD1{
public:
  Float_t time;     // Signal time in ns
  Float_t ampl;     // Signal amplitude
  Float_t pdst;     // Signal pedestal
#ifdef WAVEFORMS
  Short_t form[1024]; // Full waveform
#endif

// Build-in functions

bool ok()const{
    return time > 0;
}

bool intime()const{
    return (time > 1800 && time < 2000); // trig = 2
//    return (time > 1100 && time < 1600); // trig = 3
//    return (time > 2000 && time < 2200); // trig = 4
}

};

class PD2{
public:
  Float_t time;     // Signal time in ns
  Float_t ampl;     // Signal amplitude
  Float_t pdst;     // Signal pedestal
#ifdef WAVEFORMS
  Short_t form[1024]; // Full waveform
#endif

// Build-in functions

bool ok()const{
    return time > 0;
}

bool intime()const{
    return (time > 1950 && time < 2050); // trig = 2
//    return (time > 1100 && time < 1600); // trig = 3
//    return (time > 2000 && time < 2200); // trig = 4
} 
};

class DETECTOR{
public:
Int_t   id;         // Detector number
PD0     pd0;        // PD0
PD1     pd1;        // PD1
PD2     pd2;        // PD2

// Build-in functions

Int_t x()const{   // x position
    return id%8;
}

Int_t y()const{   // y position
    return id/8;
}
Int_t plastic(Int_t i)const{    // id of i-th plastic (i = 0..3)
    return id*4 + i;
}

bool ok(){
    return pd0.ok() || pd1.ok() || pd2.ok();
}

bool pd0pd1()const{
    return pd0.intime() && pd1.intime();
}

bool pd1pd2()const{
    return pd1.intime() && pd2.intime();
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

~BRANCH(){};

ClassDef(BRANCH, 1)
};

#undef BRANCH
#undef DETECTOR

#endif
