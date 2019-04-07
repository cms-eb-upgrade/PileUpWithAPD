#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "Riostream.h"
#include <vector>
#include <iomanip>
#include <algorithm>



TRandom      r1;
const int    nbxmax_ = 3564;
const int    nBXs_ = 100;
const double tStep_ = 25./4.;
const int    nSamples_ = 20;



class EBChannel{
  
public:
  EBChannel(double);
  ~EBChannel();
  double energySignal()   { return energySignal_; };
  double energySpike()    { return energySpike_;  };
  double timeSpike()      { return timeSpike_;    };  
  void   GenerateBX(int);
  
private:
  TH1D  *hpdf_signal_;
  TH2D  *hpdf_spike_;
  double apd_probability_;
  double energySignal_;
  double energySpike_;
  double timeSpike_;
};



EBChannel::EBChannel(double eta)
{
  char hname[120];
  int indx;

  indx = int(fabs(eta-0.05) * 10.0);
  if(indx > 13) indx = 13;
  if(indx < 0)  indx = 0;
  sprintf(hname,"pupdf_%d",indx);
  TFile *fin1 = new TFile("pileupPDFs.root");
  hpdf_signal_ = (TH1D*)fin1->Get(hname);
  hpdf_signal_->SetDirectory(0);
  fin1->Close();

  TFile *fin2 = new TFile("pdf_apd_1M.root");
  TProfile *hp = (TProfile*)fin2->Get("h01");
  indx = int(fabs(eta)/1.5 * 85.) + 1;
  if(indx < 1)  indx = 1;
  if(indx > 85) indx = 85;
  apd_probability_ = hp->GetBinContent(indx);
  hpdf_spike_ = (TH2D*)fin2->Get("h03");
  hpdf_spike_->SetDirectory(0);
  fin2->Close();

  energySignal_ = 0;
  energySpike_  = 0;
  timeSpike_    = 0;
}



EBChannel::~EBChannel()
{
  delete hpdf_signal_;
  delete hpdf_spike_;
}



void EBChannel::GenerateBX(int NPU)
{
  energySignal_ = 0;
  energySpike_  = 0;
  timeSpike_  = 0;

  double rnd1, rnd2;
  for(int i=0; i<NPU; i++){
    rnd1 = hpdf_signal_->GetRandom();
    if(rnd1 > -10.9){
      energySignal_ += pow(10, rnd1);
    }
  }
  if( r1.Rndm() < (NPU * apd_probability_) ){
    hpdf_spike_->GetRandom2(rnd1,rnd2);
    energySpike_ += pow(10, rnd1) * 1.476;    // factor 1.476 accounts for pulse shaping by CATIA vs legacy FE
    if(rnd2 < 0)   rnd2 += 25.;
    if(rnd2 > 25.) rnd2 -= 25.; 
    timeSpike_ = rnd2;
    timeSpike_ -= 9.3; // offset wrt prompt signal
  }
}




class Pulse{

 public:
  Pulse(int);
  ~Pulse();
  double norm()  const { return norm_; };
  void SetNorm(double x)   { norm_ = x; };
  double Value(double);


 private:
  TGraph *grPulse_;
  TF1    *fPulse_;
  double tMin_;
  double tLimit_;
  double tMax_;
  double norm_;
  double offset_;
};



Pulse::Pulse(int iopt)
{
  TFile *file  = new TFile("TIA_ASIC_signal.root");
  TFile *file2 = new TFile("TIA_ASIC_signal_v2.root");
  if(iopt==0){
    grPulse_ = (TGraph*)file->Get("grTIA_Signal");
    fPulse_ = new TF1("fPulse_","[0]*exp(-x/[1])",0.,1e+6);
    fPulse_->SetParameters(2.21312e-05,1e+2);
    tMin_ = -0.4;
    tLimit_ = 200.;
    tMax_ = 1e+5;
    norm_ = 164.844;
    offset_ = 0;
  }else if(iopt==1){
    grPulse_ = (TGraph*)file->Get("grTIA_Spike");
    fPulse_ = new TF1("fPulse_","[0]*exp(-pow((x-[1])/[2],2))",0.,1e+6);
    fPulse_->SetParameters(3.79488e-01,5.44116e+01,9.62779);
    tMin_ = -0.4;
    tLimit_ = 51.3750;
    tMax_ = 1e+5;
    norm_ = 0.0235474;
    offset_ = +0.95;
  }else{
    grPulse_ = (TGraph*)file2->Get("grTIA_Signal_v2");
    fPulse_ = new TF1("fPulse_","[0]*exp(-x/[1])+[2]*exp(-x/[3])",0.,1e+6);
    fPulse_->SetParameters(1.65336e+05, 5.18932e+00, 4.33376e-02, 1.00000e+02);
    tMin_ =  0.0;
    tLimit_ = 50.0;
    tMax_ = 1e+5;
    //    norm_ = 1;
    norm_ = 1./503.888;
    offset_ = 0;
  }
  file->Close();
  file2->Close();
}



Pulse::~Pulse()
{
  delete fPulse_;
}


double Pulse::Value(double t)
{
  if(t < tMin_){
    return 0.;
  }else if(t < tLimit_){
    return grPulse_->Eval(t + offset_) * norm_;
  }else if(t < tMax_){
    return fPulse_->Eval(t + offset_) * norm_;
  }else{
    return 0.;
  }
}



class EBFrames{
  
 public:
  EBFrames();
  ~EBFrames();
  bool bxIsFilled(int i)  { return orbitScheme_[i]; };
  double eSignal(int i)   { return eSignal_[i];     };
  double tSignal(int i)   { return tSignal_[i];     };
  double eAPD(int i)      { return eAPD_[i];        };
  double tAPD(int i)      { return tAPD_[i];        };
  double eSample(int i)   { return eSample_[i];     };
  void   pushBX(double,double,double,double);
  void   pushSample(double);
  
 private:
  bool         orbitScheme_[nbxmax_];
  double       eSignal_[nBXs_];
  double       tSignal_[nBXs_];
  double       eAPD_[nBXs_];
  double       tAPD_[nBXs_];
  double       eSample_[nSamples_];
};



EBFrames::EBFrames()
{
  // Filling scheme:
  // 3564 = 12 x 297 = ((81b + 8e) x 3 + 30e) x 11 + ((81b + 8e) x 2 + 119e)

  for(int i=0; i<nbxmax_; i++){
    orbitScheme_[i] = false;
  }
  int idFilled = 0;
  int idAll    = 0;
  for(int itrain=0; itrain<11; itrain++){
    for(int igroup=0; igroup<3; igroup++){
      for(int ib=0; ib<81; ib++){
	orbitScheme_[idAll] = true;
        idFilled++;
        idAll++;
      }
      for(int ib=0; ib<8; ib++){
        idAll++;
      }
    }
    for(int ib=0; ib<30; ib++){
      idAll++;
    }
  }
  for(int igroup=0; igroup<2; igroup++){
    for(int ib=0; ib<81; ib++){
      orbitScheme_[idAll] = true;
      idFilled++;
      idAll++;
    }
    for(int ib=0; ib<8; ib++){
      idAll++;
    }
  }
  for(int ib=0; ib<119; ib++){
    idAll++;
  }

  for(int i=0; i<nBXs_; i++){
    eSignal_[i] = 0;
    tSignal_[i] = 0;
    eAPD_[i] = 0;
    tAPD_[i] = 0;
  }

  for(int i=0; i<nSamples_; i++){
    eSample_[i] = 0;
  }
}



EBFrames::~EBFrames()
{
}



void EBFrames::pushBX(double eS=0, double tS=0, double eA=0, double tA=0)
{
  for(int i=1; i<nBXs_; i++){
    eSignal_[i-1] = eSignal_[i];    
    tSignal_[i-1] = tSignal_[i];
    eAPD_[i-1] = eAPD_[i];    
    tAPD_[i-1] = tAPD_[i];    
  }
  eSignal_[nBXs_-1] = eS;
  tSignal_[nBXs_-1] = tS;
  eAPD_[nBXs_-1] = eA;
  tAPD_[nBXs_-1] = tA;
}



void EBFrames::pushSample(double e=0)
{
  for(int i=1; i<nSamples_; i++){
    eSample_[i-1] = eSample_[i];    
  }
  eSample_[nSamples_-1] = e;
}
