#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "Riostream.h"
#include <algorithm>

TRandom r1;


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



void Example_1EBch_PUandSpikes(int nevents=1000000, double eta=1.45, double npu=200.)
{
  EBChannel *ch = new EBChannel(eta);

  TH1D *hS = new TH1D("hS","Energy in a crystal", 1000, 0, 10.);
  TH1D *hA = new TH1D("hA","Energy in APD",       1000, 0, 10.);
  TH1D *hT = new TH1D("hT","Time in APD",         100, -25, 25.);
  TH1D *hE = new TH1D("hE","Total energy",        1000, 0, 10.);
  
  for(int ievt=0; ievt<nevents; ievt++){
    if(ievt%10000==0){
      cout << "Event " << ievt << " out of " << nevents << endl;
    }
    int NPU = r1.Poisson(npu);
    ch->GenerateBX(NPU);
    hS->Fill( ch->energySignal() );
    hA->Fill( ch->energySpike() );
    if(ch->energySpike() > 1e-9)
      hT->Fill( ch->timeSpike() );
    hE->Fill( ch->energySignal() + ch->energySpike() );
  }

  TH1D *hP1 = (TH1D*)hE->Clone("hP1");
  TH1D *hP2 = (TH1D*)hE->Clone("hP2");
  hP1->Reset();
  hP2->Reset();
  double nNow = 0;
  for(int ib=hE->GetNbinsX()+1; ib>=0; ib--){
    nNow += hE->GetBinContent(ib);
    hP1->SetBinContent(ib, nNow);
    hP2->SetBinContent(ib, hE->GetEntries());
  }
  TGraphAsymmErrors *gr = new TGraphAsymmErrors(hP1,hP2);
  gr->SetName("gr");
  gr->SetTitle("Probability to have a hit above threshold");
  gr->GetXaxis()->SetTitle("threshold (GeV)");
  gr->GetYaxis()->SetTitle("probability");
  gr->GetXaxis()->SetLimits(hE->GetBinCenter(1),hE->GetBinCenter(hE->GetNbinsX()));
  gr->Draw("APL");
  
  TFile *fout = new TFile("output.root","recreate");
  hS->Write();
  hA->Write();
  hT->Write();
  hE->Write();
  gr->Write();
  fout->Close();
}

