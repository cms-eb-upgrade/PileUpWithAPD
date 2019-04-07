#include "PileUpWithAPD.h"


void Example02()
{
  double eta = 1.4;
  double npu = 200.;
  double rophase = 1.0;
  
  EBChannel *ch = new EBChannel(eta);
  EBFrames  *frame = new EBFrames();
  Pulse *pulseS = new Pulse(0);
  Pulse *pulseA = new Pulse(1);
  
  int bxNow  = 0;
  cout << " Press Enter to get next BX or Ctrl-C to quit..." << endl;

  while(true){

    int nOrbit  = bxNow / 3564;
    int bxLocal = bxNow % 3564;
    
    int NPU = 0;
    if(frame->bxIsFilled(bxLocal)){
      NPU = r1.Poisson(npu);
    }
  
    ch->GenerateBX(NPU);
    frame->pushBX(ch->energySignal(), 0., ch->energySpike(), ch->timeSpike());

    // ADC samples of 160 MHz (in GeV)
    for(int is=0; is<4; is++){
      double e = 0;
      double tSample = 6.25 * is;
      for(int ibx=0; ibx<100; ibx++){
	double tS = tSample + (99 - ibx) * 25.0 + rophase + frame->tSignal(ibx);
	if(tS>0){
	  e += frame->eSignal(ibx) * pulseS->Value(tS);
	}
	double tA = tSample + (99 - ibx) * 25.0 + rophase + frame->tAPD(ibx);
	if(tA>0){
	  e += frame->eAPD(ibx) * pulseA->Value(tA);
	}
      }
      frame->pushSample(e);
    }
	
    
    cout << "  Orbit=" << nOrbit
	 << "  BX=" << bxLocal
	 << "  Energy in BX = " <<  frame->eSignal(99) + frame->eAPD(99)
	 << endl;
    cout << "          Sample 0:   " << frame->eSample(16) << endl;
    cout << "          Sample 1:   " << frame->eSample(17) << endl;
    cout << "          Sample 2:   " << frame->eSample(18) << endl;
    cout << "          Sample 3:   " << frame->eSample(19) << endl;
    cin.get();
    bxNow++;
  }
}
