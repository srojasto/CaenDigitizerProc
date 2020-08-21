// The program extract the signal form binary file and average them
#include <fstream>
#include <string>
#include <iostream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TAxis.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"


// *************************************************
//        Variables and structure declarations
// *************************************************
Double_t ampRes = 0.244140625;// 1/4096;// amplitude resolution = 0.2 mV per bit (LSB) ;
const UInt_t baseLineSample=100;

struct h1Type {
    Double_t mean;
    Double_t rms;
};

struct grType {
  UInt_t n;
  Double_t * y;// = gr->GetY();
  Double_t * x;// = gr->GetY();
  UInt_t locmax;// = TMath::LocMax(n,y);
  UInt_t locmin;// = TMath::LocMax(n,y);
  Double_t vmax;// = y[locmax];
  Double_t vmin;// = y[locmax];
  Double_t tmin;// = x[locmax];
  Double_t tmax;// = x[locmax];
  Double_t tshift;// = y[locmax];
  UInt_t loctshift;
  Bool_t UnderThr = kFALSE;
};


// *************************************************
//        Function prototype
// *************************************************
std::ifstream::pos_type filesize(const char* filename);
h1Type GetBaseLine(float * tempBuffer);
grType GetGrProp(TGraph * gr, Double_t threshold, Int_t PulseSign , Double_t WStart, Double_t WEnd );


// *************************************************
//        Main program
// *************************************************
void RefTimeCorrection_v0(const Char_t * file = "wave_1.dat", Int_t sign = -1, Double_t Wstart =0, Double_t Wend = 204.8, Double_t thrs = 0){
  const UInt_t eventSz=1024;
  float  buffer[eventSz];
  FILE *ptr;

  if(sign >= 0) sign = 1;
  else sign = -1;
  printf("Sign %d\n",sign );

  // Get the size of the file
  // TODO set option for full or par of the file to be analized
  Int_t dataPointSz = sizeof buffer[0];
  cout <<  "Size of data point: " << sizeof buffer[0] <<" Bytes" << endl;
  Int_t fSize = filesize(file);
  UInt_t nEvents = fSize/(eventSz*dataPointSz);
  nEvents = 5000;
  printf("File size %i Bytes\nEvent size (%i)\nTotal Events(%i)\n\n", fSize, eventSz*dataPointSz, nEvents);

  ptr = fopen(file,"rb");  // r for read, b for binary
  //cout <<  "Size of file: " << sizeof ptr << endl;

  Double_t dt[eventSz];
  Double_t amplitude;
  float RawAmp;

  TGraph* grAv ;// = new TGraph(eventSz);
  TGraph* grRaw = new TGraph(eventSz);
  //TGraph* gr[eventSz] ;// = new TGraph(eventSz);

  Int_t nValid = 0;

  TH1 * h1BaseLineAll = new TH1D("h1BaseLineAll", "Base Line from all events;Base line level (mV);Entries",1024,0,1024);
  TH1 * h1RMSAll = new TH1D("h1RMSAll", "RMS from all events;RMS (mV);Entries",1000,0,10);

  TH2 * h2Signal = new TH2D("h2Signal", "Signal Average; time (ns);Amplitude (mV);", 1024,0,0.2*1023 ,4096*2,-4096*0.24414,4095*0.24414);
  TH2 * h2SignalShift = new TH2D("h2SignalShift", "Signal Average shifted; time (ns);Amplitude (mV);",1024*3,-0.2*1024,0.2*1023*2 ,4096*0.24414,-4095*0.24414,4096*0.24414);
  // h2RiseTime: Correlation between te start of the rise time at certain threshold and the amplitude
  TH2 * h2RiseTime = new TH2D("h2RiseTime", "Signal time rise vs amplitude; time (ns);Amplitude (mV);",512,-200,250, 2048,-500,100);
  TH1 * h1AmpMin = new TH1D("h1AmpMin", "Minimum amplitude;Amplitude (mV);Entries",2048,-500,100);
  TH1 * h1TimeMin = new TH1D("h1TimeMin", "Time of minimum;time of min. Amp. (ns);Entries",512,-200,250);
  TH1 * h1TimeShift = new TH1D("h1TimeShift", "Shifting time to the threshold;Shift time from peak (ns);Entries",512,-200,250);

  Double_t avSignal[eventSz];
  for (UInt_t i = 0; i < eventSz; i++){
    avSignal[i]=0;
    dt[i]=i*.2;
  }

  for(UInt_t event=0; event < nEvents; event++){
    grAv = new TGraph(eventSz);
    grAv -> GetYaxis() -> SetRange(-1000,1000);
    // Put one signal (1024*UInt32) in the buffer
    fread(buffer,sizeof(buffer),1,ptr);

    h1Type h1Pars = GetBaseLine(buffer);

    //Get base line of every event
    h1BaseLineAll -> Fill(h1Pars.mean);
    h1RMSAll-> Fill(h1Pars.rms);

    for (UInt_t i = 0; i < eventSz; i++){
      amplitude = buffer[i]*ampRes - h1Pars.mean;//
      grAv ->SetPoint(i,dt[i],amplitude);
    }

    // Localize the peak of the signal and its corresponding entrie number
    // Get signal properties
    grType sgProp = GetGrProp(grAv, thrs, sign, Wstart, Wend);

    // Cleaning samples
    // NOTE: Remove kFALSE to enable
    if (kFALSE &&  (sgProp.vmin < (h1Pars.mean* sign + 10)  || h1Pars.rms > 0.5) ){
    //if ( (sgProp.vmin > sign*thrs) || ( (sgProp.vmin) < sign*800) || h1Pars.rms > 1.5){
      delete grAv;
      continue; //saturated signals and noisy
    }

    nValid++;
    printf("Event number %i \r",nValid);


    h2RiseTime -> Fill(sgProp.tmin-sgProp.tshift, sgProp.vmin);
    h1AmpMin -> Fill(sgProp.vmin);
    h1TimeMin  -> Fill(sgProp.tmin);
    h1TimeShift -> Fill(sgProp.tshift);

    for (UInt_t i=0; i < sgProp.n; i++){
      h2Signal -> Fill(sgProp.x[i],sgProp.y[i]);
      // Shift graph
      sgProp.x[i] = sgProp.x[i] - sgProp.tshift;
      avSignal[i] += sgProp.y[i];
      h2SignalShift -> Fill(sgProp.x[i],sgProp.y[i]);
    }
    if(event < nEvents-1){ //keep the last graph for debuging
      delete grAv;
      }
    else{
      for (UInt_t i = 0; i < eventSz; i++){
        RawAmp = buffer[i];//&0x0000fff;
        //printf("point %d, buffer %d\n", i,buffer[i]);
        grRaw ->SetPoint(i,i,RawAmp);
      }
    }
  }
  cout << "\n\n";

  gStyle->SetOptTitle(kTRUE);
  //gStyle->SetPalette(kSolar);

  TCanvas * c2 = new TCanvas("c2","c2 title",1400,800);
  c2 -> Divide(2,1);
  c2 -> cd(1) -> SetLogz();
  h2Signal -> Draw("COLZ");
  c2 -> cd(2) ->SetLogz();
  h2SignalShift -> Draw("COLZ");

  TCanvas * c3 = new TCanvas("c3","Base line for all events",1400,800);
  c3 -> Divide (2,1);
  c3 -> cd(1) -> SetLogy();
  h1BaseLineAll -> Draw();
  c3 -> cd(2) -> SetLogy();
  h1RMSAll -> Draw();

  TCanvas * cProp = new TCanvas("cProp","Properties of signals",1400,1000);
  cProp -> Divide (2,2);

  cProp -> cd(1) -> SetLogz();
  h2RiseTime -> Draw("COLZ");
  cProp -> cd(2);
  h1AmpMin -> Draw();
  cProp -> cd(3);
  h1TimeMin  -> Draw();
  cProp -> cd(4);
  h1TimeShift -> Draw();


  TCanvas * c1 = new TCanvas("c1","c1 title",1400,800);
  grAv ->SetLineColor(kBlack);
  grAv ->SetLineWidth(1);
  grAv ->Draw("ALP*");


  TCanvas * cRaw = new TCanvas("cRaw","cRaw title",1400,800);
  grRaw ->SetLineColor(kBlack);
  grRaw ->SetLineWidth(1);
  grRaw ->Draw("ALP*");

}

h1Type GetBaseLine(float * tempBuffer){
  TH1D * h1BaseLine = new TH1D("h1BaseLine", "Base Line;Base line level (mV);Entries",1024,0,1024);
  for (UInt_t i = 0; i < baseLineSample; i++){ //Calculate base line
      h1BaseLine -> Fill(tempBuffer[i]*ampRes);
  }

  h1Type h1Result;
  h1Result.mean = h1BaseLine -> GetMean();
  h1Result.rms = h1BaseLine -> GetRMS();
  delete h1BaseLine;
  return h1Result;
}

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

grType GetGrProp(TGraph * gr, Double_t threshold = 0, Int_t PulseSign = -1, Double_t WStart =0, Double_t WEnd = 204.8){ // Graph properties
  grType graph;

  // Assign vertical (x) and horizontal (y) of the graph to pointers X and Y
  graph.n = gr ->GetN();
  graph.y = gr ->GetY();
  graph.x = gr ->GetX();

  // Locate array points for the maximum and minimum values
  graph.locmax = TMath::LocMax(graph.n,graph.y);
  graph.locmin = TMath::LocMin(graph.n,graph.y);

  // Get value of the minimum and maximum voltages and time
  graph.vmax = graph.y[graph.locmax];
  graph.vmin = graph.y[graph.locmin];
  graph.tmin = graph.x[graph.locmin];

  // Finding the location of a point over specific threshold
  for(UInt_t i = 0; i < graph.n; i++){
    if( (graph.x[i] > WStart) && (graph.x[i] < WEnd) ){

      if(graph.y[i] < threshold && PulseSign == -1){
        graph.tshift = graph.x[i];
        graph.loctshift = i;
        //printf("time %f, %f, %f\r", graph.x[i], graph.y[i], PulseSign * threshold);
        break;
      }

      if(graph.y[i] > threshold && PulseSign == 1){
        graph.tshift = graph.x[i];
        graph.loctshift = i;
        //printf("time %f, %f, %f \r", graph.x[i], graph.y[i], PulseSign * threshold);
        break;
      }
    }
    if(graph.x[i] > WStart && graph.y[i] < (PulseSign * threshold)) graph.UnderThr=kTRUE;

    // TODO handle signals below the threshold
  }

  return graph;
}
