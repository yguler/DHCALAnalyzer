//////////////////////////////////////////////////////////
//Code originally written by Yalcin Guler
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cassert>
#include <vector>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"

#include "inputDHCAL.h"
#include "histDHCAL.h"
#include "TreeAnalysis.h"
#include <TLorentzVector.h>

using namespace std;

string ANATYPE="DHCAL";

// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================

//g++ `root-config --cflags` -o DHCALAnalyzer TreeAnalysis.cc histDHCAL.cc inputDHCAL.cc anaDHCAL.cc `root-config --glibs`

class anaDHCAL {

public:
   anaDHCAL(string anaType_,string mydataset_,int isMC_,TFile * fout);
   ~anaDHCAL();

   void analyze(TreeAnalysis & ev, inputDHCAL & input, int isMC);
   void summary();

private:

   int ncall;
   string anaType_;
   string mydataset_;
   int isMC_;
//signal
  histoSet1 *h1NoCut;
};


// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================

anaDHCAL::anaDHCAL(string anaType, string mydataset,int isMC_,  TFile * fout) {
   anaType_=anaType;
   mydataset_=mydataset;
   std::cout<<"  anaDHCAL is created.  anaType= "<<anaType<<std::endl;
   ncall=0;

   string t=anaType;
//signal
   h1NoCut= new histoSet1(fout,t,"");
}

// ====================================================================
anaDHCAL::~anaDHCAL() {

}

// ====================================================================
void anaDHCAL::summary() {
  string anaType=anaType_;
  string mydataset=mydataset_;

  h1NoCut->summary(anaType,mydataset);
}

// ====================================================================

// ====================================================================
void anaDHCAL::analyze(TreeAnalysis & ev, inputDHCAL & input, int isMC_)
{
   ncall++;
   if(ncall<-1) {
      std::cout<<"yg anaDHCAL::analyze starts..."<<std::endl;
      //std::cout<<"s nhits ="<<ev.s_nhits<<std::endl;
   }

      double wt=1;
      if (ev.fNRPCClusters_NN > 4 && ev.NoiseIndex==0)
      	h1NoCut->process(ev,wt);

   return;

}   // end of analyze...

// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================

int main( int argc, char **argv ) {

  string anaType=ANATYPE;

  std::cout<<"  starting monojet ntuple analyzer ...  "<<std::endl;

  if ( argc < 3 ) {
    std::cerr << "Program need more than this parameter " << std::endl;
    std::cout << "Program need more than this parameter " << std::endl;
    std::cout <<"arg 1:  input dataset name"<<std::endl;
    std::cout <<"arg 2:  data/mc flag    (0=data, 1=mc)"<<std::endl;
    std::cout <<"arg 3:  number of events top rocess  (0 for all events)"<<std::endl;
    return 1;
  }

    
  string mydataset=argv[1];
  string foutName="DHCALResults_TrackFit/"+anaType+"_"+mydataset+".root";

  int isMC=std::stoi(argv[2]);

  int nev=std::stoi(argv[3]);
  // if(nev==0) nev=1000000000;
  if(nev<1) nev=1000000000;


  std::cout<<"anaType="<<anaType<<std::endl;
  std::cout<<"   isMC="<<isMC<<std::endl;
  std::cout<<"input dataset="<<mydataset<<std::endl;
  std::cout<<"output histogram file="<<foutName<<std::endl;
  std::cout<<"process "<<nev<<" events"<<std::endl;
//  std::cout<<"target luminosity "<<targetLumi<<std::endl;


  TFile *fout = new TFile(foutName.c_str(),"recreate");

  inputDHCAL input;

  //  setup analysis...
  anaDHCAL ana(anaType,mydataset,isMC,fout); 
  string inputFileName=input.getInput(mydataset);

  std::cout<<"InputFileName="<<inputFileName<<std::endl;
  //TFile f(inputFileName.c_str());
  TFile* f = TFile::Open( inputFileName.c_str() );
  TDirectory* dir = (TDirectory*)f->Get("DHCAL");  


  //string treeName=inputFileName+":/tree";
  //TDirectory * dir = (TDirectory*) f.Get(treeName.c_str());

  TTree * tree;
  dir->GetObject("DHCALEvent",tree);

  //   print out ntuple data structure...
  // tree->Print();

  double nentries = tree->GetEntriesFast();
  double nToUse=nentries;
//  int nToUse=nev;
//  if(nev<nentries) nToUse=nev;

  std::cout<<" Envents in nutple="<<nentries<<"  events to use="<<nToUse<<std::endl;


   //test ev(tree);
   TreeAnalysis ev(tree);
   std::cout<<"  yg calling LOOP"<<std::endl;
   // t.Loop();

   cout<<"  nev="<<nev<<endl;

   int nn=0;
   while ( ev.GetEntry(nn)){
       if(nn>(nev-1)) break;
       nn++;
       if(((nn%1)==100 && nn<1000) || ((nn%10000)==0 && nn<50000) || ((nn%100000)==0 && nn<500000) || (nn%1000000)==0)  {
          std::cout<<"nevents="<<nn<<std::endl;
       }
       // analyze event...
       ana.analyze(ev,input,isMC);
   }
         std::cout<<"\n total nevents = "<<nn<<std::endl;

  fout->Write();
  fout->Close();

  std::cout<<"======== Summary ========"<<std::endl;
  std::cout<<"    anaType: "<<anaType<<"   my dataset="<<mydataset<<"  total events="<<nn<<std::endl;
  std::cout<<"  "<<inputFileName<<std::endl;
  std::cout<<"  "<<std::endl;
  ana.summary();
  std::cout<<"========================="<<std::endl;

  return 0; 
}

