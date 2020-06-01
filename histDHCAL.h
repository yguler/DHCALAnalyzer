#ifndef histDHCAL_h
#define histDHCAL_h


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
#include "TH3D.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"

#include "TreeAnalysis.h"

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
using namespace ROOT::Math;
using namespace std;

// ====================================================================
class histoSet1{

public:
  histoSet1(TFile * fout, string anaType,  string cutname);
  ~histoSet1();

  void process(TreeAnalysis & ev, double weight);
  void summary(string anaType, string mydataset);
  vector<double>  PileUpWeights();

private:

  // Histograms with names
  std::map<std::string, TH1D*> histo1D;
  std::map<std::string, TH1D*>::iterator histo1Diter;

  std::map<std::string, TH2D*> histo2D;
  std::map<std::string, TH2D*>::iterator histo2Diter;

  std::map<std::string, TH3D*> histo3D;
  std::map<std::string, TH3D*>::iterator histo3Diter;
  
  std::map<std::string, TGraph2D*> graph2D;
  std::map<std::string, TGraph2D*>::iterator graph2Diter;
  double nevents_unweighted;
  double nevents_weighted;
  double nevents_weighted_New;
  string dirname;
  TGraph2D* grTrack = new TGraph2D();
  TGraph *gr1 = new TGraph();
};

#endif  //~histDHCAL_h
