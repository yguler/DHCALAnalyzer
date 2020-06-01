#ifndef inputDHCAL_h
#define inputDHCAL_h

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
#include "TGraph2D.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "histDHCAL.h"
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

using namespace ROOT::Math;
using namespace std;


// ====================================================================
class inputDHCAL {

public:
   inputDHCAL();
   ~inputDHCAL();

   string  getInput(const string & dataset);  // return the input file name
   //void setDebugFlag();
   bool ClusterDistValid(TreeAnalysis & ev ,int z, int ind, int Rpc);
   int GetRPCClusterSize(TreeAnalysis & ev, int z);
   bool RPCClusterSizeValid(TreeAnalysis & ev);
   int  RPCClusterHitSizeValid(TreeAnalysis & ev, double x, double y, int z);
   bool LayerCheck(TreeAnalysis & ev, int z);
   int ActiveLayerCheck(TreeAnalysis & ev, int z);
   int LayerMin(TreeAnalysis & ev);
   int LayerMax(TreeAnalysis & ev);
   int RPCid(TreeAnalysis & ev);
   double set2DTrack(TreeAnalysis & ev);
   double SetFit(TreeAnalysis & ev);
   void line(int t, const double *p, double &x, double &y, double &z);
   double distance2(double x,double y,double z, const double *p);
   bool ClusterIsolation(TreeAnalysis & ev, double x1, double y1, int z1);
   //double SumDistance2(TGaph2D * g);
  
private:

   map<string, string> mydataset;
   map<string, double> myxsec;
//   map<string, int> mynentries;
   map<string, double> mynentries;
   double ncall;

   void DHcalNtuple();   

   //TTree * tree;
   string datasetname_; 
   double targetLumi_;
   double xsection_;
   double nentries_;
   double weight_;
   TGraph2D* grTrack = new TGraph2D();

};

#endif //~inputDHCAL_h
