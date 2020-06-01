#include "inputDHCAL.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cassert>
#include <vector>
#include <time.h>

#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"
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
#include <TLorentzVector.h>
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

//g++ `root-config --cflags` -o DHCALAnalyzer TreeAnalysis.cc histDHCAL.cc inputDHCAL.cc anaDHCAL.cc  `root-config --glibs`
// ====================================================================

inputDHCAL::inputDHCAL() {
   //std::cout<<"  inputDHCAL is created "<<std::endl;
   ncall=0;

   DHcalNtuple();
}

// ====================================================================
inputDHCAL::~inputDHCAL() {

}


// ====================================================================
string inputDHCAL::getInput(const string & dataset){

  //  find and set tree name...
  //
  //  arg1: string-  data set name
  string filename=mydataset[dataset];

  // save weight...
  datasetname_=dataset;
  nentries_=mynentries[dataset];
   std::cout<<" yg inputDHCAL::getInput - dateset name="<<dataset<<std::endl;


 return filename;
}
inputDHCAL input;
bool inputDHCAL::ClusterDistValid(TreeAnalysis & ev,int z, int ind, int Rpc){
	bool cldist = false;
	bool clup =false; 
	bool cldown =false;
	int ClDistValid[55][100]={0};
		int zz=0;
		for (int j=0; j<ev.fNRPCClusters_NN; ++j) {
			double x2=ev.RPCClusters_NN_x[j]; double y2=ev.RPCClusters_NN_y[j]; int z2=((int)ev.RPCClusters_NN_z[j]);
			double x1=ev.RPCClusters_NN_x[ind]; double y1=ev.RPCClusters_NN_y[ind];
			if(z == z2 || j == ind) continue;
			bool Rpc0 = (y2>3 && y2<=28) && (x2>3 && x2<93);
			bool Rpc1 = (y2>35 && y2<=60) && (x2>3 && x2<93);
			bool Rpc2 = (y2>65 && y2<93) && (x2>3 && x2<93); 
				if ((Rpc==0 && Rpc0) || (Rpc==1 && Rpc1) || (Rpc==2 && Rpc2)){
				//cout <<z1 <<".Layer ile "<<z2<<".LAyer "<<" x2 = "<<x2<<" y2 = "<<y2<<endl;
                        		double cldist2=sqrt(pow(x1-x2,2.)+pow(y1-y2,2.)); 
					if ( cldist2<=3 ){
						clup = true;
						ClDistValid[z2][zz]=1; 
						zz++;
					}
				}  
		}
	return clup;	
}

bool inputDHCAL::ClusterIsolation(TreeAnalysis & ev, double x1, double y1, int z1){
        int hit=0;
                bool iso = false;
                for (int j=0; j<ev.fNRPCHits; ++j) {
                        double x2=ev.RPCHits_x[j]; double y2=ev.RPCHits_y[j]; int z2=((int)ev.RPCHits_z[j]);
                        if (z1==z2){
                                double dist=sqrt(pow(x1-x2,2.)+pow(y1-y2,2.));
                                if (dist>3 && dist<=7) hit++;
                                        //cout<< j<<" RPCHits_x = "<<x2<<" RPCHits_y = "<<y2<< " RPCHits_z = "<<z2<<endl;
                        }
                }
		if (hit>0) iso = true; 
                return iso;
}

int inputDHCAL::GetRPCClusterSize(TreeAnalysis & ev, int z){
	int NClPerLayers[55]={0};	
	int ncl=0;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
		NClPerLayers[z1]+=1;
        }
	ncl = NClPerLayers[z]; 
	return ncl;
}

int inputDHCAL::RPCClusterHitSizeValid(TreeAnalysis & ev, double x1, double y1, int z1){
	int hit=0;
		int zz=0;
		for (int j=0; j<ev.fNRPCHits; ++j) {
                	double x2=ev.RPCHits_x[j]; double y2=ev.RPCHits_y[j]; int z2=((int)ev.RPCHits_z[j]);
			if (z1==z2){
				double dist=sqrt(pow(x1-x2,2.)+pow(y1-y2,2.));
				if (dist<=2) hit++; 
                			//cout<< j<<" RPCHits_x = "<<x2<<" RPCHits_y = "<<y2<< " RPCHits_z = "<<z2<<endl;
			}
        	}
		return hit;
}

bool inputDHCAL::RPCClusterSizeValid(TreeAnalysis & ev){
        int NClHitsPerLayers[55]={0};
        bool cl = false;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
                NClHitsPerLayers[z1]+=1;
        }

        for(int i1=0; i1<38; ++i1){
		if (NClHitsPerLayers[i1]>1) cl=true; 
        }
        return cl;
}

bool inputDHCAL::LayerCheck(TreeAnalysis & ev, int z){
	bool lcheck=false;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
                if (z1==z) {lcheck=true; break;}
        }
	return lcheck;
}

int inputDHCAL::ActiveLayerCheck(TreeAnalysis & ev, int z){
        int lcheck=-1;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
		if (((y1>3 && y1<=28) && (x1>3 && x1<93)) && z1==z){lcheck=0; break;}
		else if (((y1>35 && y1<=60) && (x1>3 && x1<93)) && z1==z){lcheck=1; break;}		
		else if (((y1>65 && y1<93) && (x1>3 && x1<93)) && z1==z){lcheck=2; break;}
        }
        return lcheck;
}
int inputDHCAL::LayerMin(TreeAnalysis & ev){
	int min = 9999;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
                if (z1<min) min=z1;
        }
	return min;
}

int inputDHCAL::LayerMax(TreeAnalysis & ev){
	int max=-1;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
                if (z1>max) max=z1;
        }
	return max;
}

int inputDHCAL::RPCid(TreeAnalysis & ev){
        int rpcid=-1;
        for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
        	if ((y1>3 && y1<=28) && (x1>3 && x1<93)) {rpcid=0;} 
		else if ((y1>35 && y1<=60) && (x1>3 && x1<93)) {rpcid=1;}
		else if ((y1>65 && y1<93) && (x1>3 && x1<93)) {rpcid=2;}
        }
        return rpcid;
}

void inputDHCAL::line(int t, const double *p, double &x, double &y, double &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

double inputDHCAL::distance2(double x,double y,double z, const double *p) {
   // distance line point is D= | (xp-x0) cross  ux |
   // where ux is direction of line and x0 is a point in the line (like t = 0)
   XYZVector xp(x,y,z);
   XYZVector x0(p[0], p[2], 0. );
   XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
   XYZVector u = (x1-x0).Unit();
   double d2 = ((xp-x0).Cross(u)) .Mag2();
//cout <<" x = "<<x<<" y = "<<y<<" z = "<<z<<" p[0] = "<<p[0]<<" p[1] = "<<p[1]<<" p[2] = "<<p[2]<<" p[3] = "<<p[3]<<" d2 = "<<d2<<endl;
   return d2;
//https://root.cern.ch/doc/master/classTVector3.html
//https://root.cern.ch/root/htmldoc/guides/users-guide/PhysicsVectors.html
}

// ==================================================================
void inputDHCAL::DHcalNtuple(){

  string fileDirectoryRe = "root://cmsxrootd.fnal.gov//store/user/";
  string fileDirectoryM17GJets ="root://cmseos.fnal.gov//store/group/";
  string fileDirectorydatav1 ="/uscms_data/d3/ygule/";
  string fileDirectorylx ="/afs/cern.ch/work/y/yalcin/public/hgcal_March2018/ntuples/";
  string fileDirEOSDhcal ="root://eoscms.cern.ch//eos/cms/store/user/yalcin/DHCAL/data/";

////////////////
  mydataset["pionBeam_Run650357_4GeV"] = fileDirEOSDhcal+"DHCAL_650357.root";
  mydataset["MuonBeam_Run630157_32GeV"] = fileDirEOSDhcal+"DHCAL_630157.root";
  mydataset["MuonBeam_Run630158_32GeV"] = fileDirEOSDhcal+"DHCAL_630158.root";
  mydataset["MuonBeam_Run630159_32GeV"] = fileDirEOSDhcal+"DHCAL_630159.root";
  mydataset["MuonBeam_Run630160_32GeV"] = fileDirEOSDhcal+"DHCAL_630160.root";
  mydataset["MuonBeam_Run600006_32GeV"] = fileDirEOSDhcal+"DHCAL_600006.root";
  mydataset["MuonBeam_Run600007_32GeV"] = fileDirEOSDhcal+"DHCAL_600007.root";
  mydataset["MuonBeam_Run600008_32GeV"] = fileDirEOSDhcal+"DHCAL_600008.root";
  mydataset["MuonBeam_Run600009_32GeV"] = fileDirEOSDhcal+"DHCAL_600009.root";
  mydataset["MuonBeam_Run600010_32GeV"] = fileDirEOSDhcal+"DHCAL_600010.root";
  mydataset["MuonBeam_Run600011_32GeV"] = fileDirEOSDhcal+"DHCAL_600011.root";
  mydataset["MuonBeam_Run600013_32GeV"] = fileDirEOSDhcal+"DHCAL_600013.root";
  mydataset["MuonBeam_Run600015_32GeV"] = fileDirEOSDhcal+"DHCAL_600015.root";
  mydataset["MuonBeam_Run600016_32GeV"] = fileDirEOSDhcal+"DHCAL_600016.root";
  mydataset["MuonBeam_Run600018_32GeV"] = fileDirEOSDhcal+"DHCAL_600018.root";
  mydataset["MuonBeam_Run600021_32GeV"] = fileDirEOSDhcal+"DHCAL_600021.root";
  mydataset["MuonBeam_Run600022_32GeV"] = fileDirEOSDhcal+"DHCAL_600022.root";
  mydataset["MuonBeam_Run600024_32GeV"] = fileDirEOSDhcal+"DHCAL_600024.root";
  mydataset["MuonBeam_Run600025_32GeV"] = fileDirEOSDhcal+"DHCAL_600025.root";
  mydataset["MuonBeam_Run600026_32GeV"] = fileDirEOSDhcal+"DHCAL_600026.root";
  mydataset["MuonBeam_Run600111_32GeV"] = fileDirEOSDhcal+"DHCAL_600111.root";
  mydataset["MuonBeam_Run600113_32GeV"] = fileDirEOSDhcal+"DHCAL_600113.root";  
  mydataset["MuonBeam_Run600114_32GeV"] = fileDirEOSDhcal+"DHCAL_600114.root";  
  mydataset["MuonBeam_Run600116_32GeV"] = fileDirEOSDhcal+"DHCAL_600116.root";  
  mydataset["MuonBeam_Run600117_32GeV"] = fileDirEOSDhcal+"DHCAL_600117.root";  
  mydataset["MuonBeam_Run600121_32GeV"] = fileDirEOSDhcal+"DHCAL_600121.root";  
  mydataset["MuonBeam_Run600122_32GeV"] = fileDirEOSDhcal+"DHCAL_600122.root";  
  mydataset["MuonBeam_Run600124_32GeV"] = fileDirEOSDhcal+"DHCAL_600124.root";
  mydataset["MuonBeam_Run600125_32GeV"] = fileDirEOSDhcal+"DHCAL_600125.root";
  mydataset["MuonBeam_Run600127_32GeV"] = fileDirEOSDhcal+"DHCAL_600127.root";
  mydataset["MuonBeam_Run600131_32GeV"] = fileDirEOSDhcal+"DHCAL_600131.root";
  mydataset["MuonBeam_Run600133_32GeV"] = fileDirEOSDhcal+"DHCAL_600133.root";
  mydataset["MuonBeam_Run600134_32GeV"] = fileDirEOSDhcal+"DHCAL_600134.root";
  mydataset["MuonBeam_Run600135_32GeV"] = fileDirEOSDhcal+"DHCAL_600135.root";
/*
DHCAL_600008.root
DHCAL_600009.root
DHCAL_600010.root
DHCAL_600011.root
DHCAL_600013.root
DHCAL_600015.root
DHCAL_600016.root
DHCAL_600018.root
DHCAL_600021.root
DHCAL_600022.root
DHCAL_600024.root
DHCAL_600025.root
DHCAL_600026.root
DHCAL_600111.root
DHCAL_600113.root
DHCAL_600114.root
DHCAL_600116.root
DHCAL_600117.root
DHCAL_600121.root
DHCAL_600122.root
DHCAL_600124.root
DHCAL_600125.root
DHCAL_600127.root
DHCAL_600131.root
DHCAL_600133.root
DHCAL_600134.root
DHCAL_600135.root
DHCAL_600214.root
DHCAL_600215.root
DHCAL_600216.root
DHCAL_600217.root
DHCAL_610034.root
DHCAL_610047.root
DHCAL_610051.root
DHCAL_610055.root
DHCAL_610056.root
DHCAL_610057.root
DHCAL_610063.root
*/
 
}


