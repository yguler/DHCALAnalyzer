#include "histDHCAL.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cassert>
#include <vector>
#include <time.h>

#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TGraph2D.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "inputDHCAL.h"
#include "TreeAnalysis.h"
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include "TVirtualFitter.h"
#include "TSystem.h"
#include "TVector3.h"
using namespace ROOT::Math;

// g++ `root-config --cflags` -o DHCALAnalyzer TreeAnalysis.cc histDHCAL.cc inputDHCAL.cc anaDHCAL.cc  `root-config --glibs`
// ./DHCALAnalyzer pdgID211_beamMomentum100_listFTFP_BERT_EMM 1 0 
double P=28.28;
histoSet1::histoSet1(TFile * fout,string anaType,string cutname){

  dirname=anaType+"_set1_"+cutname;
  //fout->mkdir(dirname.c_str());
  //fout->cd(dirname.c_str());
  const int max_Layer = 37;

  char name[100], Title[100], name2[100], Title2[100];
  for(int i = 0; i<max_Layer+1; i++){
        string s_ee="hNRPCHits_Layer";
        string s_layer=to_string(i+1);
        string hname=s_ee+s_layer;
        string title = "distribution of RPC Hits at layer";
        string htitle = title+s_layer;
        histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),50,0.0,50);
        histo1D[hname]->SetXTitle("Numbe of Hits");
        histo1D[hname]->SetYTitle("# Entry");
        histo1D[hname]->GetYaxis()->SetTitleOffset(1.);
        histo1D[hname]->Sumw2();

        s_ee="hClusterxy_Layer";
        s_layer=to_string(i+1);
        hname=s_ee+s_layer;
        title = "Positions of Cluster at layer";
        htitle = title+s_layer;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);
        //histo2D[hname]->Sumw2();

        s_ee="hRPCHitMap_Layer";
        s_layer=to_string(i+1);
        hname=s_ee+s_layer;
        title = "Positions of RPC Hits at layer";
        htitle = title+s_layer;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);  

	s_ee="Valid_RPCClusterMap_Layer_";
	s_layer=to_string(i+1);
	hname=s_ee+s_layer;
        htitle = hname;
	histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);

        s_ee="Eff_RPCClusterMap_Layer_";
        s_layer=to_string(i+1);
        hname=s_ee+s_layer;
        htitle = hname;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);
 
        s_ee="Eff_RPCClusterNum_Layer_";
        s_layer=to_string(i+1);
        hname=s_ee+s_layer;
        htitle = hname;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);

        s_ee="Mult_HitMap_Layer_";
        s_layer=to_string(i+1);
        hname=s_ee+s_layer;
        htitle = hname;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);	

        s_ee="Eff_ExcludedRPCClusterMap_Layer_";
        s_layer=to_string(i+1);
        hname=s_ee+s_layer;
        htitle = hname;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),96,-0.5,95.5, 96, -0.5, 95.5);
        histo2D[hname]->SetXTitle("x");
        histo2D[hname]->SetYTitle("y");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);
}

  for(int i = 0; i<3; i++){
        string hname="hNHitsRPC_id"+to_string(i+1)+"_vs_Layer";
        string title = "Number of Hits vs Layer Numbers";
        string htitle = "RPC id "+to_string(i+1) + "  "+title;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),55,0.5,55.5, 30, -0.5, 29.5);
        histo2D[hname]->SetXTitle("Layer Numbers");
        histo2D[hname]->SetYTitle("Number of Hits");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);	
  }
        string hname="hNHits";
        string title = "distribution of Nhits";
        string htitle = title;
        histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),150,0.0,1500);
        histo1D[hname]->SetXTitle("Number of Hits");
        histo1D[hname]->SetYTitle("Event");
        histo1D[hname]->GetYaxis()->SetTitleOffset(1.);
        histo1D[hname]->Sumw2();

        hname="hNRPChits_vs_Layer";
        title = "Number of RPC Hits vs Layer Numbers";
        htitle = title;
        histo2D[hname]=new TH2D(hname.c_str(),htitle.c_str(),55,0.5,55.5, 30, -0.5, 29.5);
        histo2D[hname]->SetXTitle("Layer Numbers");
        histo2D[hname]->SetYTitle("Number of Hits");
        histo2D[hname]->GetYaxis()->SetTitleOffset(1.);
        //histo2D[hname]->Sumw2();

	histo1D["RPC0Efficiency"]=new TH1D("RPC0Efficiency","RPC0 Efficiency",38,0.5,38.5);
	histo1D["RPC0InEfficiency"]=new TH1D("RPC0InEfficiency","RPC0 InEfficiency",38,0.5,38.5);
        histo1D["RPC1Efficiency"]=new TH1D("RPC1Efficiency","RPC1 Efficiency",38,0.5,38.5);
        histo1D["RPC1InEfficiency"]=new TH1D("RPC1InEfficiency","RPC1 InEfficiency",38,0.5,38.5);
        histo1D["RPC2Efficiency"]=new TH1D("RPC2Efficiency","RPC2 Efficiency",38,0.5,38.5);
        histo1D["RPC2InEfficiency"]=new TH1D("RPC2InEfficiency","RPC2 InEfficiency",38,0.5,38.5);
    
	histo1D["RPC0ActiveLayerNum"]=new TH1D("RPC0ActiveLayerNum","RPC0 Active Layer Num",38,0.5,38.5);
        histo1D["RPC1ActiveLayerNum"]=new TH1D("RPC1ActiveLayerNum","RPC1 Active Layer Num",38,0.5,38.5);
        histo1D["RPC2ActiveLayerNum"]=new TH1D("RPC2ActiveLayerNum","RPC2 Active Layer Num",38,0.5,38.5);
	histo1D["RPC0PadMultiplicity"]=new TH1D("RPC0PadMultiplicity","RPC0 Pad Multiplicity",38,0.5,38.5);
	histo1D["RPC1PadMultiplicity"]=new TH1D("RPC1PadMultiplicity","RPC1 Pad Multiplicity",38,0.5,38.5);
	histo1D["RPC2PadMultiplicity"]=new TH1D("RPC2PadMultiplicity","RPC2 Pad Multiplicity",38,0.5,38.5);

        hname="ExcludedRPCCluster3D";
        htitle = hname;
        histo3D[hname]=new TH3D(hname.c_str(),htitle.c_str(),38,0.5,38.5,96,-0.5,95.5, 96,-0.5, 95.5);
        histo3D[hname]->SetXTitle("z");
        histo3D[hname]->SetYTitle("x");
	histo3D[hname]->SetYTitle("y");
        histo3D[hname]->GetYaxis()->SetTitleOffset(1.);

        histo1D["hDelta"] = new TH1D("hDelta","",10000,0,100);
        histo2D["hDeltaVsZ"]  = new TH2D("hDeltaVsZ","",50,0,50,1000,0,10);
        histo2D["hRes"]  = new TH2D("hRes","",100,-10,10,100,-10,10);

	histo1D["GoodBadHit"] = new TH1D("GoodBadHit","",3,0,3);
        histo1D["first3LayerCheck"] = new TH1D("first3LayerCheck","",3,0,3);
	histo1D["last3LayerCheck"] = new TH1D("last3LayerCheck","",3,0,3);
	histo1D["totalCut"] = new TH1D("totalCut","",3,0,3);
        histo1D["DeltaHit"] = new TH1D("DeltaHit","",3,0,3);
	TGraphErrors *gr1 = new TGraphErrors();
	//gr1->Write("Fit");
}

histoSet1::~histoSet1(){
}

void histoSet1::summary(string anaType, string mydataset){
    // std::cout<<" yg  "<<std::endl;
}

bool first = true;

struct SumDistance2 {
   // the TGraph is a data member of the object
   TGraph2D * fGraph;

   SumDistance2(TGraph2D * g) : fGraph(g) {}
   inputDHCAL input;
   // implementation of the function to be minimized
   double operator() (const double * par) {

      assert(fGraph    != 0);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      double * dist;
      //cout<<" npoints = "<<npoints<<endl;
      for (int i  = 0; i < npoints; ++i) {
         double d = input.distance2(x[i],y[i],z[i],par);
         sum += d;
	//dist[z]=d;
#ifdef DEBUG
        //if (first) 
	std::cout << "point " << i << "\t"
            << x[i] << "\t\t"
            << y[i] << "\t\t"
            << z[i] << "\t\t"
	    << d << "\t\t"
            << std::sqrt(d) << std::endl;
#endif
      }
      if (first)
         //std::cout << "Total Initial distance square = " << sum << std::endl;
      first = false;
      return sum;
   }

};

void histoSet1::process(TreeAnalysis & ev, double weight){

   // double wt=1.0;
   double wt=weight;
   uint32_t event;
   uint32_t run;
   double beamEnergy;
   double trueBeamEnergy;
   run = ev.RunNo;
   beamEnergy =ev.Energy;
   inputDHCAL input;

	string hname="hNHits";
	histo1D[hname]->Fill(ev.fNRPCHits);
	
	int NClHitsPerLayers[55]={0};
        int NRPCHitsPerLayers[55]={0}; 
	
	for (int i=0; i<ev.fNRPCHits; ++i) {
		double x=ev.RPCHits_x[i]; double y=ev.RPCHits_y[i]; int z=((int)ev.RPCHits_z[i]);
		//cout<< i<<" RPCHits_x = "<<x<<" RPCHits_y = "<<y<< " RPCHits_z = "<<z<<endl;
        	string s_ee="hRPCHitMap_Layer" +to_string(z+1);
        	string hname=s_ee;
		histo2D[hname]->Fill(x,y);
		NRPCHitsPerLayers[z]+=1;
	}
	int NLHits[55][3]={0};
	for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
		double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; int z1=((int)ev.RPCClusters_NN_z[i]);
			string s="hClusterxy_Layer"+to_string(z1+1);
			string hname=s;
			histo2D[hname]->Fill(x1,y1);
			NClHitsPerLayers[z1]+=1;
                	int rpcid=-1;
                	if((y1>3 && y1<=28) && (x1>3 && x1<93)) rpcid=0;
                	else if((y1>35 && y1<=60) && (x1>3 && x1<93)) rpcid=1;
                	else rpcid=2;	
			NLHits[z1][rpcid]+=1;
		//cout<< i<<" RPCClusters_NN_x = "<<x1<<" RPCClusters_NN_y = "<<y1<< " RPCClusters_NN_z = "<<z1<<endl;
	}
	for(int i1=0;i1<38;i1++){
                string s="hNRPCHits_Layer"+to_string(i1+1);
                string hname=s;
                histo1D[hname]->Fill(NRPCHitsPerLayers[i1]);
        	hname="hNRPChits_vs_Layer";
		histo2D[hname]->Fill(i1+1,NRPCHitsPerLayers[i1]);
		for (int j=0; j<3; j++){
			hname="hNHitsRPC_id"+to_string(j+1)+"_vs_Layer";
			histo2D[hname]->Fill(i1+1,NLHits[i1][j]);
		}
	}
	int hitCheck=0;
	int zL[38]={0};
	int zL2[38]={0};
	for (int j=0; j<38; ++j) {zL[j]=0; zL2[j]=0;}
		
		for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
			double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i];
			int z1=((int)ev.RPCClusters_NN_z[i]);
			if (input.LayerCheck(ev,z1)){ 
				hitCheck = input.RPCClusterHitSizeValid(ev, x1, y1, z1);
				if (hitCheck>3) zL[z1]+=1;
				zL2[z1]=hitCheck;
				hitCheck=0;
			}
		}
	bool hitValid=true;
	bool deltaHit=false;
	for (int j=0; j<38; ++j){
		if(zL[j]>0 && zL[j+1]>0) hitValid=false;
		if(zL2[j]>0 && zL2[j+1]>zL2[j] && zL2[j+2]> zL2[j+1] && zL2[j+3]> zL2[j+2]) deltaHit=true;
	}
	bool first3Layer = false; bool last3Layer = false; 
	first3Layer = (input.LayerCheck(ev,0) && input.LayerCheck(ev,1) && input.LayerCheck(ev,2));
	last3Layer  = (input.LayerCheck(ev,35) && input.LayerCheck(ev,36) && input.LayerCheck(ev,37));
	histo1D["GoodBadHit"]->Fill(hitValid);
	histo1D["first3LayerCheck"]->Fill(first3Layer);
	histo1D["last3LayerCheck"]->Fill(last3Layer);
	bool totalcut = false;
	totalcut = (hitValid && first3Layer && last3Layer);
	histo1D["totalCut"]->Fill(totalcut);
	histo1D["DeltaHit"]->Fill(!deltaHit);
	if (hitValid){
				TGraph *gr0 = new TGraph();
				TGraph *gr1 = new TGraph();
				TGraph *gr2 = new TGraph();
				TGraph2D *gr3d0 = new TGraph2D();
				TGraph2D *gr3d1 = new TGraph2D();
				TGraph2D *gr3d2 = new TGraph2D();
				double rpc0ActiveLayer[38][2]={0}; double rpc1ActiveLayer[38][2]={0}; double rpc2ActiveLayer[38][2]={0};
				double f0dx=0.0; double f0dy=0.0; double l0dx=0.0; double l0dy=0.0; int rpcid = -1; int i0 =0; int i1= 0; int i2=0;
				double f1dx=0.0; double f1dy=0.0; double l1dx=0.0; double l1dy=0.0;
				double f2dx=0.0; double f2dy=0.0; double l2dx=0.0; double l2dy=0.0;
				bool rpc0 = false; bool rpc1 = false; bool rpc2 = false; 
				int RPC0 = -1; int RPC1 = -1; int RPC2 = -1;
				double rpc0mapX[38]={0}; double rpc0mapY[38]={0}; double rpc1mapX[38]={0}; double rpc1mapY[38]={0};
				double rpc2mapX[38]={0}; double rpc2mapY[38]={0};
                                if (input.ActiveLayerCheck(ev,0)==0 && input.ActiveLayerCheck(ev,1)==0 && input.ActiveLayerCheck(ev,2)==0 &&
                                   /*input.ActiveLayerCheck(ev,35)==0 &&*/ input.ActiveLayerCheck(ev,36)==0 && input.ActiveLayerCheck(ev,37)==0) RPC0=0;
                                if (input.ActiveLayerCheck(ev,0)==1 && input.ActiveLayerCheck(ev,1)==1 && input.ActiveLayerCheck(ev,2)==1 &&
                                   input.ActiveLayerCheck(ev,35)==1 && input.ActiveLayerCheck(ev,36)==1 && input.ActiveLayerCheck(ev,37)==1) RPC1=1;
                                if (input.ActiveLayerCheck(ev,0)==2 && input.ActiveLayerCheck(ev,1)==2 && input.ActiveLayerCheck(ev,2)==2 &&
                                   input.ActiveLayerCheck(ev,35)==2 && input.ActiveLayerCheck(ev,36)==2 && input.ActiveLayerCheck(ev,37)==2) RPC2=2;
        			for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                			double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i]; 
					int z1=((int)ev.RPCClusters_NN_z[i]);
					if((y1>3 && y1<=28) && (x1>3 && x1<93) && RPC0==0){
						if (z1==12) continue;
						if(input.ClusterDistValid(ev,z1,i,RPC0)){
							gr3d0->SetPoint(i0,double(x1), y1, z1);
							++i0; rpc0 = true;
							rpc0ActiveLayer[z1][0]=x1;
							rpc0ActiveLayer[z1][1]=y1;
					//cout<< i0<<" rpc0 "<<" RPCClusters x = "<<x1<<" RPCClusters y = "<<y1<< " RPCClusters z = "<<z1<<endl;
						}
					}
                                        if((y1>35 && y1<=60) && (x1>3 && x1<93) && RPC1==1){
						if(input.ClusterDistValid(ev,z1,i,RPC1)){
							gr3d1->SetPoint(i1,double(x1), y1, z1);
                                                	++i1; rpc1 = true;
							rpc1ActiveLayer[z1][0]=x1;
							rpc1ActiveLayer[z1][1]=y1;
					//cout<< i1<<" rpc1 "<<" RPCClusters x = "<<x1<<" RPCClusters y = "<<y1<< " RPCClusters z = "<<z1<<endl;
						}
                                        }
                                        if((y1>65 && y1<93) && (x1>3 && x1<93) && RPC2==2){
						if(input.ClusterDistValid(ev,z1,i,RPC2)){
							gr3d2->SetPoint(i2,double(x1), y1, z1);
                                                	++i2; rpc2 = true;
							rpc2ActiveLayer[z1][0]=x1;
							rpc2ActiveLayer[z1][1]=x1;
						//cout<< i2<<" rpc2 "<<" RPCClusters x = "<<x1<<" RPCClusters y = "<<y1<< " RPCClusters z = "<<z1<<endl;
						}
                                        }
					//cout<< i<<" RPCClusters x = "<<x1<<" RPCClusters y = "<<y1<< " RPCClusters z = "<<z1<<endl;
			    	}
				if (rpc0){
                                ROOT::Fit::Fitter  fitter;
                                SumDistance2 sdist(gr3d0);
                                #ifdef __CINT__
                                  ROOT::Math::Functor fcn(&sdist, 4,"SumDistance2");
                                #else
                                  ROOT::Math::Functor fcn(sdist,4);
                                #endif
                                double pStart[4] = {1,1,1,1};
                                fitter.SetFCN(fcn,pStart,0,true);
                                // set step sizes different than default ones (0.3 times parameter values)
                                for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
                                           bool ok = fitter.FitFCN();
                                           if (!ok) {
                                              Error("line3Dfit","Line3D Fit failed");
                                           }
                                const ROOT::Fit::FitResult & result = fitter.Result();
                                const double * parFit = result.GetParams();
  				double * x = gr3d0->GetX();
  				double * y = gr3d0->GetY();
  				double * z = gr3d0->GetZ();
  				int npoints = gr3d0->GetN();
				double dist[npoints]={-1.};
				double totalDistance=0;
				bool goodFit = true;
  				for (int i  = 0; i < npoints; ++i) {
					int z1=z[i];
					if(z1==12) continue;
    					XYZVector xp(x[i],y[i],z[i]);
    					XYZVector x0(parFit[0], parFit[2], 0. );
    					XYZVector x1(parFit[0] + parFit[1], parFit[2] + parFit[3], 1. );
    					XYZVector u = (x1-x0).Unit();
    					XYZVector cross = ((xp-x0).Cross(u));
    					double d2 = cross.Mag2();
    					Double_t d = TMath::Sqrt(d2);
					dist[i]=d;
					totalDistance+=d;
					if (z[i]==0) {f0dx=x[i]; f0dy=y[i];}
					if (z[i]==37) {l0dx=x[i]; l0dy=y[i];}
//cout<<i <<" x = "<<x[i]<<" y = "<<y[i]<<" z = "<<z[i]<<" parFit0 = "<<parFit[0]<<" parFit2 = "<<parFit[2]<<" d2 = "<<d2<<" sqrtd2 = "<<d<<endl;
    					histo1D["hDelta"]->Fill(d);
    					histo2D["hDeltaVsZ"]->Fill(z[i],d);
    					histo2D["hRes"]->Fill(cross.Y(),cross.Z());
				}
				double dr = sqrt(pow(f0dx-l0dx,2) +pow(f0dy-l0dy,2));
				int effTrue[38]={0}; int effFalse[38]={0}; bool eff=false;
				double mapX[38]={0}; double mapY[38]={0}; int Nhit[38]={0};
				double LXYdistNhits[38][4]={-1.};
				bool multiLayer = false; int ind=-1; 		
				for (int i  = 0; i < npoints; ++i) {
					int z1 = z[i];
					ind=i;
					for (int jj  = 0; jj < npoints; ++jj) {
						int z2 = z[jj];
						if (z1==z2 && i!=jj) {
							ind = ( dist[i]<dist[jj])? i:jj; break;
						}
					}
							LXYdistNhits[z1][0]=x[ind];
							LXYdistNhits[z1][1]=y[ind];
							LXYdistNhits[z1][2]=dist[ind];
							LXYdistNhits[z1][3]=input.RPCClusterHitSizeValid(ev, x[ind], y[ind], z1);
				//cout<<i<<" ind = "<<ind<<" z1 = "<<z1<<" x1 = "<<x[ind]<<" y1 = "<<y[ind]<<endl;
				}
                		//TVector3 InitVec(parFit[0],parFit[2],0.);
                		//TVector3 DirVec(parFit[1],parFit[3],1.);
                		//TVector3 UnitDirVec=DirVec.Unit();
                                for (int j=0; j<ev.fNRPCClusters_NN; ++j) {
					if (!goodFit) break;
					double dr = sqrt(pow(f0dx-l0dx,2) +pow(f0dy-l0dy,2));
                                        double x1=ev.RPCClusters_NN_x[j]; double y1=ev.RPCClusters_NN_y[j];
                                        int z1=((int)ev.RPCClusters_NN_z[j]);
					if (z1==12) continue;
                                        if ((y1>3 && y1<=28) && (x1>3 && x1<93)){
                                                if ((result.Chi2()/npoints) < 1 && (dr/npoints)<0.5 ){
                                                        for (int i  = 0; i < 38; ++i) {
                                                                if (x1==LXYdistNhits[i][0] && y1==LXYdistNhits[i][1] && z1==i){
									histo1D["RPC0ActiveLayerNum"]->Fill(z1+1,1);
									string hname="Valid_RPCClusterMap_Layer_"+to_string(z1+1);
									histo2D[hname]->Fill(x1,y1);
									
									if(LXYdistNhits[i][2]<=2 && LXYdistNhits[i][2]!=-1){
										string hname ="Eff_RPCClusterMap_Layer_" + to_string(z1+1);
										string hname2="Eff_RPCClusterNum_Layer_"+to_string(z1+1);
										string hname3="Mult_HitMap_Layer_"+to_string(z1+1);
                                                				histo1D["RPC0Efficiency"]->Fill(z1+1,1);
                                                				histo2D[hname]->Fill(x1,y1);
                                                				histo2D[hname2]->Fill(x1,y1,1);
                                                				histo2D[hname3]->Fill(x1,y1,LXYdistNhits[i][3]);
                                                				histo1D["RPC0PadMultiplicity"]->Fill(z1+1,LXYdistNhits[i][3]);	
					//cout<<j<<" efficient = "<<" dist = "<<LXYdistNhits[i][2]<<" Layer = "<<z1<<" Nhits = "<<LXYdistNhits[i][3]<<endl;
									}
								}
                                                	}
						}
						for (int i  = 0; i < 38; ++i) {
							if((result.Chi2()/npoints)>1 || (dr/npoints)>0.5 || LXYdistNhits[i][2]>2) { 
								if (z1==i){
									histo1D["RPC0InEfficiency"]->Fill(i+1);
								//cout<<i<<" Inefficient = "<<" dist = "<<LXYdistNhits[i][2]<<" Layer = "<<z1<<endl;
								}
							}
						}
					}
                                }
				}// rpc0 ended

				if(rpc1){
                                ROOT::Fit::Fitter  fitter;
                                SumDistance2 sdist(gr3d1);
                                #ifdef __CINT__
                                  ROOT::Math::Functor fcn(&sdist, 4,"SumDistance2");
                                #else
                                  ROOT::Math::Functor fcn(sdist,4);
                                #endif
                                double pStart[4] = {1,1,1,1};
                                fitter.SetFCN(fcn,pStart,0,true);
                                // set step sizes different than default ones (0.3 times parameter values)
                                for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
                                           bool ok = fitter.FitFCN();
                                           if (!ok) {
                                              Error("line3Dfit","Line3D Fit failed");
                                           }
                                const ROOT::Fit::FitResult & result = fitter.Result();
                                const double * parFit = result.GetParams();
                                double * x = gr3d1->GetX();
                                double * y = gr3d1->GetY();
                                double * z = gr3d1->GetZ();
                                int npoints = gr3d1->GetN();
                                double dist[npoints]={-1.};
				bool goodFit = true;
                                for (int i  = 0; i < npoints; ++i) {
					int z1=z[i];
                                        XYZVector xp(x[i],y[i],z[i]);
                                        XYZVector x0(parFit[0], parFit[2], 0. );
                                        XYZVector x1(parFit[0] + parFit[1], parFit[2] + parFit[3], 1. );
                                        XYZVector u = (x1-x0).Unit();
                                        XYZVector cross = ((xp-x0).Cross(u));
                                        double d2 = cross.Mag2();
                                        Double_t d = TMath::Sqrt(d2);
                                        dist[i]=d;
					if (z[i]==0) {f1dx=x[i]; f1dy=y[i];}
					if (z[i]==37) {l1dx=x[i]; l1dy=y[i];}
//cout<<i <<" x = "<<x[i]<<" y = "<<y[i]<<" z = "<<z[i]<<" parFit0 = "<<parFit[0]<<" parFit2 = "<<parFit[2]<<" d2 = "<<d2<<" sqrtd2 = "<<d<<endl;
                                }
                                int effTrue[38]={0}; int effFalse[38]={0}; bool eff=false;
				double mapX[38]={0}; double mapY[38]={0}; int Nhit[38]={0};
				int ind=-1; int t=0; int zl[100]={0};
				bool single=true;
                                double LXYdistNhits[38][4]={-1.};
                                bool multiLayer = false;
                                for (int i  = 0; i < npoints; ++i) {
                                        int z1 = z[i];
                                        ind=i;
                                        for (int jj  = 0; jj < npoints; ++jj) {
                                                int z2 = z[jj];
                                                if (z1==z2 && i!=jj) {
                                                        ind = ( dist[i]<dist[jj])? i:jj; break;
                                                }
                                        }
                                                        LXYdistNhits[z1][0]=x[ind];
                                                        LXYdistNhits[z1][1]=y[ind];
                                                        LXYdistNhits[z1][2]=dist[ind];
                                                        LXYdistNhits[z1][3]=input.RPCClusterHitSizeValid(ev, x[ind], y[ind], z1);
					//cout<<i<<" ind = "<<ind<<" z1 = "<<z1<<" x1 = "<<x[ind]<<" y1 = "<<y[ind]<<endl;
                                }
                                for (int j=0; j<ev.fNRPCClusters_NN; ++j) {
					double dr = sqrt(pow(f1dx-l1dx,2) +pow(f1dy-l1dy,2));
                                        double x1=ev.RPCClusters_NN_x[j]; double y1=ev.RPCClusters_NN_y[j];
                                        int z1=((int)ev.RPCClusters_NN_z[j]);
					int tut=-1;
                                        if ((y1>35 && y1<=60) && (x1>3 && x1<93)){
                                                if ((result.Chi2()/npoints) < 1 && (dr/npoints)<0.5 ){
                                                        for (int i  = 0; i < 38; ++i) {
                                                                if (x1==LXYdistNhits[i][0] && y1==LXYdistNhits[i][1] && z1==i){
                                                                        histo1D["RPC1ActiveLayerNum"]->Fill(z1+1,1);
                                                                        string hname="Valid_RPCClusterMap_Layer_"+to_string(z1+1);
                                                                        histo2D[hname]->Fill(LXYdistNhits[i][0],LXYdistNhits[i][1]);
                                                                        if(LXYdistNhits[i][2]<=2 && LXYdistNhits[i][2]!=-1){
                                                                                string hname ="Eff_RPCClusterMap_Layer_" + to_string(z1+1);
                                                                                string hname2="Eff_RPCClusterNum_Layer_"+to_string(z1+1);
                                                                                string hname3="Mult_HitMap_Layer_"+to_string(z1+1);
                                                                                histo1D["RPC1Efficiency"]->Fill(z1+1,1);
                                                                                histo2D[hname]->Fill(x1,y1);
                                                                                histo2D[hname2]->Fill(x1,y1,1);
                                                                                histo2D[hname3]->Fill(x1,y1,LXYdistNhits[i][3]);
                                                                                histo1D["RPC1PadMultiplicity"]->Fill(z1+1,LXYdistNhits[i][3]);
                                        //cout<<i<<" efficient  "<<" dist = "<<LXYdistNhits[i][2]<<" Layer = "<<z1<<" Nhits = "<<LXYdistNhits[i][3]<<endl;
                                                                        }
                                                                }
                                                        }
                                                }
                                        	for (int i  = 0; i < 38; ++i) {
                                                	if((result.Chi2()/npoints)>1 || (dr/npoints)>0.5 || LXYdistNhits[i][2]>2) {
                                                        	if ( z1==i){
                                                                	effFalse[i]=1; eff=false;
									histo1D["RPC1InEfficiency"]->Fill(i+1);
                                                          //cout<<i<<" Ineff = "<<effFalse[i]<<" dist = "<<LXYdistNhits[i][2]<<" Layer = "<<z1<<endl;
                                                        	}
                                                	}
                                        	}
                                        }
                                }
                                }	// rpc1
                                if(rpc2){
                                ROOT::Fit::Fitter  fitter;
                                SumDistance2 sdist(gr3d2);
                                #ifdef __CINT__
                                  ROOT::Math::Functor fcn(&sdist, 4,"SumDistance2");
                                #else
                                  ROOT::Math::Functor fcn(sdist,4);
                                #endif
                                double pStart[4] = {1,1,1,1};
                                fitter.SetFCN(fcn,pStart,0,true);
                                // set step sizes different than default ones (0.3 times parameter values)
                                for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
                                           bool ok = fitter.FitFCN();
                                           if (!ok) {
                                              Error("line3Dfit","Line3D Fit failed");
                                           }
                                const ROOT::Fit::FitResult & result = fitter.Result();
                                const double * parFit = result.GetParams();
                                double * x = gr3d2->GetX();
                                double * y = gr3d2->GetY();
                                double * z = gr3d2->GetZ();
                                int npoints = gr3d2->GetN();
                                double dist[npoints]={-1.};
				bool goodFit = true;
                                for (int i  = 0; i < npoints; ++i) {
					int z1=z[i];
                                        XYZVector xp(x[i],y[i],z[i]);
                                        XYZVector x0(parFit[0], parFit[2], 0. );
                                        XYZVector x1(parFit[0] + parFit[1], parFit[2] + parFit[3], 1. );
                                        XYZVector u = (x1-x0).Unit();
                                        XYZVector cross = ((xp-x0).Cross(u));
                                        double d2 = cross.Mag2();
                                        Double_t d = TMath::Sqrt(d2);
                                        dist[i]=d;
					if (z[i]==0) {f2dx=x[i]; f2dy=y[i];}
					if (z[i]==37) {l2dx=x[i]; l2dy=y[i];}
//cout<<i <<" x = "<<x[i]<<" y = "<<y[i]<<" z = "<<z[i]<<" parFit0 = "<<parFit[0]<<" parFit2 = "<<parFit[2]<<" d2 = "<<d2<<" sqrtd2 = "<<d<<endl;
                                }
                                        int effTrue[38]={0}; int effFalse[38]={0}; bool eff=false;
					double mapX[38]={0}; double mapY[38]={0}; int Nhit[38]={0};
                                double LXYdistNhits[38][4]={-1.};
                                bool multiLayer = false; int ind=-1;
                                for (int i  = 0; i < npoints; ++i) {
                                        int z1 = z[i];
                                        ind=i;
                                        for (int jj  = 0; jj < npoints; ++jj) {
                                                int z2 = z[jj];
                                                if (z1==z2 && i!=jj) {
                                                        cout<<" multi layer = "<<z1<<endl;
                                                        ind = ( dist[i]<dist[jj])? i:jj; break;
                                                }
                                        }
                                                        LXYdistNhits[z1][0]=x[ind];
                                                        LXYdistNhits[z1][1]=y[ind];
                                                        LXYdistNhits[z1][2]=dist[ind];
                                                        LXYdistNhits[z1][3]=input.RPCClusterHitSizeValid(ev, x[ind], y[ind], z1);
						//cout<<i<<" ind = "<<ind<<" z1 = "<<z1<<" x1 = "<<x[ind]<<" y1 = "<<y[ind]<<endl;
                                }
                                for (int j=0; j<ev.fNRPCClusters_NN; ++j) {
					if (!goodFit) break;
					double dr = sqrt(pow(f2dx-l2dx,2) +pow(f2dy-l2dy,2));
                                        double x1=ev.RPCClusters_NN_x[j]; double y1=ev.RPCClusters_NN_y[j];
                                        int z1=((int)ev.RPCClusters_NN_z[j]);
                                        if ((y1>65 && y1<93) && (x1>3 && x1<93)){
                                                if ((result.Chi2()/npoints) < 1 && (dr/npoints)<0.5 ){
                                                        for (int i  = 0; i < 38; ++i) {
                                                                if (x1==LXYdistNhits[i][0] && y1==LXYdistNhits[i][1] && z1==i){
                                                                        histo1D["RPC2ActiveLayerNum"]->Fill(z1+1,1);
                                                                        string hname="Valid_RPCClusterMap_Layer_"+to_string(z1+1);
                                                                        histo2D[hname]->Fill(LXYdistNhits[i][0],LXYdistNhits[i][1]);
                                                                        if(LXYdistNhits[i][2]<=2 && LXYdistNhits[i][2]!=-1){
                                                                                string hname ="Eff_RPCClusterMap_Layer_" + to_string(z1+1);
                                                                                string hname2="Eff_RPCClusterNum_Layer_"+to_string(z1+1);
                                                                                string hname3="Mult_HitMap_Layer_"+to_string(z1+1);
                                                                                histo1D["RPC2Efficiency"]->Fill(z1+1,1);
                                                                                histo2D[hname]->Fill(x1,y1);
                                                                                histo2D[hname2]->Fill(x1,y1,1);
                                                                                histo2D[hname3]->Fill(x1,y1,LXYdistNhits[i][3]);
                                                                                histo1D["RPC2PadMultiplicity"]->Fill(z1+1,LXYdistNhits[i][3]);
                                //cout<<i<<" eff = "<<effTrue[z1]<<" dist = "<<LXYdistNhits[i][2]<<" Layer = "<<z1<<" Nhits = "<<LXYdistNhits[i][3]<<endl;
                                                                        }
                                                                }
                                                        }
                                                }
                                                for (int i  = 0; i < 38; ++i) {
                                                        if((result.Chi2()/npoints)>1 || (dr/npoints)>0.5 || LXYdistNhits[i][2]>2) {
                                                                if ( z1==i){
                                                                effFalse[i]=1; eff=false;
								histo1D["RPC2InEfficiency"]->Fill(i+1);
                                                                //cout<<i<<" Ineff = "<<effFalse[i]<<" dist = "<<LXYdistNhits[i][2]<<" Layer = "<<z1<<endl;
                                                                }
                                                        }
                                                }
                                        }
                                }
                                } 	//rpc2	
        }

	if (!hitValid || !input.LayerCheck(ev,0) || !input.LayerCheck(ev,1) || !input.LayerCheck(ev,2) || 
			!input.LayerCheck(ev,(input.LayerMax(ev)-2)) || !input.LayerCheck(ev,(input.LayerMax(ev)-1)) || 
			!input.LayerCheck(ev,input.LayerMax(ev)) /*|| !input.ClusterDistValid(ev)*/ || input.RPCClusterSizeValid(ev)){
			
                                for (int i=0; i<ev.fNRPCClusters_NN; ++i) {
                                        double x1=ev.RPCClusters_NN_x[i]; double y1=ev.RPCClusters_NN_y[i];
                                        int z1=((int)ev.RPCClusters_NN_z[i]);
                                        string s_ee="Eff_ExcludedRPCClusterMap_Layer_" + to_string(z1+1);
                                        string hname=s_ee;
                                        histo2D[hname]->Fill(x1,y1);
					hname="ExcludedRPCCluster3D";
					histo3D[hname]->Fill(z1,x1,y1);
				//cout<<"Delta Hit DHCAL"<<"   "<<i+1<<"   "<<"energy" <<"   "<<"5"<<"   "<<"xyz"<<"   "<<x1<<"   "<<y1<<"   "<<z1<<endl;
				}
				//cout<<"DHCAL"<<"   "<<"-999"<<"   "<<"done"<<endl;
	}
	if (deltaHit){
        	for (int i=0; i<ev.fNRPCHits; ++i) {
                	double x=ev.RPCHits_x[i]; double y=ev.RPCHits_y[i]; int z=((int)ev.RPCHits_z[i]);
			//cout<<"DHCALdetaHit"<<"   "<<i+1<<"   "<<"energy" <<"   "<<"5"<<"   "<<"xyz"<<"   "<<x<<"   "<<y<<"   "<<z<<endl;
        	}
		//cout<<"DHCAL"<<"   "<<"-999"<<"   "<<"done"<<endl;
	}
	//std::cout<<" yg  end"<<std::endl;

   return;
}

