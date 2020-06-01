//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 22 19:49:33 2020 by ROOT version 5.34/38
// from TTree DHCAL/DHCAL
// found on file: DHCAL_650357.root
//////////////////////////////////////////////////////////

#ifndef TreeAnalysis_h
#define TreeAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "vector"
// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxRPCHits = 9926;
   const Int_t kMaxRPCClusters_NN = 2061;
   const Int_t kMaxTCMTHits = 1;

class TreeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //DHCALEvent      *DHCALEvent;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           fNRPCHits;
   Int_t           fNRPCClusters_NN;
   Int_t           fNTCMTHits;
   Int_t           RunNo;
   Int_t           Energy;
   Int_t           Cerenkov;
   Int_t           EventTime;
   Int_t           InteractionLayer;
   Int_t           NOutOfTimeHits;
   Int_t           NoiseIndex;
   Int_t           MaxBoardOcc;
   Int_t           MaxBoardFrameOcc;
   Int_t           ASICOccFlag;
   Int_t           ASICOccFlagID;
   Int_t           NHitsLateral2cm;
   Int_t           NHitsLateral2cmHV;
   Int_t           RPCHits_;
   UInt_t          RPCHits_fUniqueID[kMaxRPCHits];   //[RPCHits_]
   UInt_t          RPCHits_fBits[kMaxRPCHits];   //[RPCHits_]
   Int_t           RPCHits_x[kMaxRPCHits];   //[RPCHits_]
   Int_t           RPCHits_y[kMaxRPCHits];   //[RPCHits_]
   Int_t           RPCHits_z[kMaxRPCHits];   //[RPCHits_]
   Int_t           RPCHits_t[kMaxRPCHits];   //[RPCHits_]
   Int_t           RPCHits_clid_NN[kMaxRPCHits];   //[RPCHits_]
   Int_t           RPCClusters_NN_;
   UInt_t          RPCClusters_NN_fUniqueID[kMaxRPCClusters_NN];   //[RPCClusters_NN_]
   UInt_t          RPCClusters_NN_fBits[kMaxRPCClusters_NN];   //[RPCClusters_NN_]
   Double_t        RPCClusters_NN_x[kMaxRPCClusters_NN];   //[RPCClusters_NN_]
   Double_t        RPCClusters_NN_y[kMaxRPCClusters_NN];   //[RPCClusters_NN_]
   Double_t        RPCClusters_NN_z[kMaxRPCClusters_NN];   //[RPCClusters_NN_]
   //vector<int>     RPCClusters_NN_idlist[kMaxRPCClusters_NN];
   Bool_t          RPCClusters_NN_valid[kMaxRPCClusters_NN];   //[RPCClusters_NN_]
   Int_t           TCMTHits_;
   UInt_t          TCMTHits_fUniqueID[kMaxTCMTHits];   //[TCMTHits_]
   UInt_t          TCMTHits_fBits[kMaxTCMTHits];   //[TCMTHits_]
   Int_t           TCMTHits_x[kMaxTCMTHits];   //[TCMTHits_]
   Int_t           TCMTHits_y[kMaxTCMTHits];   //[TCMTHits_]
   Int_t           TCMTHits_z[kMaxTCMTHits];   //[TCMTHits_]
   Double_t        TCMTHits_E[kMaxTCMTHits];   //[TCMTHits_]

   // List of branches
   TBranch        *b_DHCALEvent_fUniqueID;   //!
   TBranch        *b_DHCALEvent_fBits;   //!
   TBranch        *b_DHCALEvent_fNRPCHits;   //!
   TBranch        *b_DHCALEvent_fNRPCClusters_NN;   //!
   TBranch        *b_DHCALEvent_fNTCMTHits;   //!
   TBranch        *b_DHCALEvent_RunNo;   //!
   TBranch        *b_DHCALEvent_Energy;   //!
   TBranch        *b_DHCALEvent_Cerenkov;   //!
   TBranch        *b_DHCALEvent_EventTime;   //!
   TBranch        *b_DHCALEvent_InteractionLayer;   //!
   TBranch        *b_DHCALEvent_NOutOfTimeHits;   //!
   TBranch        *b_DHCALEvent_NoiseIndex;   //!
   TBranch        *b_DHCALEvent_MaxBoardOcc;   //!
   TBranch        *b_DHCALEvent_MaxBoardFrameOcc;   //!
   TBranch        *b_DHCALEvent_ASICOccFlag;   //!
   TBranch        *b_DHCALEvent_ASICOccFlagID;   //!
   TBranch        *b_DHCALEvent_NHitsLateral2cm;   //!
   TBranch        *b_DHCALEvent_NHitsLateral2cmHV;   //!
   TBranch        *b_DHCALEvent_RPCHits_;   //!
   TBranch        *b_RPCHits_fUniqueID;   //!
   TBranch        *b_RPCHits_fBits;   //!
   TBranch        *b_RPCHits_x;   //!
   TBranch        *b_RPCHits_y;   //!
   TBranch        *b_RPCHits_z;   //!
   TBranch        *b_RPCHits_t;   //!
   TBranch        *b_RPCHits_clid_NN;   //!
   TBranch        *b_DHCALEvent_RPCClusters_NN_;   //!
   TBranch        *b_RPCClusters_NN_fUniqueID;   //!
   TBranch        *b_RPCClusters_NN_fBits;   //!
   TBranch        *b_RPCClusters_NN_x;   //!
   TBranch        *b_RPCClusters_NN_y;   //!
   TBranch        *b_RPCClusters_NN_z;   //!
//   TBranch        *b_RPCClusters_NN_idlist;   //!
   TBranch        *b_RPCClusters_NN_valid;   //!
   TBranch        *b_DHCALEvent_TCMTHits_;   //!
   TBranch        *b_TCMTHits_fUniqueID;   //!
   TBranch        *b_TCMTHits_fBits;   //!
   TBranch        *b_TCMTHits_x;   //!
   TBranch        *b_TCMTHits_y;   //!
   TBranch        *b_TCMTHits_z;   //!
   TBranch        *b_TCMTHits_E;   //!

   TreeAnalysis(TTree *tree=0);
   virtual ~TreeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeAnalysis_cxx
TreeAnalysis::TreeAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DHCAL_600006.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/user/yalcin/DHCAL/data/DHCAL_600006.root");
      }
      f->GetObject("DHCAL",tree);

   }
   Init(tree);
}

TreeAnalysis::~TreeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_DHCALEvent_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_DHCALEvent_fBits);
   fChain->SetBranchAddress("fNRPCHits", &fNRPCHits, &b_DHCALEvent_fNRPCHits);
   fChain->SetBranchAddress("fNRPCClusters_NN", &fNRPCClusters_NN, &b_DHCALEvent_fNRPCClusters_NN);
   fChain->SetBranchAddress("fNTCMTHits", &fNTCMTHits, &b_DHCALEvent_fNTCMTHits);
   fChain->SetBranchAddress("RunNo", &RunNo, &b_DHCALEvent_RunNo);
   fChain->SetBranchAddress("Energy", &Energy, &b_DHCALEvent_Energy);
   fChain->SetBranchAddress("Cerenkov", &Cerenkov, &b_DHCALEvent_Cerenkov);
   fChain->SetBranchAddress("EventTime", &EventTime, &b_DHCALEvent_EventTime);
   fChain->SetBranchAddress("InteractionLayer", &InteractionLayer, &b_DHCALEvent_InteractionLayer);
   fChain->SetBranchAddress("NOutOfTimeHits", &NOutOfTimeHits, &b_DHCALEvent_NOutOfTimeHits);
   fChain->SetBranchAddress("NoiseIndex", &NoiseIndex, &b_DHCALEvent_NoiseIndex);
   fChain->SetBranchAddress("MaxBoardOcc", &MaxBoardOcc, &b_DHCALEvent_MaxBoardOcc);
   fChain->SetBranchAddress("MaxBoardFrameOcc", &MaxBoardFrameOcc, &b_DHCALEvent_MaxBoardFrameOcc);
   fChain->SetBranchAddress("ASICOccFlag", &ASICOccFlag, &b_DHCALEvent_ASICOccFlag);
   fChain->SetBranchAddress("ASICOccFlagID", &ASICOccFlagID, &b_DHCALEvent_ASICOccFlagID);
   fChain->SetBranchAddress("NHitsLateral2cm", &NHitsLateral2cm, &b_DHCALEvent_NHitsLateral2cm);
   fChain->SetBranchAddress("NHitsLateral2cmHV", &NHitsLateral2cmHV, &b_DHCALEvent_NHitsLateral2cmHV);
   fChain->SetBranchAddress("RPCHits", &RPCHits_, &b_DHCALEvent_RPCHits_);
   fChain->SetBranchAddress("RPCHits.fUniqueID", RPCHits_fUniqueID, &b_RPCHits_fUniqueID);
   fChain->SetBranchAddress("RPCHits.fBits", RPCHits_fBits, &b_RPCHits_fBits);
   fChain->SetBranchAddress("RPCHits.x", RPCHits_x, &b_RPCHits_x);
   fChain->SetBranchAddress("RPCHits.y", RPCHits_y, &b_RPCHits_y);
   fChain->SetBranchAddress("RPCHits.z", RPCHits_z, &b_RPCHits_z);
   fChain->SetBranchAddress("RPCHits.t", RPCHits_t, &b_RPCHits_t);
   fChain->SetBranchAddress("RPCHits.clid_NN", RPCHits_clid_NN, &b_RPCHits_clid_NN);
   fChain->SetBranchAddress("RPCClusters_NN", &RPCClusters_NN_, &b_DHCALEvent_RPCClusters_NN_);
   fChain->SetBranchAddress("RPCClusters_NN.fUniqueID", RPCClusters_NN_fUniqueID, &b_RPCClusters_NN_fUniqueID);
   fChain->SetBranchAddress("RPCClusters_NN.fBits", RPCClusters_NN_fBits, &b_RPCClusters_NN_fBits);
   fChain->SetBranchAddress("RPCClusters_NN.x", RPCClusters_NN_x, &b_RPCClusters_NN_x);
   fChain->SetBranchAddress("RPCClusters_NN.y", RPCClusters_NN_y, &b_RPCClusters_NN_y);
   fChain->SetBranchAddress("RPCClusters_NN.z", RPCClusters_NN_z, &b_RPCClusters_NN_z);
//   fChain->SetBranchAddress("RPCClusters_NN.idlist", RPCClusters_NN_idlist, &b_RPCClusters_NN_idlist);
   fChain->SetBranchAddress("RPCClusters_NN.valid", RPCClusters_NN_valid, &b_RPCClusters_NN_valid);
   fChain->SetBranchAddress("TCMTHits", &TCMTHits_, &b_DHCALEvent_TCMTHits_);
   fChain->SetBranchAddress("TCMTHits.fUniqueID", &TCMTHits_fUniqueID, &b_TCMTHits_fUniqueID);
   fChain->SetBranchAddress("TCMTHits.fBits", &TCMTHits_fBits, &b_TCMTHits_fBits);
   fChain->SetBranchAddress("TCMTHits.x", &TCMTHits_x, &b_TCMTHits_x);
   fChain->SetBranchAddress("TCMTHits.y", &TCMTHits_y, &b_TCMTHits_y);
   fChain->SetBranchAddress("TCMTHits.z", &TCMTHits_z, &b_TCMTHits_z);
   fChain->SetBranchAddress("TCMTHits.E", &TCMTHits_E, &b_TCMTHits_E);
   Notify();
}

Bool_t TreeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeAnalysis_cxx
