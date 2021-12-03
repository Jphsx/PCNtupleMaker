//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 26 14:53:52 2020 by ROOT version 6.14/09
// from TTree Events/Events
// found on file: defaultout_numEvent100_64.root
//////////////////////////////////////////////////////////

#ifndef PCNtuple_h
#define PCNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TH2D.h>
// Headers needed by this particular selector


class PCNtuple : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UInt_t> run = {fReader, "run"};
   TTreeReaderValue<UInt_t> luminosityBlock = {fReader, "luminosityBlock"};
   TTreeReaderValue<ULong64_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> nConv = {fReader, "nConv"};
   TTreeReaderArray<Float_t> Conv_EoverP = {fReader, "Conv_EoverP"};
   TTreeReaderArray<Float_t> Conv_EoverPrefittedTracks = {fReader, "Conv_EoverPrefittedTracks"};
   TTreeReaderArray<Float_t> Conv_Tk0_chi2 = {fReader, "Conv_Tk0_chi2"};
   TTreeReaderArray<Float_t> Conv_Tk0_eta = {fReader, "Conv_Tk0_eta"};
   TTreeReaderArray<Float_t> Conv_Tk0_ndof = {fReader, "Conv_Tk0_ndof"};
   TTreeReaderArray<Float_t> Conv_Tk0_normalizedChi2 = {fReader, "Conv_Tk0_normalizedChi2"};
   TTreeReaderArray<Float_t> Conv_Tk0_phi = {fReader, "Conv_Tk0_phi"};
   TTreeReaderArray<Float_t> Conv_Tk0_pt = {fReader, "Conv_Tk0_pt"};
   TTreeReaderArray<Float_t> Conv_Tk1_chi2 = {fReader, "Conv_Tk1_chi2"};
   TTreeReaderArray<Float_t> Conv_Tk1_eta = {fReader, "Conv_Tk1_eta"};
   TTreeReaderArray<Float_t> Conv_Tk1_ndof = {fReader, "Conv_Tk1_ndof"};
   TTreeReaderArray<Float_t> Conv_Tk1_normalizedChi2 = {fReader, "Conv_Tk1_normalizedChi2"};
   TTreeReaderArray<Float_t> Conv_Tk1_phi = {fReader, "Conv_Tk1_phi"};
   TTreeReaderArray<Float_t> Conv_Tk1_pt = {fReader, "Conv_Tk1_pt"};
   TTreeReaderArray<Float_t> Conv_dPhiTracksAtVtx = {fReader, "Conv_dPhiTracksAtVtx"};
   TTreeReaderArray<Float_t> Conv_distOfMinimumApproach = {fReader, "Conv_distOfMinimumApproach"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_Tk0 = {fReader, "Conv_dlClosestHitToVtx_Tk0"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_Tk1 = {fReader, "Conv_dlClosestHitToVtx_Tk1"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_err_Tk0 = {fReader, "Conv_dlClosestHitToVtx_err_Tk0"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_err_Tk1 = {fReader, "Conv_dlClosestHitToVtx_err_Tk1"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_sig_Tk0 = {fReader, "Conv_dlClosestHitToVtx_sig_Tk0"};
   TTreeReaderArray<Float_t> Conv_dlClosestHitToVtx_sig_Tk1 = {fReader, "Conv_dlClosestHitToVtx_sig_Tk1"};
   TTreeReaderArray<Float_t> Conv_dxy = {fReader, "Conv_dxy"};
   TTreeReaderArray<Float_t> Conv_dz = {fReader, "Conv_dz"};
   TTreeReaderArray<Float_t> Conv_lxy = {fReader, "Conv_lxy"};
   TTreeReaderArray<Float_t> Conv_lz = {fReader, "Conv_lz"};
   TTreeReaderArray<Float_t> Conv_pairCotThetaSeparation = {fReader, "Conv_pairCotThetaSeparation"};
   TTreeReaderArray<Float_t> Conv_pairInvariantMass = {fReader, "Conv_pairInvariantMass"};
   TTreeReaderArray<Float_t> Conv_pairMomentum_Px = {fReader, "Conv_pairMomentum_Px"};
   TTreeReaderArray<Float_t> Conv_pairMomentum_Py = {fReader, "Conv_pairMomentum_Py"};
   TTreeReaderArray<Float_t> Conv_pairMomentum_Pz = {fReader, "Conv_pairMomentum_Pz"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_E = {fReader, "Conv_refittedPair4Momentum_E"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_M = {fReader, "Conv_refittedPair4Momentum_M"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_Px = {fReader, "Conv_refittedPair4Momentum_Px"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_Py = {fReader, "Conv_refittedPair4Momentum_Py"};
   TTreeReaderArray<Float_t> Conv_refittedPair4Momentum_Pz = {fReader, "Conv_refittedPair4Momentum_Pz"};
   TTreeReaderArray<Float_t> Conv_refittedPairMomentum_Px = {fReader, "Conv_refittedPairMomentum_Px"};
   TTreeReaderArray<Float_t> Conv_refittedPairMomentum_Py = {fReader, "Conv_refittedPairMomentum_Py"};
   TTreeReaderArray<Float_t> Conv_refittedPairMomentum_Pz = {fReader, "Conv_refittedPairMomentum_Pz"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_X_Tk0 = {fReader, "Conv_tracksInnerPosition_X_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_X_Tk1 = {fReader, "Conv_tracksInnerPosition_X_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Y_Tk0 = {fReader, "Conv_tracksInnerPosition_Y_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Y_Tk1 = {fReader, "Conv_tracksInnerPosition_Y_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Z_Tk0 = {fReader, "Conv_tracksInnerPosition_Z_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksInnerPosition_Z_Tk1 = {fReader, "Conv_tracksInnerPosition_Z_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Px_Tk0 = {fReader, "Conv_tracksPin_Px_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Px_Tk1 = {fReader, "Conv_tracksPin_Px_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Py_Tk0 = {fReader, "Conv_tracksPin_Py_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Py_Tk1 = {fReader, "Conv_tracksPin_Py_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Pz_Tk0 = {fReader, "Conv_tracksPin_Pz_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPin_Pz_Tk1 = {fReader, "Conv_tracksPin_Pz_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Px_Tk0 = {fReader, "Conv_tracksPout_Px_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Px_Tk1 = {fReader, "Conv_tracksPout_Px_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Py_Tk0 = {fReader, "Conv_tracksPout_Py_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Py_Tk1 = {fReader, "Conv_tracksPout_Py_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Pz_Tk0 = {fReader, "Conv_tracksPout_Pz_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksPout_Pz_Tk1 = {fReader, "Conv_tracksPout_Pz_Tk1"};
   TTreeReaderArray<Float_t> Conv_tracksSigned_d0_Tk0 = {fReader, "Conv_tracksSigned_d0_Tk0"};
   TTreeReaderArray<Float_t> Conv_tracksSigned_d0_Tk1 = {fReader, "Conv_tracksSigned_d0_Tk1"};
   TTreeReaderArray<Float_t> Conv_vtx_X = {fReader, "Conv_vtx_X"};
   TTreeReaderArray<Float_t> Conv_vtx_Y = {fReader, "Conv_vtx_Y"};
   TTreeReaderArray<Float_t> Conv_vtx_Z = {fReader, "Conv_vtx_Z"};
   TTreeReaderArray<Float_t> Conv_vtx_chi2 = {fReader, "Conv_vtx_chi2"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_00 = {fReader, "Conv_vtx_cov_00"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_01 = {fReader, "Conv_vtx_cov_01"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_02 = {fReader, "Conv_vtx_cov_02"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_11 = {fReader, "Conv_vtx_cov_11"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_12 = {fReader, "Conv_vtx_cov_12"};
   TTreeReaderArray<Float_t> Conv_vtx_cov_22 = {fReader, "Conv_vtx_cov_22"};
   TTreeReaderArray<Float_t> Conv_vtx_ndof = {fReader, "Conv_vtx_ndof"};
   TTreeReaderArray<Float_t> Conv_vtx_normalizedChi2 = {fReader, "Conv_vtx_normalizedChi2"};
   TTreeReaderArray<Float_t> Conv_zOfPrimaryVertexFromTracks = {fReader, "Conv_zOfPrimaryVertexFromTracks"};
   TTreeReaderArray<Float_t> Conv_Tk0_dPtRel = {fReader, "Conv_Tk0_dPtRel"};
   TTreeReaderArray<Float_t> Conv_Tk0_dR = {fReader, "Conv_Tk0_dR"};
   TTreeReaderArray<Float_t> Conv_Tk1_dPtRel = {fReader, "Conv_Tk1_dPtRel"};
   TTreeReaderArray<Float_t> Conv_Tk1_dR = {fReader, "Conv_Tk1_dR"};
   TTreeReaderArray<Float_t> Conv_vtxdl = {fReader, "Conv_vtxdl"};
   TTreeReaderArray<Int_t> Conv_Tk0_algo = {fReader, "Conv_Tk0_algo"};
   TTreeReaderArray<Int_t> Conv_Tk0_charge = {fReader, "Conv_Tk0_charge"};
   TTreeReaderArray<Int_t> Conv_Tk0_found = {fReader, "Conv_Tk0_found"};
   TTreeReaderArray<Int_t> Conv_Tk0_lost = {fReader, "Conv_Tk0_lost"};
   TTreeReaderArray<Int_t> Conv_Tk0_quality = {fReader, "Conv_Tk0_quality"};
   TTreeReaderArray<Int_t> Conv_Tk1_algo = {fReader, "Conv_Tk1_algo"};
   TTreeReaderArray<Int_t> Conv_Tk1_charge = {fReader, "Conv_Tk1_charge"};
   TTreeReaderArray<Int_t> Conv_Tk1_found = {fReader, "Conv_Tk1_found"};
   TTreeReaderArray<Int_t> Conv_Tk1_lost = {fReader, "Conv_Tk1_lost"};
   TTreeReaderArray<Int_t> Conv_Tk1_quality = {fReader, "Conv_Tk1_quality"};
   TTreeReaderArray<Int_t> Conv_algo = {fReader, "Conv_algo"};
   TTreeReaderArray<Int_t> Conv_nHitsBeforeVtx_Tk0 = {fReader, "Conv_nHitsBeforeVtx_Tk0"};
   TTreeReaderArray<Int_t> Conv_nHitsBeforeVtx_Tk1 = {fReader, "Conv_nHitsBeforeVtx_Tk1"};
   TTreeReaderArray<Int_t> Conv_nSharedHits = {fReader, "Conv_nSharedHits"};
   TTreeReaderArray<Int_t> Conv_nTracks = {fReader, "Conv_nTracks"};
   TTreeReaderArray<Int_t> Conv_Tk0_Idx = {fReader, "Conv_Tk0_Idx"};
   TTreeReaderArray<Int_t> Conv_Tk1_Idx = {fReader, "Conv_Tk1_Idx"};
   TTreeReaderArray<Int_t> Conv_convVtxIdx = {fReader, "Conv_convVtxIdx"};
   TTreeReaderArray<Bool_t> Conv_isConverted = {fReader, "Conv_isConverted"};
   TTreeReaderValue<UInt_t> nPV = {fReader, "nPV"};
   TTreeReaderArray<Float_t> PV_X = {fReader, "PV_X"};
   TTreeReaderArray<Float_t> PV_Y = {fReader, "PV_Y"};
   TTreeReaderArray<Float_t> PV_Z = {fReader, "PV_Z"};
   TTreeReaderArray<Float_t> PV_chi2 = {fReader, "PV_chi2"};
   TTreeReaderArray<Float_t> PV_cov_00 = {fReader, "PV_cov_00"};
   TTreeReaderArray<Float_t> PV_cov_01 = {fReader, "PV_cov_01"};
   TTreeReaderArray<Float_t> PV_cov_02 = {fReader, "PV_cov_02"};
   TTreeReaderArray<Float_t> PV_cov_11 = {fReader, "PV_cov_11"};
   TTreeReaderArray<Float_t> PV_cov_12 = {fReader, "PV_cov_12"};
   TTreeReaderArray<Float_t> PV_cov_22 = {fReader, "PV_cov_22"};
   TTreeReaderArray<Float_t> PV_ndof = {fReader, "PV_ndof"};
   TTreeReaderArray<Float_t> PV_normalizedChi2 = {fReader, "PV_normalizedChi2"};
   TTreeReaderValue<UInt_t> nSimTrk = {fReader, "nSimTrk"};
   TTreeReaderArray<Float_t> SimTrk_charge = {fReader, "SimTrk_charge"};
   TTreeReaderArray<Float_t> SimTrk_eta = {fReader, "SimTrk_eta"};
   TTreeReaderArray<Float_t> SimTrk_phi = {fReader, "SimTrk_phi"};
   TTreeReaderArray<Float_t> SimTrk_pt = {fReader, "SimTrk_pt"};
   TTreeReaderArray<Int_t> SimTrk_pdgId = {fReader, "SimTrk_pdgId"};
   TTreeReaderArray<Int_t> SimTrk_simvtx_Idx = {fReader, "SimTrk_simvtx_Idx"};
   TTreeReaderArray<Int_t> SimTrk_trackId = {fReader, "SimTrk_trackId"};
   TTreeReaderValue<UInt_t> nSimVtx = {fReader, "nSimVtx"};
   TTreeReaderArray<Float_t> SimVtx_tof = {fReader, "SimVtx_tof"};
   TTreeReaderArray<Float_t> SimVtx_x = {fReader, "SimVtx_x"};
   TTreeReaderArray<Float_t> SimVtx_y = {fReader, "SimVtx_y"};
   TTreeReaderArray<Float_t> SimVtx_z = {fReader, "SimVtx_z"};
   TTreeReaderArray<Int_t> SimVtx_processType = {fReader, "SimVtx_processType"};
   TTreeReaderArray<Int_t> SimVtx_simtrk_parent_tid = {fReader, "SimVtx_simtrk_parent_tid"};


   PCNtuple(TTree * /*tree*/ =0) { }
   virtual ~PCNtuple() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(PCNtuple,0);

   TH1D* nPC=0;
   TH1D* vtxdl=0;
   TH2D* xy10=0;
   TH2D* xy25=0;
   TH1D* r10=0;
   TH1D* r25=0;
  
   //photon eff stuff
  //these are pure sim
   TH2D* xy10s=0;
   TH2D* xy25s=0;
   TH1D* r10s=0;
   TH1D* r25s=0;
  



   //these are reco matched to sims
   TH2D* xy10sim=0;
   TH2D* xy10sim_cut=0;
   TH2D* xy25sim=0;
   TH2D* xy25sim_cut=0;
   TH1D* r10sim=0;
   TH1D* r10sim_cut=0;
   TH1D* r25sim=0;
   TH1D* r25sim_cut=0;
   TH1D* numsimphot_raw=0; // number over all events
  

   TH1D* numsimconv_raw=0;
   TH1D* numsimconv_cut=0;//numer over all events

   TH1D* nump14=0;
   TH1D* numchild=0;
 
   TH1D* numprimarypart=0;
   TH1D* numprimaryq=0;
   TH1D* numprimaryphot=0;
   TH1D* numsecondaryvtx=0;
   TH1D* secondarypt=0;
   TH1D* secondaryeta=0;

   TH1D* nsimphot_raw=0;
   TH1D* nsimconv_raw=0;//number per event
   TH1D* nsimconv_cut=0;
   TH1D* simphotpt_raw=0;
   //TH1D* simphotpt_cut=0;
   TH1D* simconvpt_raw=0;
   TH1D* simconvpt_cut=0;
   TH1D* simconvR_raw=0;
   TH1D* simconvR_cut=0;
   TH2D* simconvptR_cut=0;


 ///  TH1D* numpair=0;
 //  TH1D* numtrip=0;
//   int n_sim_pc=0;
//   int n_sim_pc_rec=0;
//   int n_sim_pc_rec_cut=0;

//eff counters for 
//denom
   TH1D* p14=0;
   TH1D* npc=0;
//nearest vtx is p14 from mask and passes dl cut
   TH1D* recomatchp14=0;
//matches any vtx with dl
   TH1D* recomatchany=0;
//matches bg
  TH1D* recomatchbg=0;
//matches p14 and dl and both tracks
   TH1D* recomatchp14alltrk=0;
//matches p14 adn dl and at least one track
   TH1D* recomatchp14onetrk=0;
//duplicate match counter 2 recos match to the same p14 within dl
   TH1D* recomatchp14dup=0;
  
//nearest vtx is p14 from mask and passes dl cut
   TH1D* p14_pt=0;
   TH1D* npc_pt=0;
   TH1D* recomatchp14_pt=0;
//matches any vtx with dl
   TH1D* recomatchany_pt=0;
//matches bg
  TH1D* recomatchbg_pt=0;
//matches p14 and dl and both tracks
   TH1D* recomatchp14alltrk_pt=0;
//matches p14 adn dl and at least one track
   TH1D* recomatchp14onetrk_pt=0;
//duplicate match counter 2 recos match to the same p14 within dl
 //  TH1D* recomatchp14dup=0;

//nearest vtx is p14 from mask and passes dl cut
   TH1D* p14_r=0;
   TH1D* npc_r=0;
   TH1D* recomatchp14_r=0;
//matches any vtx with dl
   TH1D* recomatchany_r=0;
//matches bg
  TH1D* recomatchbg_r=0;
//matches p14 and dl and both tracks
   TH1D* recomatchp14alltrk_r=0;
//matches p14 adn dl and at least one track
   TH1D* recomatchp14onetrk_r=0;
//duplicate match counter 2 recos match to the same p14 within dl
 ///  TH1D* recomatchp14dup_r=0;

 
//nearest vtx is p14 from mask and passes dl cut
   TH2D* p14_ptr=0;
   TH2D* npc_ptr=0;
   TH2D* recomatchp14_ptr=0;
//matches any vtx with dl
   TH2D* recomatchany_ptr=0;
//matches bg
  TH2D* recomatchbg_ptr=0;
//matches p14 and dl and both tracks
   TH2D* recomatchp14alltrk_ptr=0;
//matches p14 adn dl and at least one track
   TH2D* recomatchp14onetrk_ptr=0;
//duplicate match counter 2 recos match to the same p14 within dl
 ///  TH1D* recomatchp14dup_r=0;

 
//newest eff stuff...
   TH1D* totnpc_raw=0;
   TH1D* totnpc_cut=0;
   TH2D* xypc_cut=0;	
   TH1D* rpc_cut=0;
   TH2D* xypc_raw=0;
   TH1D* rpc_raw=0;
   //1D num and denoms
   TH1D* ptdenom=0;
   TH1D* ptnum=0;
   TH1D* rdenom=0;
   TH1D* rnum=0; 

   //2d num and denoms 
   TH2D* ptrnum=0;
   TH2D* ptrdenom=0;

   //count mismatching
   TH1D* totsig_cut=0;
   TH1D* totbg_cut=0;
   TH1D* totunmatch_cut=0;
   TH1D* totsig_raw=0;
   TH1D* totbg_raw=0;
   TH1D* totunmatch_raw=0;

   TH1D* thetacheck=0;
int ctr=0;

};

#endif

#ifdef PCNtuple_cxx
void PCNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t PCNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef PCNtuple_cxx
