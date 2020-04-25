//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 25 12:41:47 2020 by ROOT version 6.14/09
// from TTree Events/Events
// found on file: ../OutputFiles/Jpsimc2017_numEvent5000.root
//////////////////////////////////////////////////////////

#ifndef convselMC_h
#define convselMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class convselMC : public TSelector {
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
   TTreeReaderArray<Bool_t> Conv_isConverted = {fReader, "Conv_isConverted"};
   TTreeReaderValue<UInt_t> nGenPart = {fReader, "nGenPart"};
   TTreeReaderArray<Float_t> GenPart_eta = {fReader, "GenPart_eta"};
   TTreeReaderArray<Float_t> GenPart_mass = {fReader, "GenPart_mass"};
   TTreeReaderArray<Float_t> GenPart_phi = {fReader, "GenPart_phi"};
   TTreeReaderArray<Float_t> GenPart_pt = {fReader, "GenPart_pt"};
   TTreeReaderArray<Float_t> GenPart_vtx_X = {fReader, "GenPart_vtx_X"};
   TTreeReaderArray<Float_t> GenPart_vtx_Y = {fReader, "GenPart_vtx_Y"};
   TTreeReaderArray<Float_t> GenPart_vtx_Z = {fReader, "GenPart_vtx_Z"};
   TTreeReaderArray<Int_t> GenPart_daughter0Idx = {fReader, "GenPart_daughter0Idx"};
   TTreeReaderArray<Int_t> GenPart_daughter1Idx = {fReader, "GenPart_daughter1Idx"};
   TTreeReaderArray<Int_t> GenPart_daughter2Idx = {fReader, "GenPart_daughter2Idx"};
   TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {fReader, "GenPart_genPartIdxMother"};
   TTreeReaderArray<Int_t> GenPart_nDaughter = {fReader, "GenPart_nDaughter"};
   TTreeReaderArray<Int_t> GenPart_pdgId = {fReader, "GenPart_pdgId"};
   TTreeReaderArray<Int_t> GenPart_status = {fReader, "GenPart_status"};
   TTreeReaderArray<Int_t> GenPart_statusFlags = {fReader, "GenPart_statusFlags"};
   TTreeReaderArray<Bool_t> GenPart_isConvertedPhoton = {fReader, "GenPart_isConvertedPhoton"};
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


   convselMC(TTree * /*tree*/ =0) { }
   virtual ~convselMC() { }
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

   ClassDef(convselMC,0);

};

#endif

#ifdef convselMC_cxx
void convselMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t convselMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


double dR(float eta1, float phi1, float eta2, float phi2) {
    
    float p1 = phi1;
    float p2 = phi2;
    float e1 = eta1;
    float e2 = eta2;
    auto dp = std::abs(p1 - p2);
    if (dp > float(M_PI))
      dp -= float(2 * M_PI);
    return sqrt ((e1 - e2) * (e1 - e2) + dp * dp);
  }
#endif // #ifdef convselMC_cxx
