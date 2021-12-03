// -*- C++ -*-
//
// Package:    MaterialStudy/PCNtupleMaker
// Class:      convTrkSimTrkMatcher
// 
/**\class convTrkSimTrkMatcher convTrkSimTrkMatcher.cc MaterialStudy/PCNtupleMaker/plugins/convTrkSimTrkMatcher.cc

 Description: [one line class summary]
convTrkSimTrkMatcher
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Justin Anguiano
//         Created:  Thu, 23 Apr 2020 20:14:08 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
//
// class declaration
//

//class convTrkSimTrkMatcher : public edm::stream::EDProducer<> {
class convTrkSimTrkMatcher : public edm::EDProducer {
   public:
      explicit convTrkSimTrkMatcher(const edm::ParameterSet&);
      ~convTrkSimTrkMatcher();

     // static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:

     ///test remove 
     /*
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
     */

      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<edm::View<reco::Conversion>> convs_;
      const edm::EDGetTokenT<edm::View<SimTrack>> simtrk_;
      //int trkCandIdx_;    
      //float maxdr_;

       
	
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
convTrkSimTrkMatcher::convTrkSimTrkMatcher(const edm::ParameterSet& iConfig): convs_(consumes<edm::View<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("src"))), simtrk_(consumes<edm::View<SimTrack>>(iConfig.getParameter<edm::InputTag>("simsrc"))) 
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
//  produces<std::vector<edm::RefToBase<reco::Track>>>("convTrks");
 //  produces<std::vector<reco::Track>>("convTrks");
 
  // produces<std::vector<int>>("convVtxIdx");
  //  produces<std::vector<float>>("vtxdl");
	produces<edm::ValueMap<int> >("Tk0Idx");
	produces<edm::ValueMap<int> >("Tk1Idx");
	produces<edm::ValueMap<float> >("Tk0dR");
	produces<edm::ValueMap<float> >("Tk1dR");
	produces<edm::ValueMap<float> >("Tk0dPtRel");
	produces<edm::ValueMap<float> >("Tk1dPtRel");
	
 //now do what ever other initialization is needed
  
}


convTrkSimTrkMatcher::~convTrkSimTrkMatcher()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
convTrkSimTrkMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
  Handle<View<reco::Conversion> > convs;
  iEvent.getByToken(convs_, convs);

  Handle<View<SimTrack> > simtrk;
  iEvent.getByToken(simtrk_, simtrk);
 // std::cout<<"entered produce"<<std::endl;

  //LogDebug("Trace") >>"entered produce";
  //LogDebug("Values") >>"MyOtherStuff vector has size">>otherStuffs->size();
  //LogDebug("Trace") >>"exiting produce";
  
  ///edm::RefToBase<reco::Track> cTrack;
  edm::RefToBase<reco::Track> cTrk0;
  edm::RefToBase<reco::Track> cTrk1;
  // reco::Track cTrk0;
  // reco::Track cTrk1;
  

  //std::vector<reco::Track> tracks;
  //std::unique_ptr<std::vector<edm::RefToBase<reco::Track>> > tracks( new std::vector<edm::RefToBase<reco::Track>> );
 // std::unique_ptr<std::vector<reco::Track>> tracks( new std::vector<reco::Track> );

//  std::unique_ptr<std::vector<int>> simmatchidx( new std::vector<int> );
//  std::unique_ptr<std::vector<float>> simmatchdr( new std::vector<float>);
  std::vector<int> tk0matchidx;
  std::vector<int> tk1matchidx;
  std::vector<float> tk0matchdr;
  std::vector<float> tk1matchdr;
  std::vector<float> tk0dptrel;
  std::vector<float> tk1dptrel;

  View<reco::Conversion>::const_iterator conv, endconvs=convs->end();
  View<SimTrack>::const_iterator sim, endsim=simtrk->end();
  //View<reco::Track>::const_iterator trk;
  double pt0, eta0, phi0, pt1, eta1, phi1;
  double spt, seta, sphi;
  int q0,q1,sq;
  double dr,mindr,minidx,dptrel,mindptrel ;
  double maxdptrel = 0.4;
   mindptrel = 999;
  std::vector<double> in_range_vtx;
int i=0; 
  for (conv = convs->begin(); conv != endconvs; ++conv){
	

	//loop and check tk0
	cTrk0 = conv->tracks()[0];
	pt0 = cTrk0->pt();
	eta0 = cTrk0->eta();
	phi0 = cTrk0->phi();
	q0 = cTrk0->charge();
	mindr = 999;
	dr = 9999;
	minidx = -1;
	i=0;	
	for( sim = simtrk->begin(); sim != endsim; ++sim){
		spt = sim->momentum().Pt();
		seta = sim->momentum().Eta();
		sphi = sim->momentum().Phi();
		sq = int(sim->charge());
		if(sq != q0) continue;
		dr = deltaR(eta0, phi0, seta, sphi);
		dptrel = abs((spt-pt0))/(spt);
		if(dr < mindr && dptrel < maxdptrel){
			mindr = dr;
			minidx = i;
			mindptrel = dptrel;
		}	
		i++;
	}	
	tk0matchidx.push_back(minidx);
	tk0matchdr.push_back(mindr);		
	tk0dptrel.push_back(mindptrel);


	//loop and check tk1
	cTrk1 = conv->tracks()[1];
	pt1 = cTrk1->pt();
	eta1 = cTrk1->eta();
	phi1 = cTrk1->phi();
	q1 = cTrk1->charge();
	mindr = 999;
	dr = 9999;
	minidx = -1;
	i=0;	
	for( sim = simtrk->begin(); sim != endsim; ++sim){
		
		spt = sim->momentum().Pt();
		seta = sim->momentum().Eta();
		sphi = sim->momentum().Phi();
		sq = int(sim->charge());
		if(sq != q1) continue;
		dr = deltaR(eta1, phi1, seta, sphi);
		dptrel = abs((spt-pt1))/(spt);
		if(dr < mindr && dptrel < maxdptrel){
			mindr = dr;
			minidx = i;
			mindptrel = dptrel;
		}	
		i++;
	}	
	tk1matchidx.push_back(minidx);
	tk1matchdr.push_back(mindr);		
	tk1dptrel.push_back(mindptrel);


  }
  //iEvent.put(std::move(tracks),"convTrks"); 
  //
  std::unique_ptr<edm::ValueMap<int> > map0idx( new edm::ValueMap<int>());
  std::unique_ptr<edm::ValueMap<float> >map0dr( new edm::ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> >map0dpt( new edm::ValueMap<float>());
  edm::ValueMap<int>::Filler tk0filler1(*map0idx);
  edm::ValueMap<float>::Filler tk0filler2(*map0dr); 
  edm::ValueMap<float>::Filler tk0filler3(*map0dpt);
  tk0filler1.insert( convs, tk0matchidx.begin(), tk0matchidx.end());
  tk0filler1.fill();
  tk0filler2.insert( convs, tk0matchdr.begin(), tk0matchdr.end());
  tk0filler2.fill();
  tk0filler3.insert( convs, tk0dptrel.begin(), tk0dptrel.end());
  tk0filler3.fill();

  std::unique_ptr<edm::ValueMap<int> > map1idx( new edm::ValueMap<int>());
  std::unique_ptr<edm::ValueMap<float> >map1dr( new edm::ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> >map1dpt( new edm::ValueMap<float>());
  edm::ValueMap<int>::Filler tk1filler1(*map1idx);
  edm::ValueMap<float>::Filler tk1filler2(*map1dr); 
  edm::ValueMap<float>::Filler tk1filler3(*map1dpt);
  tk1filler1.insert( convs, tk1matchidx.begin(), tk1matchidx.end());
  tk1filler1.fill();
  tk1filler2.insert( convs, tk1matchdr.begin(), tk1matchdr.end());
  tk1filler2.fill();
  tk1filler3.insert( convs, tk1dptrel.begin(), tk1dptrel.end());
  tk1filler3.fill();


  iEvent.put(std::move(map0idx), "Tk0Idx");
  iEvent.put(std::move(map0dr), "Tk0dR");
  iEvent.put(std::move(map0dpt), "Tk0dPtRel");
 
  iEvent.put(std::move(map1idx), "Tk1Idx");
  iEvent.put(std::move(map1dr), "Tk1dR");
  iEvent.put(std::move(map1dpt), "Tk1dPtRel");
 


//  iEvent.put(std::move(simmatchidx),"convVtxIdx");
//  iEvent.put(std::move(simmatchdr),"vtxdl");



}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
//convTrkSimTrkMatcher::beginStream(edm::StreamID)
convTrkSimTrkMatcher::beginJob()
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
//convTrkSimTrkMatcher::endStream() {
convTrkSimTrkMatcher::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
convTrkSimTrkMatcher::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
convTrkSimTrkMatcher::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
convTrkSimTrkMatcher::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
convTrkSimTrkMatcher::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
convTrkSimTrkMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(convTrkSimTrkMatcher);
