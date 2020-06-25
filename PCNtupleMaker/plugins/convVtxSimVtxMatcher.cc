// -*- C++ -*-
//
// Package:    MaterialStudy/PCNtupleMaker
// Class:      convVtxSimVtxMatcher
// 
/**\class convVtxSimVtxMatcher convVtxSimVtxMatcher.cc MaterialStudy/PCNtupleMaker/plugins/convVtxSimVtxMatcher.cc

 Description: [one line class summary]
convVtxSimVtxMatcher
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
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "DataFormats/Common/interface/ValueMap.h"
//
// class declaration
//

//class convVtxSimVtxMatcher : public edm::stream::EDProducer<> {
class convVtxSimVtxMatcher : public edm::EDProducer {
   public:
      explicit convVtxSimVtxMatcher(const edm::ParameterSet&);
      ~convVtxSimVtxMatcher();

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
      const edm::EDGetTokenT<edm::View<SimVertex>> simvtx_;
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
convVtxSimVtxMatcher::convVtxSimVtxMatcher(const edm::ParameterSet& iConfig): convs_(consumes<edm::View<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("src"))), simvtx_(consumes<edm::View<SimVertex>>(iConfig.getParameter<edm::InputTag>("simsrc"))) 
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
	produces<edm::ValueMap<int> >("convVtxIdx");
	produces<edm::ValueMap<float> >("vtxdl");  
 //now do what ever other initialization is needed
  
}


convVtxSimVtxMatcher::~convVtxSimVtxMatcher()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
convVtxSimVtxMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  Handle<View<SimVertex> > simvtx;
  iEvent.getByToken(simvtx_, simvtx);
 // std::cout<<"entered produce"<<std::endl;

  //LogDebug("Trace") >>"entered produce";
  //LogDebug("Values") >>"MyOtherStuff vector has size">>otherStuffs->size();
  //LogDebug("Trace") >>"exiting produce";
  
  ///edm::RefToBase<reco::Track> cTrack;
 // reco::Track cTrack;
  
  //std::vector<reco::Track> tracks;
  //std::unique_ptr<std::vector<edm::RefToBase<reco::Track>> > tracks( new std::vector<edm::RefToBase<reco::Track>> );
 // std::unique_ptr<std::vector<reco::Track>> tracks( new std::vector<reco::Track> );

//  std::unique_ptr<std::vector<int>> simmatchidx( new std::vector<int> );
//  std::unique_ptr<std::vector<float>> simmatchdr( new std::vector<float>);
  std::vector<int> simmatchidx;
  std::vector<float> simmatchdr;

  View<reco::Conversion>::const_iterator conv, endconvs=convs->end();
  View<SimVertex>::const_iterator sim, endsim=simvtx->end();
  //View<reco::Track>::const_iterator trk;
  double c_x, c_y, c_z;
  double s_x, s_y, s_z;
  double dr,mindr,minidx;
  std::vector<double> in_range_vtx;
int i=0; 
  for (conv = convs->begin(); conv != endconvs; ++conv){
	
		
	//get conv vtx
	c_x = conv->conversionVertex().x();
	c_y = conv->conversionVertex().y();
	c_z = conv->conversionVertex().z();

	mindr = 999;
	minidx = -1;
	i=0;
	for( sim = simvtx->begin(); sim != endsim; ++sim){
		s_x = sim->position().x();
		s_y = sim->position().y();
		s_z = sim->position().z();
		dr = sqrt( (s_x - c_x)*(s_x - c_x) + (s_y - c_y)*(s_y - c_y) + (s_z - c_z)*(s_z - c_z) );
		if(dr < mindr){
			mindr = dr;
			minidx = i;	
		}	
		i++;
	}	

	simmatchidx.push_back(minidx);
	simmatchdr.push_back(mindr);		

//	conv->tracks();
//	for( trk = conv->tracks()->begin(); trk != conv->Tracks()->end(); ++trk){
	//cTrack = convs->tracks()[trkCandIdx_];	
//	tracks->push_back(cTrack);
//	if(i==0) std::cout<<cTrack.pt()<<" "<<conv->EoverP()<<std::endl;
//	i++;

//	}

  }
  //iEvent.put(std::move(tracks),"convTrks"); 
  //
  std::unique_ptr<edm::ValueMap<int> > mapidx( new edm::ValueMap<int>());
  std::unique_ptr<edm::ValueMap<float> >mapdr( new edm::ValueMap<float>());
  edm::ValueMap<int>::Filler filler1(*mapidx);
  edm::ValueMap<float>::Filler filler2(*mapdr); 
  filler1.insert( convs, simmatchidx.begin(), simmatchidx.end());
  filler1.fill();
  filler2.insert( convs, simmatchdr.begin(), simmatchdr.end());
  filler2.fill();
  iEvent.put(std::move(mapidx), "convVtxIdx");
  iEvent.put(std::move(mapdr), "vtxdl");
 

//  iEvent.put(std::move(simmatchidx),"convVtxIdx");
//  iEvent.put(std::move(simmatchdr),"vtxdl");



}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
//convVtxSimVtxMatcher::beginStream(edm::StreamID)
convVtxSimVtxMatcher::beginJob()
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
//convVtxSimVtxMatcher::endStream() {
convVtxSimVtxMatcher::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
convVtxSimVtxMatcher::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
convVtxSimVtxMatcher::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
convVtxSimVtxMatcher::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
convVtxSimVtxMatcher::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
convVtxSimVtxMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(convVtxSimVtxMatcher);
