// -*- C++ -*-
//
// Package:    MaterialStudy/PCNtupleMaker
// Class:      convTrkProducer
// 
/**\class convTrkProducer convTrkProducer.cc MaterialStudy/PCNtupleMaker/plugins/convTrkProducer.cc

 Description: [one line class summary]

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
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//

class convTrkProducer : public edm::stream::EDProducer<> {
   public:
      explicit convTrkProducer(const edm::ParameterSet&);
      ~convTrkProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<edm::View<reco::Conversion>> convs_;
      int trkCandIdx_;    

       
	
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
convTrkProducer::convTrkProducer(const edm::ParameterSet& iConfig): convs_(consumes<edm::View<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("src"))),
trkCandIdx_(iConfig.getParameter<int>("trkCandIdx")) 
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
   produces<std::vector<reco::Track>>("convTrks");
   //now do what ever other initialization is needed
  
}


convTrkProducer::~convTrkProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
convTrkProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  
  ///edm::RefToBase<reco::Track> cTrack;
  reco::Track cTrack;
  //std::vector<reco::Track> tracks;
  //std::unique_ptr<std::vector<edm::RefToBase<reco::Track>> > tracks( new std::vector<edm::RefToBase<reco::Track>> );
  std::unique_ptr<std::vector<reco::Track>> tracks( new std::vector<reco::Track> );

  View<reco::Conversion>::const_iterator conv, endconvs=convs->end();
  //View<reco::Track>::const_iterator trk;

  for (conv = convs->begin(); conv != endconvs; ++conv){
	
	conv->tracks();
//	for( trk = conv->tracks()->begin(); trk != conv->Tracks()->end(); ++trk){
	//cTrack = convs->tracks()[trkCandIdx_];	
	tracks->push_back(cTrack);
		

//	}

  }
  iEvent.put(std::move(tracks),"convTrks"); 




}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
convTrkProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
convTrkProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
convTrkProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
convTrkProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
convTrkProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
convTrkProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
convTrkProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(convTrkProducer);
