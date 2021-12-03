// -*- C++ -*-
//
// Package:    MaterialStudy/PCNtupleMaker
// Class:      lumiprod
// 
/**\class lumiprod lumiprod.cc MaterialStudy/PCNtupleMaker/plugins/lumiprod.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Justin Anguiano
//         Created:  Fri, 09 Jul 2021 17:02:22 GMT
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

//
// class declaration
//

class lumiprod : public edm::stream::EDProducer<> {
   public:
      explicit lumiprod(const edm::ParameterSet&);
      ~lumiprod();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions){
		edm::ParameterSetDescription desc;
		desc.add<edm::InputTag>("lumiColl")->setComment("lumi collection");
		descriptions.add("lumip",desc);
	}

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
     
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
lumiprod::lumiprod(const edm::ParameterSet& iConfig)
{
   //register your products
   produces<double>( "evtTime" );//.setBranchAlias( "evtTime");
//   produces<LumiCollection>("LumiColl");
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


lumiprod::~lumiprod()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
lumiprod::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   double timeStamp_local = iEvent.time().unixTime();

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
    
    auto timeStamp_local_ptr = std::make_unique<double>(timeStamp_local);
    iEvent.put( std::move(timeStamp_local_ptr) ,"evtTime" ); 
 
  //std::auto_ptr<LumiCollection> LC( new LumiCollection );
  //LumiCollection* LCbase = new LumiCollection();
 // auto LCp = std::make_unique<LumiCollection>();
 // LCp->time = timeStamp_local;
  //iEvent.put( std::make_unique<LumiCollection>(LC), "LumiColl");
 // iEvent.put( std::move(LCp), "LumiColl");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
lumiprod::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
lumiprod::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
lumiprod::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
lumiprod::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
lumiprod::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
lumiprod::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
lumiprod::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}*/

//define this as a plug-in
DEFINE_FWK_MODULE(lumiprod);
