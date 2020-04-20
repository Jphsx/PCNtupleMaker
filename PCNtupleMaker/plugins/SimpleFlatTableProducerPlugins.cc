#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"


#include "DataFormats/EgammaCandidates/interface/Conversion.h"
typedef SimpleFlatTableProducer<reco::Conversion> SimpleConversionTableProducer;

#include "DataFormats/VertexReco/interface/Vertex.h"
typedef SimpleFlatTableProducer<reco::Vertex> SimpleVertexTableProducer;

//#include "DataFormats/Candidate/interface/Candidate.h"
//typedef SimpleFlatTableProducer<reco::Candidate> SimpleCandidateFlatTableProducer;

//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//typedef EventSingletonSimpleFlatTableProducer<GenEventInfoProduct> SimpleGenEventFlatTableProducer;

//#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
//typedef EventSingletonSimpleFlatTableProducer<HTXS::HiggsClassification> SimpleHTXSFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(SimpleCandidateFlatTableProducer);
//DEFINE_FWK_MODULE(SimpleGenEventFlatTableProducer);
//EFINE_FWK_MODULE(SimpleHTXSFlatTableProducer);
DEFINE_FWK_MODULE(SimpleConversionTableProducer);
DEFINE_FWK_MODULE(SimpleVertexTableProducer);
