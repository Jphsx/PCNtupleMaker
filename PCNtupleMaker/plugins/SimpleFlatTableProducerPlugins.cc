#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"


#include "DataFormats/EgammaCandidates/interface/Conversion.h"
typedef SimpleFlatTableProducer<reco::Conversion> SimpleConversionTableProducer;

#include "DataFormats/VertexReco/interface/Vertex.h"
typedef SimpleFlatTableProducer<reco::Vertex> SimpleVertexTableProducer;

#include "DataFormats/Candidate/interface/Candidate.h"
typedef SimpleFlatTableProducer<reco::Candidate> SimpleCandidateFlatTableProducer;

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
typedef SimpleFlatTableProducer<reco::GenParticle> SimpleGenParticleFlatTableProducer;

#include "SimDataFormats/Track/interface/SimTrack.h"
typedef SimpleFlatTableProducer<SimTrack> SimpleSimTrackFlatTableProducer;

#include "SimDataFormats/Vertex/interface/SimVertex.h"
typedef SimpleFlatTableProducer<SimVertex> SimpleSimVertexFlatTableProducer;

typedef SimpleFlatTableProducer<int> SimpleIntTableProducer;

typedef SimpleFlatTableProducer<float> SimpleFloatTableProducer;

//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//typedef EventSingletonSimpleFlatTableProducer<GenEventInfoProduct> SimpleGenEventFlatTableProducer;

//#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
//typedef EventSingletonSimpleFlatTableProducer<HTXS::HiggsClassification> SimpleHTXSFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimpleCandidateFlatTableProducer);
DEFINE_FWK_MODULE(SimpleIntTableProducer);
DEFINE_FWK_MODULE(SimpleFloatTableProducer);
//DEFINE_FWK_MODULE(SimpleGenEventFlatTableProducer);
//EFINE_FWK_MODULE(SimpleHTXSFlatTableProducer);
DEFINE_FWK_MODULE(SimpleConversionTableProducer);
DEFINE_FWK_MODULE(SimpleVertexTableProducer);
DEFINE_FWK_MODULE(SimpleSimTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleSimVertexFlatTableProducer);
//DEFINE_FWK_MODULE(SimpleGenParticleTableProducer);

