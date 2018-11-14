#include "RecoEgamma/EgammaElectronProducers/plugins/LowPtGsfElectronCoreProducer.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <map>

using namespace reco ;

LowPtGsfElectronCoreProducer::LowPtGsfElectronCoreProducer( const edm::ParameterSet& config )
  : GsfElectronCoreBaseProducer(config)
{
  superClusters_ = consumes<reco::SuperClusterCollection>(config.getParameter<edm::InputTag>("superClusters"));
  superClusterRefs_ = consumes< edm::ValueMap<reco::SuperClusterRef> >(config.getParameter<edm::InputTag>("superClusters"));
}

LowPtGsfElectronCoreProducer::~LowPtGsfElectronCoreProducer() {}

void LowPtGsfElectronCoreProducer::produce( edm::Event& event, const edm::EventSetup& setup ) {

  // Output collection
  auto electrons = std::make_unique<GsfElectronCoreCollection>();

  // Init
  GsfElectronCoreBaseProducer::initEvent(event,setup) ;
  if ( !useGsfPfRecTracks_ ) { edm::LogError("useGsfPfRecTracks_ is (redundantly) set to False!"); }
  if ( !gsfPfRecTracksH_.isValid() ) { edm::LogError("gsfPfRecTracks handle is invalid!"); }
  if ( !gsfTracksH_.isValid() ) { edm::LogError("gsfTracks handle is invalid!"); }

  edm::Handle<reco::SuperClusterCollection> superClusters;
  event.getByToken(superClusters_,superClusters);
  if ( !superClusters.isValid() ) { edm::LogError("Problem with superClusters handle"); }

  edm::Handle< edm::ValueMap<reco::SuperClusterRef> > superClusterRefs;
  event.getByToken(superClusterRefs_,superClusterRefs);
  if ( !superClusterRefs.isValid() ) { edm::LogError("Problem with superClusterRefs handle"); }

  // Create ElectronCore objects
  for ( size_t ipfgsf = 0; ipfgsf < gsfPfRecTracksH_->size(); ++ipfgsf ) {

    // Refs to GSF(PF) objects and SC
    reco::GsfPFRecTrackRef pfgsf(gsfPfRecTracksH_, ipfgsf);
    reco::GsfTrackRef gsf = pfgsf->gsfTrackRef();
    const reco::SuperClusterRef sc = (*superClusterRefs)[pfgsf];

    // Construct and keep ElectronCore if GSF(PF) track and SC are present
    GsfElectronCore* core = new GsfElectronCore(gsf);
    if ( core->ecalDrivenSeed() ) { delete core; return; }

    // Add GSF(PF) track information
    GsfElectronCoreBaseProducer::fillElectronCore(core);

    // Add super cluster
    core->setSuperCluster(sc);

    // Store
    electrons->push_back(*core);

  }

  event.put(std::move(electrons));

}

void LowPtGsfElectronCoreProducer::fillDescription( edm::ParameterSetDescription& desc )
{
  desc.add<edm::InputTag>("gsfPfRecTracks",edm::InputTag("lowPtGsfElePfGsfTracks")) ;
  desc.add<edm::InputTag>("gsfTracks",edm::InputTag("lowPtGsfEleGsfTracks")) ;
  desc.add<edm::InputTag>("ctfTracks",edm::InputTag("generalTracks")) ;
  desc.add<bool>("useGsfPfRecTracks",true) ;
}
