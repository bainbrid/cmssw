#ifndef LowPtGsfElectronSeedProducer_h
#define LowPtGsfElectronSeedProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
class BaseParticlePropagator;

class LowPtGsfElectronSeedProducer : public edm::stream::EDProducer<> 
{
  
 public:
  
  explicit LowPtGsfElectronSeedProducer( const edm::ParameterSet& );
  ~LowPtGsfElectronSeedProducer() override ;
  
  void beginRun( edm::Run const&, edm::EventSetup const& );
  void endRun( edm::Run const&, edm::EventSetup const& );
  
  void produce( edm::Event&, const edm::EventSetup& );
  
  static void fillDescriptions( edm::ConfigurationDescriptions& );
  
 private:
  
  reco::PreId create_preid( const reco::TrackRef&, 
			    const edm::Handle<reco::PFClusterCollection>& );
  
  math::XYZVector field_;
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  edm::EDGetTokenT<reco::PFClusterCollection> clusters_;

};

#endif // LowPtGsfElectronSeedProducer_h
