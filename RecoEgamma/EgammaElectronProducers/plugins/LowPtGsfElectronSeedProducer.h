#ifndef RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedProducer_h
#define RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"

class LowPtGsfElectronSeedProducer final : public edm::stream::EDProducer< edm::GlobalCache<lowptgsfeleseed::HeavyObjectCache> >
{
  
 public:
  
  explicit LowPtGsfElectronSeedProducer( const edm::ParameterSet&, 
					 const lowptgsfeleseed::HeavyObjectCache* );

  ~LowPtGsfElectronSeedProducer() override;
  
  static std::unique_ptr<lowptgsfeleseed::HeavyObjectCache> 
    initializeGlobalCache( const edm::ParameterSet& conf ) {
    return std::make_unique<lowptgsfeleseed::HeavyObjectCache>(lowptgsfeleseed::HeavyObjectCache(conf));
  }
  
  static void globalEndJob( lowptgsfeleseed::HeavyObjectCache const* ) {}
  
  void beginLuminosityBlock( edm::LuminosityBlock const&, 
			     edm::EventSetup const& ) override;

  void produce( edm::Event&, const edm::EventSetup& ) override;
  
  static void fillDescription( edm::ParameterSetDescription& );
  
 private:

  void propagate_to_ecal( reco::PreId&,
			  const reco::TrackRef&,
			  const edm::Handle<reco::PFClusterCollection>& );
  
  bool light_gsf_tracking( reco::PreId&,
			   const reco::TrackRef&,
			   const reco::ElectronSeed&,
			   const edm::EventSetup& );

  bool bdt_decision( reco::PreId& );
  
  edm::ESHandle<MagneticField> field_;
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  edm::EDGetTokenT<reco::PFClusterCollection> clusters_;
  std::string fitter_;
  std::string smoother_;
  std::string builder_;
  double threshold_;
  bool pass_through_;
  
};

#endif // RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedProducer_h
