#ifndef RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedProducer_h
#define RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
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
    
  void withPfTracks( edm::Handle<reco::PFRecTrackCollection>&,
		     edm::Handle<reco::PFClusterCollection>& ecalClusters,
		     edm::Handle<reco::PFClusterCollection>& hcalClusters,
		     reco::ElectronSeedCollection&,
		     reco::PreIdCollection& ecal, 
		     reco::PreIdCollection& hcal,
		     const edm::EventSetup& );

  void withKfTracks( edm::Handle<reco::TrackCollection>&,
		     edm::Handle<reco::PFClusterCollection>& ecal,
		     reco::ElectronSeedCollection&,
		     reco::PreIdCollection&,
		     const edm::EventSetup& );

  void propagatePfTrackToEcal( reco::PreId&,
			       const reco::PFRecTrackRef&,
			       const edm::Handle<reco::PFClusterCollection>&, 
			       std::vector<int>& matchedClusters );
  
  void propagatePfTrackToHcal( reco::PreId&,
			       const reco::PFRecTrackRef&,
			       const edm::Handle<reco::PFClusterCollection>&, 
			       std::vector<int>& matchedClusters );
  
  void propagatePfTrackToCalo( reco::PreId&,
			       const reco::PFRecTrackRef&,
			       const edm::Handle<reco::PFClusterCollection>&, 
			       std::vector<int>& matchedClusters, 
			       bool ecal );
  
  void propagateKfTrackToEcal( reco::PreId&,
			       const reco::TrackRef&,
			       const edm::Handle<reco::PFClusterCollection>& ecalclusters );

  bool lightGsfTracking( reco::PreId&,
			 const reco::TrackRef&,
			 const reco::ElectronSeed&,
			 const edm::EventSetup& );
  
  edm::ESHandle<MagneticField> field_;
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  edm::EDGetTokenT<reco::PFClusterCollection> ecalClusters_;
  edm::EDGetTokenT<reco::PFClusterCollection> hcalClusters_;
  edm::EDGetTokenT<reco::PFRecTrackCollection> pfTracks_;
  std::string fitter_;
  std::string smoother_;
  std::string builder_;
  bool passThrough_;
  bool usePfTracks_;
  
  class Debug {
  public:
    int success = 0;
    int track = 0;
    int gsf = 0;
    void init() {
      success = 0;
      track = 0;
      gsf = 0;
    }
    void print() {
      std::cout << "Debug:" << std::endl
		<< " success: " << success << std::endl
  		<< " track: " << track << std::endl
		<< " gsf: " << gsf << std::endl;
    }
  } debug_;


};

#endif // RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedProducer_h
