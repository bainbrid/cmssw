#ifndef LowPtGsfElectronSeedProducer_h
#define LowPtGsfElectronSeedProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"

class LowPtGsfElectronSeedProducer final : public edm::stream::EDProducer< edm::GlobalCache<lowpt::HeavyObjectCache> >
{
  
 public:
  
  explicit LowPtGsfElectronSeedProducer( const edm::ParameterSet&, 
					 const lowpt::HeavyObjectCache* );

  ~LowPtGsfElectronSeedProducer() override;
  
  static std::unique_ptr<lowpt::HeavyObjectCache> 
    initializeGlobalCache( const edm::ParameterSet& conf ) {
    return std::unique_ptr<lowpt::HeavyObjectCache>(new lowpt::HeavyObjectCache(conf));
  }
  
  static void globalEndJob( lowpt::HeavyObjectCache const* ) {}
  
  void beginLuminosityBlock( edm::LuminosityBlock const&, 
			     edm::EventSetup const& ) override;

  void produce( edm::Event&, const edm::EventSetup& ) override;
  
  static void fillDescriptions( edm::ConfigurationDescriptions& );
  
 private:

  typedef reco::PFClusterCollection PFClusters;
  
  bool propagate_to_ecal( reco::PreId&,
			  const reco::TrackRef&,
			  const edm::Handle<PFClusters>& );
  
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
  bool pass_through_;

  // For BDT
  float trk_pt;
  float trk_eta;
  float trk_phi;
  float trk_p;
  float trk_charge;
  float trk_nhits;
  float trk_high_purity;
  float trk_inp;
  float trk_outp;
  float trk_chi2red;
  float preid_trk_ecal_Deta;
  float preid_trk_ecal_Dphi;
  float preid_e_over_p;
  float preid_gsf_dpt;
  float preid_trk_gsf_chiratio;
  float preid_gsf_chi2red;
  
};

#endif // LowPtGsfElectronSeedProducer_h
