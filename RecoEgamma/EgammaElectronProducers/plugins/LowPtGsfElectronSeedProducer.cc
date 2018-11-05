// -*- C++ -*-
//
// Package:    ElectronProducers
// Class:      LowPtGsfElectronSeedProducer
//
/**\class LowPtGsfElectronSeedProducer RecoEgamma/ElectronProducers/src/LowPtGsfElectronSeedProducer.cc
 Description: EDProducer of ElectronSeed objects
 Implementation:
     <Notes on implementation>
*/
// Original Author:  Robert Bainbridge

#include "LowPtGsfElectronSeedProducer.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TMath.h"

LowPtGsfElectronSeedProducer::LowPtGsfElectronSeedProducer( const edm::ParameterSet& iConfig ) :
  tracks_()
{
  tracks_ = consumes<reco::TrackCollection>( iConfig.getParameter<edm::InputTag>("tracks") );
  clusters_ = consumes<reco::PFClusterCollection>( iConfig.getParameter<edm::InputTag>("clusters") );
  produces<reco::ElectronSeedCollection>();
  produces<reco::PreIdCollection>();
}

LowPtGsfElectronSeedProducer::~LowPtGsfElectronSeedProducer() {}

void LowPtGsfElectronSeedProducer::beginRun( edm::Run const&, 
					     edm::EventSetup const& es ) 
{
  edm::ESHandle<MagneticField> field;
  es.get<IdealMagneticFieldRecord>().get(field);
  field_ = field->inTesla(GlobalPoint(0,0,0));
}

void LowPtGsfElectronSeedProducer::endRun( edm::Run const&, 
					   edm::EventSetup const& ) {}

void LowPtGsfElectronSeedProducer::produce( edm::Event& iEvent, 
					    const edm::EventSetup& iSetup ) 
{
  
  auto seeds = std::make_unique<reco::ElectronSeedCollection>();
  auto preids = std::make_unique<reco::PreIdCollection>();
  
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tracks_, tracks);
  
  edm::Handle<reco::PFClusterCollection> clusters;
  iEvent.getByToken(clusters_,clusters);

  for ( unsigned int itrk = 0; itrk < tracks.product()->size(); itrk++ ) {

    reco::TrackRef ref(tracks,itrk);
    if ( !(ref->quality(reco::TrackBase::qualityByName("highPurity"))) ) { continue; }

    //@@ BDT logic here to identify electron seeds based on tracks and other inputs

    reco::ElectronSeed seed(*(ref->seedRef()));
    seed.setCtfTrack(ref);
    seeds->push_back(seed);

    reco::PreId preid = create_preid(ref,clusters);
    preids->push_back(preid);

  }

  iEvent.put(std::move(seeds));
  iEvent.put(std::move(preids));

}

reco::PreId LowPtGsfElectronSeedProducer::create_preid( const reco::TrackRef& ref,
							const edm::Handle<reco::PFClusterCollection>& clusters )
{

  reco::PreId preid;
  preid.setTrack(ref);

  // Propagate track to ECAL surface
  float energy = sqrt( pow(0.000511,2.) + ref->outerMomentum().Mag2() );
  
  XYZTLorentzVector mom = XYZTLorentzVector( ref->outerMomentum().x(),
					     ref->outerMomentum().y(),
					     ref->outerMomentum().z(),
					     energy );
  XYZTLorentzVector pos = XYZTLorentzVector( ref->outerPosition().x(),
					     ref->outerPosition().y(),
					     ref->outerPosition().z(),
					     0. );
  
  BaseParticlePropagator particle( RawParticle(mom,pos), 0, 0, field_.z() );
  particle.setCharge(ref->charge());
  particle.propagateToEcalEntrance(false);
  
  if ( particle.getSuccess() != 0 ) { return preid; }

  // ECAL entry point for track
  GlobalPoint ecal_pos(particle.vertex().x(),
		       particle.vertex().y(),
		       particle.vertex().z());
  // Preshower limit
  bool below_ps = pow(ecal_pos.z(),2.) > pow(2.50746495928f,2.)*ecal_pos.perp2();

  // Store info 
  struct Info {
    reco::PFClusterRef clu = reco::PFClusterRef();
    float min_dr = 1.e6;
    float dr = 1.e6;
    float deta = 1.e6;
    float dphi = 1.e6;
    math::XYZPoint shower_pos = math::XYZPoint(0.,0.,0.);
  } info;

  // Iterate through clusters 
  for ( unsigned int iclu = 0; iclu < clusters.product()->size(); iclu++ ) {
    reco::PFClusterRef clu(clusters,iclu);

    // Correct ecal_pos for shower depth 
    double shower_depth = reco::PFCluster::getDepthCorrection(clu->correctedEnergy(),
							      below_ps,
							      false);
    GlobalPoint shower_pos = ecal_pos + 
      GlobalVector(particle.momentum().x(),
		   particle.momentum().y(),
		   particle.momentum().z()).unit() * shower_depth;

    // Determine deta, dphi, dr
    float deta = std::abs( clu->positionREP().eta() - shower_pos.eta() );
    float dphi = std::abs( clu->positionREP().phi() - shower_pos.phi() );
    if ( dphi > float(TMath::Pi())) { dphi -= float(TMath::TwoPi()); }
    float dr = std::sqrt( std::pow(dphi,2.f) + std::pow(deta,2.f) );
    
    // Find nearest cluster
    if ( dr < info.min_dr ) {
      info.min_dr = dr;
      info.clu = clu;
      info.dr = dr;
      info.deta = deta;
      info.dphi = dphi;
      info.shower_pos = shower_pos;
    }
  
  }

  // Populate PreId object
  math::XYZPoint point( ecal_pos.x(),
			ecal_pos.y(),
			ecal_pos.z() );
  float ep = info.clu->correctedEnergy() / std::sqrt( ref->innerMomentum().mag2() );
  preid.setECALMatchingProperties( info.clu,
				   point,
				   info.shower_pos,
				   info.deta,
				   info.dphi,
				   0.f, // chieta
				   0.f, // chiphi
				   0.f, // chi2
				   ep );
  
  return preid;
}

void LowPtGsfElectronSeedProducer::fillDescriptions( edm::ConfigurationDescriptions& descriptions ) 
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks",edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("clusters",edm::InputTag("particleFlowClusterECAL"));
  descriptions.add("lowPtGsfElectronSeeds",desc);
}
