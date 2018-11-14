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
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkClonerImpl.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TMath.h"

LowPtGsfElectronSeedProducer::LowPtGsfElectronSeedProducer( const edm::ParameterSet& conf, 
							    const lowpt::HeavyObjectCache* ) :
  tracks_()
{
  tracks_ = consumes<reco::TrackCollection>( conf.getParameter<edm::InputTag>("tracks") );
  clusters_ = consumes<reco::PFClusterCollection>( conf.getParameter<edm::InputTag>("clusters") );
  fitter_ = conf.getParameter<std::string>("Fitter");
  smoother_ = conf.getParameter<std::string>("Smoother");
  builder_ = conf.getParameter<std::string>("TTRHBuilder");
  threshold_ = conf.getParameter<double>("Threshold");
  pass_through_ = conf.getParameter<bool>("PassThrough");
  produces<reco::ElectronSeedCollection>();
  produces<reco::PreIdCollection>();
}

LowPtGsfElectronSeedProducer::~LowPtGsfElectronSeedProducer() {}

void LowPtGsfElectronSeedProducer::beginLuminosityBlock( edm::LuminosityBlock const&, 
							 edm::EventSetup const& setup ) 
{
  setup.get<IdealMagneticFieldRecord>().get(field_);
}

void LowPtGsfElectronSeedProducer::produce( edm::Event& event, 
					    const edm::EventSetup& setup ) 
{

  auto seeds = std::make_unique<reco::ElectronSeedCollection>();
  auto preids = std::make_unique<reco::PreIdCollection>();
  
  edm::Handle<reco::TrackCollection> tracks;
  event.getByToken(tracks_, tracks);
  
  edm::Handle<reco::PFClusterCollection> clusters;
  event.getByToken(clusters_,clusters);

  for ( unsigned int itrk = 0; itrk < tracks.product()->size(); itrk++ ) {

    reco::TrackRef track_ref(tracks,itrk);
    if ( !(track_ref->quality(reco::TrackBase::qualityByName("highPurity"))) ) { continue; }

    // Create ElectronSeed 
    reco::ElectronSeed seed( *(track_ref->seedRef()) );
    seed.setCtfTrack(track_ref);
    
    // Create PreId
    reco::PreId preid;

    // Add track ref to PreId
    preid.setTrack(track_ref);

    // Add variables related to track-ECAL matching to PreId
    propagate_to_ecal(preid,track_ref,clusters);

    // Add variables related to GSF tracks to PreId
    light_gsf_tracking(preid,track_ref,seed,setup); 

    // Decision based on BDT 
    bool bdt = pass_through_ || bdt_decision(preid);

    // Store PreId
    preids->push_back(preid);

    // 
    if ( !bdt ) { continue; }
    
    seeds->push_back(seed);

  }

  event.put(std::move(seeds));
  event.put(std::move(preids));

}

bool LowPtGsfElectronSeedProducer::propagate_to_ecal( reco::PreId& preid,
						      const reco::TrackRef& track_ref,
						      const edm::Handle<reco::PFClusterCollection>& clusters )
{

  // Propagate 'electron' to ECAL surface
  float energy = sqrt( pow(0.000511,2.) + track_ref->outerMomentum().Mag2() );
  XYZTLorentzVector mom = XYZTLorentzVector( track_ref->outerMomentum().x(),
					     track_ref->outerMomentum().y(),
					     track_ref->outerMomentum().z(),
					     energy );
  XYZTLorentzVector pos = XYZTLorentzVector( track_ref->outerPosition().x(),
					     track_ref->outerPosition().y(),
					     track_ref->outerPosition().z(),
					     0. );
  math::XYZVector field(field_->inTesla(GlobalPoint(0,0,0)));
  BaseParticlePropagator particle( RawParticle(mom,pos), 0, 0, field.z() );
  particle.setCharge(track_ref->charge());
  particle.propagateToEcalEntrance(false);
  if ( particle.getSuccess() != 0 ) { return false; }
  
  // ECAL entry point for track
  GlobalPoint ecal_pos(particle.vertex().x(),
		       particle.vertex().y(),
		       particle.vertex().z());
  // Preshower limit
  bool below_ps = pow(ecal_pos.z(),2.) > pow(2.50746495928f,2.)*ecal_pos.perp2();

  // Store info for PreId
  struct Info {
    reco::PFClusterRef clu_ref = reco::PFClusterRef();
    float min_dr = 1.e6;
    float dr = 1.e6;
    float deta = 1.e6;
    float dphi = 1.e6;
    math::XYZPoint shower_pos = math::XYZPoint(0.,0.,0.);
  } info;
  
  // Iterate through ECAL clusters 
  for ( unsigned int iclu = 0; iclu < clusters.product()->size(); iclu++ ) {
    reco::PFClusterRef clu_ref(clusters,iclu);

    // Correct ecal_pos for shower depth 
    double shower_depth = reco::PFCluster::getDepthCorrection(clu_ref->correctedEnergy(),
							      below_ps,
							      false);
    GlobalPoint shower_pos = ecal_pos + 
      GlobalVector(particle.momentum().x(),
		   particle.momentum().y(),
		   particle.momentum().z()).unit() * shower_depth;

    // Determine deta, dphi, dr
    float deta = std::abs( clu_ref->positionREP().eta() - shower_pos.eta() );
    float dphi = std::abs( clu_ref->positionREP().phi() - shower_pos.phi() );
    if ( dphi > float(TMath::Pi())) { dphi -= float(TMath::TwoPi()); }
    float dr = std::sqrt( std::pow(dphi,2.f) + std::pow(deta,2.f) );
    
    // Find nearest ECAL cluster
    if ( dr < info.min_dr ) {
      info.min_dr = dr;
      info.clu_ref = clu_ref;
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
  float ep = info.clu_ref->correctedEnergy() / std::sqrt( track_ref->innerMomentum().mag2() );

  // Set PreId content
  preid.setECALMatchingProperties( info.clu_ref,
				   point,
				   info.shower_pos,
				   info.deta,
				   info.dphi,
				   0.f, // chieta
				   0.f, // chiphi
				   track_ref->normalizedChi2(), // chi2
				   ep );
  
  return true;
  
}

bool LowPtGsfElectronSeedProducer::light_gsf_tracking( reco::PreId& preid,
						       const reco::TrackRef& track_ref,
						       const reco::ElectronSeed& seed,
						       const edm::EventSetup& setup )
{
  
  edm::ESHandle<TrajectoryFitter> fitter;
  setup.get<TrajectoryFitter::Record>().get(fitter_,fitter);
  std::unique_ptr<TrajectoryFitter> fitter_ptr = fitter->clone();

  edm::ESHandle<TrajectorySmoother> smoother;
  setup.get<TrajectoryFitter::Record>().get(smoother_,smoother);
  std::unique_ptr<TrajectorySmoother> smoother_ptr;
  smoother_ptr.reset(smoother->clone());

  edm::ESHandle<TransientTrackingRecHitBuilder> builder;
  setup.get<TransientRecHitRecord>().get(builder_,builder);
  TkClonerImpl hit_cloner = static_cast<TkTransientTrackingRecHitBuilder const*>(builder.product())->cloner();
  fitter_ptr->setHitCloner(&hit_cloner);
  smoother_ptr->setHitCloner(&hit_cloner);

  Trajectory::ConstRecHitContainer hits;
  for ( unsigned int ihit = 0; ihit < track_ref->recHitsSize(); ++ihit ) {
    hits.push_back( track_ref->recHit(ihit)->cloneSH() );
  }
  
  GlobalVector gv( track_ref->innerMomentum().x(),
		   track_ref->innerMomentum().y(),
		   track_ref->innerMomentum().z() );
  GlobalPoint gp( track_ref->innerPosition().x(),
		  track_ref->innerPosition().y(),
		  track_ref->innerPosition().z() );

  GlobalTrajectoryParameters gtps( gp,
				   gv,
				   track_ref->charge(),
				   &*field_ );

  TrajectoryStateOnSurface tsos( gtps,
				 track_ref->innerStateCovariance(),
				 *hits[0]->surface() );

  // Track fitted and smoothed under electron hypothesis
  Trajectory traj1 = fitter_ptr->fitOne( seed, hits, tsos );
  if ( !traj1.isValid() ) { return false; }
  Trajectory traj2 = smoother_ptr->trajectory(traj1);
  if ( !traj2.isValid() ) {  return false; }
      
  // Set PreId content
  float chi2_ratio = traj2.chiSquared() / track_ref->chi2();
  float gsf_reduced_chi2 = chi2_ratio * track_ref->normalizedChi2();
  float pt_out = traj2.firstMeasurement().updatedState().globalMomentum().perp();
  float pt_in = traj2.lastMeasurement().updatedState().globalMomentum().perp();
  float gsf_dpt = ( pt_in > 0 ) ? fabs( pt_out - pt_in ) / pt_in : 0.;
  preid.setTrackProperties(gsf_reduced_chi2,chi2_ratio,gsf_dpt);

  return true;

}
  
bool LowPtGsfElectronSeedProducer::bdt_decision( reco::PreId& preid )
{
  
  trk_pt = preid.trackRef()->pt();
  trk_eta = preid.trackRef()->eta();
  trk_phi = preid.trackRef()->phi();
  trk_p = preid.trackRef()->p();
  trk_charge = preid.trackRef()->charge();
  trk_nhits = preid.trackRef()->found();
  trk_high_purity = preid.trackRef()->quality(reco::TrackBase::qualityByName("highPurity"));
  trk_inp = sqrt( preid.trackRef()->innerMomentum().Mag2() );
  trk_outp = sqrt( preid.trackRef()->outerMomentum().Mag2() );
  trk_chi2red = preid.trackRef()->normalizedChi2();
  preid_trk_ecal_Deta = preid.geomMatching()[0];
  preid_trk_ecal_Dphi = preid.geomMatching()[1];
  preid_e_over_p = preid.eopMatch();
  preid_gsf_dpt = preid.dpt();
  preid_trk_gsf_chiratio = preid.chi2Ratio();
  preid_gsf_chi2red = preid.gsfChi2();
  
  float features[16] = { 
    trk_pt,
    trk_eta,
    trk_phi,
    trk_p,
    trk_charge,
    trk_nhits,
    trk_high_purity,
    trk_inp,
    trk_outp,
    trk_chi2red,
    preid_trk_ecal_Deta,
    preid_trk_ecal_Dphi,
    preid_e_over_p,
    preid_gsf_dpt,
    preid_trk_gsf_chiratio,
    preid_gsf_chi2red
  };
  float output = globalCache()->gbr->GetClassifier(features);
  preid.setMVA( output > threshold_, output );
  
  return output > threshold_;
  
}

void LowPtGsfElectronSeedProducer::fillDescription( edm::ParameterSetDescription& desc ) 
{
  desc.add<edm::InputTag>("tracks",edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("clusters",edm::InputTag("particleFlowClusterECAL"));
  desc.add<std::string>("Fitter","GsfTrajectoryFitter_forPreId");
  desc.add<std::string>("Smoother","GsfTrajectorySmoother_forPreId");
  desc.add<std::string>("TTRHBuilder","WithAngleAndTemplate");
  desc.add<std::string>("Weights","RecoEgamma/EgammaElectronProducers/data/BDT.xml");
  desc.add<double>("Threshold",-0.616);
  desc.add<bool>("PassThrough",false);
}
