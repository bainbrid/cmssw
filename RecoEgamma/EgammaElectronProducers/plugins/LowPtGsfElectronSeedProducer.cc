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
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
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

//////////////////////////////////////////////////////////////////////////////////////////
//
LowPtGsfElectronSeedProducer::LowPtGsfElectronSeedProducer( const edm::ParameterSet& conf, 
							    const lowptgsfeleseed::HeavyObjectCache* ) :
  tracks_(),
  ecalClusters_(),
  hcalClusters_(),
  pfTracks_(),
  fitter_(""),
  smoother_(""),
  builder_(""),
  passThrough_(false),
  usePfTracks_(false)
{
  tracks_ = consumes<reco::TrackCollection>( conf.getParameter<edm::InputTag>("tracks") );
  ecalClusters_ = consumes<reco::PFClusterCollection>( conf.getParameter<edm::InputTag>("ecalClusters") );
  hcalClusters_ = consumes<reco::PFClusterCollection>( conf.getParameter<edm::InputTag>("hcalClusters") );
  pfTracks_ = consumes<reco::PFRecTrackCollection>( conf.getParameter<edm::InputTag>("pfTracks") );
  fitter_ = conf.getParameter<std::string>("Fitter");
  smoother_ = conf.getParameter<std::string>("Smoother");
  builder_ = conf.getParameter<std::string>("TTRHBuilder");
  passThrough_ = conf.getParameter<bool>("PassThrough");
  usePfTracks_ = conf.getParameter<bool>("UsePfTracks");
  produces<reco::ElectronSeedCollection>();
  produces<reco::PreIdCollection>();
  produces<reco::PreIdCollection>("HCAL");
}

//////////////////////////////////////////////////////////////////////////////////////////
//
LowPtGsfElectronSeedProducer::~LowPtGsfElectronSeedProducer() {}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::beginLuminosityBlock( edm::LuminosityBlock const&, 
							 edm::EventSetup const& setup ) 
{
  setup.get<IdealMagneticFieldRecord>().get(field_);
}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::produce( edm::Event& event, 
					    const edm::EventSetup& setup ) 
{
  
  // Debug info
  debug_.init();

  auto seeds = std::make_unique<reco::ElectronSeedCollection>();
  auto ecalPreIds = std::make_unique<reco::PreIdCollection>();
  auto hcalPreIds = std::make_unique<reco::PreIdCollection>();
    
  edm::Handle<reco::PFClusterCollection> ecalClusters;
  event.getByToken(ecalClusters_,ecalClusters);

  edm::Handle<reco::PFClusterCollection> hcalClusters;
  if ( usePfTracks_ ) { event.getByToken(hcalClusters_,hcalClusters); }

  edm::Handle<reco::TrackCollection> kfTracks;
  event.getByToken(tracks_, kfTracks);

  edm::Handle<reco::PFRecTrackCollection> pfTracks;
  if ( usePfTracks_ ) { event.getByToken(pfTracks_,pfTracks); }

  if ( usePfTracks_ ) {  
    withPfTracks(pfTracks,ecalClusters,hcalClusters,*seeds,*ecalPreIds,*hcalPreIds,setup); 
  } else { 
    withKfTracks(kfTracks,ecalClusters,*seeds,*ecalPreIds,setup);
  }

  event.put(std::move(seeds));
  event.put(std::move(ecalPreIds));
  event.put(std::move(hcalPreIds),"HCAL");

  debug_.print();

}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::withPfTracks( edm::Handle<reco::PFRecTrackCollection>& pfTracks,
						 edm::Handle<reco::PFClusterCollection>& ecalClusters,
						 edm::Handle<reco::PFClusterCollection>& hcalClusters,
						 reco::ElectronSeedCollection& seeds,
						 reco::PreIdCollection& ecalPreIds,
						 reco::PreIdCollection& hcalPreIds,
						 const edm::EventSetup& setup ) 
{

  std::vector<int> matchedEcalClusters;
  std::vector<int> matchedHcalClusters;

  for ( unsigned int itrk = 0; itrk < pfTracks.product()->size(); itrk++ ) {

    reco::PFRecTrackRef pfTrackRef(pfTracks,itrk);
    reco::TrackRef trackRef = pfTrackRef->trackRef();
    if ( !(trackRef->quality(reco::TrackBase::qualityByName("highPurity"))) ) { continue; }
    debug_.track++;

    // Create ElectronSeed 
    reco::ElectronSeed seed( *(trackRef->seedRef()) );
    seed.setCtfTrack(trackRef);
    
    // Create PreIds
    reco::PreId ecalPreId;
    reco::PreId hcalPreId;

    // Add track ref to PreId
    ecalPreId.setTrack(trackRef);
    hcalPreId.setTrack(trackRef);

    // Add variables related to track-ECAL and track-HCAL matching to PreId
    propagatePfTrackToEcal(ecalPreId,pfTrackRef,ecalClusters,matchedEcalClusters);
    propagatePfTrackToHcal(hcalPreId,pfTrackRef,hcalClusters,matchedHcalClusters);

    // Add variables related to GSF tracks to PreId
    lightGsfTracking(ecalPreId,trackRef,seed,setup); 

    // Decision based on BDT 
    bool unbiased = globalCache()->evalUnbiased(ecalPreId,hcalPreId);
    bool ptbiased = globalCache()->evalPtBiased(ecalPreId,hcalPreId);
    bool result = passThrough_ || ( unbiased || ptbiased );

    // Store PreId
    ecalPreIds.push_back(ecalPreId);
    hcalPreIds.push_back(hcalPreId);

    // If fails BDT, do not store seed
    if ( !result ) { continue; }
    
    seeds.push_back(seed);

  }

}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::withKfTracks( edm::Handle<reco::TrackCollection>& kfTracks,
						 edm::Handle<reco::PFClusterCollection>& clusters,
						 reco::ElectronSeedCollection& seeds,
						 reco::PreIdCollection& preIds,
						 const edm::EventSetup& setup ) 
{

  for ( unsigned int itrk = 0; itrk < kfTracks.product()->size(); itrk++ ) {

    reco::TrackRef trackRef(kfTracks,itrk);
    if ( !(trackRef->quality(reco::TrackBase::qualityByName("highPurity"))) ) { continue; }
    debug_.track++;

    // Create ElectronSeed 
    reco::ElectronSeed seed( *(trackRef->seedRef()) );
    seed.setCtfTrack(trackRef);
    
    // Create PreId
    reco::PreId preId;

    // Add track ref to PreId
    preId.setTrack(trackRef);

    // Add variables related to track-ECAL matching to PreId
    propagateKfTrackToEcal(preId,trackRef,clusters);

    // Add variables related to GSF tracks to PreId
    lightGsfTracking(preId,trackRef,seed,setup); 

    // Decision based on BDT 
    reco::PreId hcal;
    bool unbiased = globalCache()->evalUnbiased(preId,hcal);
    bool ptbiased = globalCache()->evalPtBiased(preId,hcal);
    bool result = passThrough_ || ( unbiased || ptbiased );

    // Store PreId
    preIds.push_back(preId);

    // If fails BDT, do not store seed
    if ( !result ) { continue; }
    
    seeds.push_back(seed);

  }

}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::propagatePfTrackToEcal( reco::PreId& preId,
							   const reco::PFRecTrackRef& pfTrackRef,
							   const edm::Handle<reco::PFClusterCollection>& clusters,
							   std::vector<int>& matchedClusters )
{
  propagatePfTrackToCalo(preId,pfTrackRef,clusters,matchedClusters,true);
}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::propagatePfTrackToHcal( reco::PreId& preId,
							   const reco::PFRecTrackRef& pfTrackRef,
							   const edm::Handle<reco::PFClusterCollection>& clusters,
							   std::vector<int>& matchedClusters )
{
  propagatePfTrackToCalo(preId,pfTrackRef,clusters,matchedClusters,false);
}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::propagatePfTrackToCalo( reco::PreId& preId,
							   const reco::PFRecTrackRef& pfTrackRef,
							   const edm::Handle<reco::PFClusterCollection>& clusters,
							   std::vector<int>& matchedClusters,
							   bool ecal )
{

  // Store info for PreId
  struct Info {
    reco::PFClusterRef cluRef = reco::PFClusterRef();
    float dr2min = 1.e6;
    float deta = 1.e6;
    float dphi = 1.e6;
    math::XYZPoint showerPos = math::XYZPoint(0.,0.,0.);
  } info;
  
  // Find closest "seed cluster" to KF track extrapolated to ECAL (or HCAL)
  reco::PFTrajectoryPoint point;
  if ( ecal ) { point = pfTrackRef->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::ECALShowerMax); }
  else        { point = pfTrackRef->extrapolatedPoint(reco::PFTrajectoryPoint::LayerType::HCALEntrance); }

  if ( point.isValid() ) {
    if ( ecal ) debug_.success++;

    Info info;
    for ( unsigned int iclu = 0; iclu < clusters.product()->size(); iclu++ ) {

      if ( std::find( matchedClusters.begin(), matchedClusters.end(), iclu ) == matchedClusters.end() ) {
	reco::PFClusterRef cluRef(clusters,iclu);

	// Determine deta, dphi, dr
	float deta = cluRef->positionREP().eta() - point.positionREP().eta();
	float dphi = reco::deltaPhi( cluRef->positionREP().phi(), point.positionREP().phi() );
	float dr2 = reco::deltaR2( cluRef->positionREP(), point.positionREP() );

	if ( dr2 < info.dr2min ) {
	  info.dr2min = dr2;
	  info.cluRef = cluRef;
	  info.deta = deta;
	  info.dphi = dphi;
	  info.showerPos = point.position();
	}

      }
    }

    // Set PreId content if match found
    if ( info.dr2min < 1.e5 ) { 
      float ep = info.cluRef->correctedEnergy() / std::sqrt( pfTrackRef->trackRef()->innerMomentum().mag2() );
      preId.setECALMatchingProperties( info.cluRef,
				       point.position(), // ECAL or HCAL surface
				       info.showerPos, // 
				       info.deta,
				       info.dphi,
				       0.f, // chieta
				       0.f, // chiphi
				       pfTrackRef->trackRef()->normalizedChi2(), // chi2
				       ep );
    }

  } // clusters

}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::propagateKfTrackToEcal( reco::PreId& preId,
							   const reco::TrackRef& trackRef,
							   const edm::Handle<reco::PFClusterCollection>& ecalClusters )
{

  // Propagate 'electron' to ECAL surface
  float energy = sqrt( pow(0.000511,2.) + trackRef->outerMomentum().Mag2() );
  XYZTLorentzVector mom = XYZTLorentzVector( trackRef->outerMomentum().x(),
					     trackRef->outerMomentum().y(),
					     trackRef->outerMomentum().z(),
					     energy );
  XYZTLorentzVector pos = XYZTLorentzVector( trackRef->outerPosition().x(),
					     trackRef->outerPosition().y(),
					     trackRef->outerPosition().z(),
					     0. );
  math::XYZVector field(field_->inTesla(GlobalPoint(0,0,0)));
  BaseParticlePropagator particle( RawParticle(mom,pos), 0, 0, field.z() );
  particle.setCharge(trackRef->charge());
  particle.propagateToEcalEntrance(false);
  if ( particle.getSuccess() == 0 ) { return; }
  debug_.success++;
  
  // ECAL entry point for track
  GlobalPoint ecal_pos(particle.vertex().x(),
		       particle.vertex().y(),
		       particle.vertex().z());
  // Preshower limit
  bool below_ps = pow(ecal_pos.z(),2.) > pow(2.50746495928f,2.)*ecal_pos.perp2();

  // Store info for PreId
  struct Info {
    reco::PFClusterRef cluRef = reco::PFClusterRef();
    float dr2min = 1.e6;
    float deta = 1.e6;
    float dphi = 1.e6;
    math::XYZPoint showerPos = math::XYZPoint(0.,0.,0.);
  } info;
  
  // Iterate through ECAL clusters 
  for ( unsigned int iclu = 0; iclu < ecalClusters.product()->size(); iclu++ ) {
    reco::PFClusterRef cluRef(ecalClusters,iclu);

    // Correct ecal_pos for shower depth 
    double shower_depth = reco::PFCluster::getDepthCorrection(cluRef->correctedEnergy(),
							      below_ps,
							      false);
    GlobalPoint showerPos = ecal_pos + 
      GlobalVector(particle.momentum().x(),
		   particle.momentum().y(),
		   particle.momentum().z()).unit() * shower_depth;

    // Determine deta, dphi, dr
    float deta = std::abs( cluRef->positionREP().eta() - showerPos.eta() );
    float dphi = std::abs( reco::deltaPhi( cluRef->positionREP().phi(), showerPos.phi() ));
    float dr2 = reco::deltaR2( cluRef->positionREP(), showerPos );

    // Find nearest ECAL cluster
    if ( dr2 < info.dr2min ) {
      info.dr2min = dr2;
      info.cluRef = cluRef;
      info.deta = deta;
      info.dphi = dphi;
      info.showerPos = showerPos;
    }
  
  }

  // Populate PreId object
  math::XYZPoint point( ecal_pos.x(),
			ecal_pos.y(),
			ecal_pos.z() );

  // Set PreId content
  preId.setECALMatchingProperties( info.cluRef,
				   point,
				   info.showerPos,
				   info.deta,
				   info.dphi,
				   0.f, // chieta
				   0.f, // chiphi
				   trackRef->normalizedChi2(), // chi2
				   info.cluRef->correctedEnergy() / std::sqrt( trackRef->innerMomentum().mag2() ) ); // E/p
  
}

//////////////////////////////////////////////////////////////////////////////////////////
//
bool LowPtGsfElectronSeedProducer::lightGsfTracking( reco::PreId& preId,
						     const reco::TrackRef& trackRef,
						     const reco::ElectronSeed& seed,
						     const edm::EventSetup& setup )
{
  
  edm::ESHandle<TrajectoryFitter> fitter;
  setup.get<TrajectoryFitter::Record>().get(fitter_,fitter);
  std::unique_ptr<TrajectoryFitter> fitterPtr = fitter->clone();

  edm::ESHandle<TrajectorySmoother> smoother;
  setup.get<TrajectoryFitter::Record>().get(smoother_,smoother);
  std::unique_ptr<TrajectorySmoother> smootherPtr;
  smootherPtr.reset(smoother->clone());

  edm::ESHandle<TransientTrackingRecHitBuilder> builder;
  setup.get<TransientRecHitRecord>().get(builder_,builder);
  TkClonerImpl hitCloner = static_cast<TkTransientTrackingRecHitBuilder const*>(builder.product())->cloner();
  fitterPtr->setHitCloner(&hitCloner);
  smootherPtr->setHitCloner(&hitCloner);

  Trajectory::ConstRecHitContainer hits;
  for ( unsigned int ihit = 0; ihit < trackRef->recHitsSize(); ++ihit ) {
    hits.push_back( trackRef->recHit(ihit)->cloneSH() );
  }
  
  GlobalVector gv( trackRef->innerMomentum().x(),
		   trackRef->innerMomentum().y(),
		   trackRef->innerMomentum().z() );
  GlobalPoint gp( trackRef->innerPosition().x(),
		  trackRef->innerPosition().y(),
		  trackRef->innerPosition().z() );

  GlobalTrajectoryParameters gtps( gp,
				   gv,
				   trackRef->charge(),
				   &*field_ );

  TrajectoryStateOnSurface tsos( gtps,
				 trackRef->innerStateCovariance(),
				 *hits[0]->surface() );

  // Track fitted and smoothed under electron hypothesis
  Trajectory traj1 = fitterPtr->fitOne( seed, hits, tsos );
  if ( !traj1.isValid() ) { return false; }
  Trajectory traj2 = smootherPtr->trajectory(traj1);
  if ( !traj2.isValid() ) {  return false; }

  debug_.gsf++;
      
  // Set PreId content
  float chi2Ratio = traj2.chiSquared() / trackRef->chi2();
  float gsfReducedChi2 = chi2Ratio * trackRef->normalizedChi2();
  float ptOut = traj2.firstMeasurement().updatedState().globalMomentum().perp();
  float ptIn = traj2.lastMeasurement().updatedState().globalMomentum().perp();
  float gsfDpt = ( ptIn > 0 ) ? fabs( ptOut - ptIn ) / ptIn : 0.;
  preId.setTrackProperties(gsfReducedChi2,chi2Ratio,gsfDpt);

  return true;

}

//////////////////////////////////////////////////////////////////////////////////////////
//
void LowPtGsfElectronSeedProducer::fillDescription( edm::ParameterSetDescription& desc ) 
{
  desc.add<edm::InputTag>("tracks",edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("pfTracks",edm::InputTag("lowPtGsfElePfTracks"));
  desc.add<edm::InputTag>("ecalClusters",edm::InputTag("particleFlowClusterECAL"));
  desc.add<edm::InputTag>("hcalClusters",edm::InputTag("particleFlowClusterHCAL"));
  desc.add<std::string>("Fitter","GsfTrajectoryFitter_forPreId");
  desc.add<std::string>("Smoother","GsfTrajectorySmoother_forPreId");
  desc.add<std::string>("TTRHBuilder","WithAngleAndTemplate");
  desc.add<edm::FileInPath>("UnbiasedWeights",edm::FileInPath("RecoEgamma/EgammaElectronProducers/data/BDT.xml"));
  desc.add<edm::FileInPath>("PtBiasedWeights",edm::FileInPath("RecoEgamma/EgammaElectronProducers/data/BDT.xml"));
  desc.add<double>("PtBiasedThreshold",1.);
  desc.add<double>("UnbiasedThreshold",1.);
  desc.add<bool>("PassThrough",false);
  desc.add<bool>("UsePfTracks",false);
}
