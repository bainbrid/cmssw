#include "RecoEgamma/EgammaElectronProducers/plugins/LowPtGsfElectronProducer.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/GsfElectronAlgo.h"
#include <iostream>

using namespace reco;

LowPtGsfElectronProducer::LowPtGsfElectronProducer( const edm::ParameterSet& cfg, 
						    const gsfAlgoHelpers::HeavyObjectCache* hoc )
  : GsfElectronBaseProducer(cfg,hoc)
{}

LowPtGsfElectronProducer::~LowPtGsfElectronProducer()
{}

void LowPtGsfElectronProducer::produce( edm::Event& event, const edm::EventSetup& setup )
{
  auto electrons = std::make_unique<GsfElectronCollection>();
  edm::Handle<reco::GsfElectronCoreCollection> coreElectrons;
  event.getByToken(inputCfg_.gsfElectronCores,coreElectrons);
  for ( unsigned int ii=0; ii < coreElectrons->size(); ++ii ) {
    const GsfElectronCoreRef ref = edm::Ref<GsfElectronCoreCollection>(coreElectrons,ii);
    GsfElectron* ele = new GsfElectron(ref);
    const GsfTrackRef& gsf = ref->gsfTrack();
    ele->setP4(GsfElectron::P4_FROM_SUPER_CLUSTER,
	       Candidate::LorentzVector(gsf->px(),gsf->py(),gsf->pz(),0.511E-3),
	       0,
	       true);
    LogTrace("GsfElectronAlgo") 
      << "[LowPtGsfElectronProducer::produce]"
      << " Constructed new electron with energy " 
      << ele->p4().e();
    electrons->push_back(*ele) ;
  }
  event.put(std::move(electrons));
}
