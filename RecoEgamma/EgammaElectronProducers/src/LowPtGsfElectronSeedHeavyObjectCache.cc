#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"
#include "RecoEgamma/EgammaTools/interface/GBRForestTools.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Reader.h"
#include <string>

namespace lowptgsfeleseed {

  ////////////////////////////////////////////////////////////////////////////////
  //
  HeavyObjectCache::HeavyObjectCache( const edm::ParameterSet& conf ) {
    edm::FileInPath unbiased( conf.getParameter<edm::FileInPath>("UnbiasedWeights") );
    edm::FileInPath ptbiased( conf.getParameter<edm::FileInPath>("PtBiasedWeights") );
    unbiased_ = GBRForestTools::createGBRForest(unbiased);
    ptbiased_ = GBRForestTools::createGBRForest(ptbiased);
    unbiasedThreshold_ = conf.getParameter<double>("UnbiasedThreshold");
    ptbiasedThreshold_ = conf.getParameter<double>("PtBiasedThreshold");
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  bool HeavyObjectCache::evalUnbiased( reco::PreId& ecal, reco::PreId& hcal ) const {

    float features[16] = { 
      (float)ecal.trackRef()->pt(),
      (float)ecal.trackRef()->eta(),
      (float)ecal.trackRef()->phi(),
      (float)ecal.trackRef()->p(),
      (float)ecal.trackRef()->charge(),
      (float)ecal.trackRef()->found(),
      (float)ecal.trackRef()->quality(reco::TrackBase::qualityByName("highPurity")),
      (float)sqrt( ecal.trackRef()->innerMomentum().Mag2() ),
      (float)sqrt( ecal.trackRef()->outerMomentum().Mag2() ),
      (float)ecal.trackRef()->normalizedChi2(),
      (float)ecal.geomMatching()[0],
      (float)ecal.geomMatching()[1],
      (float)ecal.eopMatch(),
      (float)ecal.dpt(),
      (float)ecal.chi2Ratio(),
      (float)ecal.gsfChi2()
    };
    float output = unbiased_->GetClassifier(features);
    bool pass = output > unbiasedThreshold_;
    ecal.setMVA( pass, output );
    return pass;

  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  bool HeavyObjectCache::evalPtBiased( reco::PreId& ecal, reco::PreId& hcal ) const {

    float features[16] = { 
      (float)ecal.trackRef()->pt(),
      (float)ecal.trackRef()->eta(),
      (float)ecal.trackRef()->phi(),
      (float)ecal.trackRef()->p(),
      (float)ecal.trackRef()->charge(),
      (float)ecal.trackRef()->found(),
      (float)ecal.trackRef()->quality(reco::TrackBase::qualityByName("highPurity")),
      (float)sqrt( ecal.trackRef()->innerMomentum().Mag2() ),
      (float)sqrt( ecal.trackRef()->outerMomentum().Mag2() ),
      (float)ecal.trackRef()->normalizedChi2(),
      (float)ecal.geomMatching()[0],
      (float)ecal.geomMatching()[1],
      (float)ecal.eopMatch(),
      (float)ecal.dpt(),
      (float)ecal.chi2Ratio(),
      (float)ecal.gsfChi2()
    };
    float output = ptbiased_->GetClassifier(features);
    bool pass = output > ptbiasedThreshold_;
    ecal.setMVA( pass, output );
    return pass;

  }

}
