#ifndef RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedHeavyObjectCache_h
#define RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedHeavyObjectCache_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include <memory>

namespace reco { class PreId; }

namespace lowptgsfeleseed {
  class HeavyObjectCache {

  public:

    HeavyObjectCache( const edm::ParameterSet& );

    bool evalUnbiased( reco::PreId& ecal, reco::PreId& hcal ) const;
    bool evalPtBiased( reco::PreId& ecal, reco::PreId& hcal ) const;

  private:

    std::unique_ptr<const GBRForest> unbiased_;
    std::unique_ptr<const GBRForest> ptbiased_;
    float unbiasedThreshold_;
    float ptbiasedThreshold_;

  };
}

#endif // RecoEgamma_EgammaElectronProducers_LowPtGsfElectronSeedHeavyObjectCache_h
