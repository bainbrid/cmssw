import FWCore.ParameterSet.Config as cms
# Low pT Electron value maps
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronSeedValueMaps_cfi import lowPtGsfElectronSeedValueMaps

# Low pT Electron value maps, rekeyed by reco::GsfElectron
rekeyLowPtGsfElectronSeedValueMaps = lowPtGsfElectronSeedValueMaps.clone(
    Rekey=cms.bool(True),
    gsfElectrons=cms.InputTag("lowPtGsfElectrons"),
    floatValueMaps=cms.VInputTag(["lowPtGsfElectronSeedValueMaps:unbiased",
                                  "lowPtGsfElectronSeedValueMaps:ptbiased"]),
)
