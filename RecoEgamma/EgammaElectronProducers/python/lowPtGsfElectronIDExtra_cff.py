import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaElectronProducers.defaultLowPtGsfElectronIDExtra_cfi import defaultLowPtGsfElectronIDExtra

lowPtGsfElectronIDExtra = defaultLowPtGsfElectronIDExtra.clone(
    ModelNames = cms.vstring(['']),
    ModelWeights = cms.vstring([
            'RecoEgamma/ElectronIdentification/data/LowPtElectrons/RunII_Autumn18_LowPtElectrons_mva_id_2020Jun25_depth13_ntrees2000.xml.gz',  
            ]),
    ModelThresholds = cms.vdouble([-10.])
    )
