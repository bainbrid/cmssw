import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaElectronProducers.defaultLowPtGsfElectronID_cfi import defaultLowPtGsfElectronID

lowPtGsfElectronID = defaultLowPtGsfElectronID.clone(
    ModelNames = cms.vstring(['']),
    ModelWeights = cms.vstring([
            'RecoEgamma/ElectronIdentification/data/LowPtElectrons/RunII_Autumn18_LowPtElectrons_mva_id.xml.gz',
            ]),
    ModelThresholds = cms.vdouble([-10.])
    )

from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
run2_miniAOD_UL.toModify(
    lowPtGsfElectronID,
    rho = "fixedGridRhoFastjetAll",
    ModelWeights = ["RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root"],
    ModelThresholds = [-99.],
    Version = "V1",
)

from Configuration.Eras.Modifier_bParking_cff import bParking
from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel
(bParking | run2_miniAOD_devel).toModify(
    lowPtGsfElectronID,
    ModelWeights = ["RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Nov28.root"],
)
