import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.slimming.slimmedLowPtElectrons_cfi import slimmedLowPtElectrons
from PhysicsTools.PatAlgos.slimming.lowPtGsfLinks_cfi import lowPtGsfLinks
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cfi import lowPtGsfElectronID

# Clone to run ID V2 on patLowPtElectrons
patLowPtGsfElectronID = lowPtGsfElectronID.clone(
    electrons='patLowPtElectrons',
    rho='fixedGridRhoFastjetAll',
    ModelWeights=cms.vstring(['RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root']),
    Version=cms.string('V1'),
)

# Task
slimmedLowPtElectronsTask = cms.Task(
    lowPtGsfLinks,
    patLowPtGsfElectronID,
    slimmedLowPtElectrons,
)
