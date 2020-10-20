import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.slimming.slimmedLowPtElectrons_cfi import slimmedLowPtElectrons

# Clone to run ID V2 on patLowPtElectrons
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cfi import lowPtGsfElectronID
patLowPtGsfElectronID = lowPtGsfElectronID.clone(
    electrons='patLowPtElectrons',
    rho='fixedGridRhoFastjetAll',
    ModelWeights=cms.vstring(['RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root']),
    Version=cms.string('V1'),
)

# Task
slimmedLowPtElectronsTask = cms.Task(
    slimmedLowPtElectrons,
    patLowPtGsfElectronID
)
