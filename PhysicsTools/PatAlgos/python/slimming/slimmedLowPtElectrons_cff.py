import FWCore.ParameterSet.Config as cms

# Clone with config for userFloat("ID")
from PhysicsTools.PatAlgos.slimming.slimmedLowPtElectrons_cfi import slimmedLowPtElectrons
slimmedLowPtElectrons_ = slimmedLowPtElectrons.clone(
    modifyElectrons = cms.bool(True),
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(
            cms.PSet(
                electron_config = cms.PSet(
                    ID = cms.InputTag("patLowPtGsfElectronID"),
                    electronSrc = cms.InputTag("patLowPtElectrons"),
                ),
                modifierName = cms.string('EGExtraInfoModifierFromFloatValueMaps'),
                photon_config = cms.PSet()
            ),
        )
    )
)

# Add userFloat("ID") only for run2_miniAOD_UL
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
run2_miniAOD_UL.toReplaceWith(slimmedLowPtElectrons,slimmedLowPtElectrons_)

# Links
from PhysicsTools.PatAlgos.slimming.lowPtGsfLinks_cfi import lowPtGsfLinks

# Task
slimmedLowPtElectronsTask = cms.Task(
    slimmedLowPtElectrons,
    lowPtGsfLinks,
)

# Clone to run ID V2 on patLowPtElectrons
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cfi import lowPtGsfElectronID
patLowPtGsfElectronID = lowPtGsfElectronID.clone(
    electrons='patLowPtElectrons',
    rho='fixedGridRhoFastjetAll',
    ModelWeights=cms.vstring(['RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Sept15.root']),
    Version=cms.string('V1'),
)

# Add ID module to task only for run2_miniAOD_UL
slimmedLowPtElectronsTask_ = cms.Task(slimmedLowPtElectronsTask.copy(),patLowPtGsfElectronID)
run2_miniAOD_UL.toReplaceWith(slimmedLowPtElectronsTask,slimmedLowPtElectronsTask_)
