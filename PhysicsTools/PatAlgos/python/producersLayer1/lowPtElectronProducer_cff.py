import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.lowPtGsfLinks_cfi import lowPtGsfLinks

sourceElectrons = cms.InputTag("lowPtGsfElectrons")

lowPtElectronMatch = electronMatch.clone(
   src = sourceElectrons,
)

patLowPtElectrons = patElectrons.clone(

    # Input collection
    electronSource = sourceElectrons,

    # MC matching
    addGenMatch      = cms.bool(True),
    embedGenMatch    = cms.bool(True),
    genParticleMatch = cms.InputTag("lowPtElectronMatch"),

    # Electron ID
    addElectronID = cms.bool(True),
    electronIDSources = cms.PSet(
        unbiased = cms.InputTag("rekeyLowPtGsfElectronSeedValueMaps:unbiased"),
        ptbiased = cms.InputTag("rekeyLowPtGsfElectronSeedValueMaps:ptbiased"),
    ),

    # User data to add
    userData = cms.PSet(
        userClasses = cms.PSet( src = cms.VInputTag('') ),
        userFloats = cms.PSet( src = cms.VInputTag('') ),
        userInts = cms.PSet( src = cms.VInputTag('') ),
        userCands = cms.PSet( src = cms.VInputTag('lowPtGsfLinks') ),
        userFunctions = cms.vstring(),
        userFunctionLabels = cms.vstring()
    ),

    # Embedding of RECO/AOD items
    embedGsfElectronCore        = cms.bool(True),
    embedPFCandidate            = cms.bool(False),
    embedGsfTrack               = cms.bool(True),
    embedTrack                  = cms.bool(True),
    embedSuperCluster           = cms.bool(True),
    embedSeedCluster            = cms.bool(True),
    embedBasicClusters          = cms.bool(True),
    embedPreshowerClusters      = cms.bool(False),
    embedRecHits                = cms.bool(False),
    embedPflowSuperCluster      = cms.bool(False),
    embedPflowBasicClusters     = cms.bool(False),
    embedPflowPreshowerClusters = cms.bool(False),

    # Mini-Isolation 
    computeMiniIso    = cms.bool(True), # See: PhysicsTools/PatUtils/src/MiniIsolation.cc
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    miniIsoParamsE    = cms.vdouble(0.2, 0.2, 10.0, 0.0, 0.015, 0.015, 0.08, 0.0, 0.0),
    miniIsoParamsB    = cms.vdouble(0.2, 0.2, 10.0, 0.0, 0.000, 0.000, 0.00, 0.0, 0.0),

    # Miscellaneous flags
    useParticleFlow         = cms.bool(False),
    usePfCandidateMultiMap  = cms.bool(False),
    addMVAVariables         = cms.bool(False),
    addEfficiencies         = cms.bool(False),
    addResolutions          = cms.bool(False),
    embedHighLevelSelection = cms.bool(False),
    addPFClusterIso         = cms.bool(False),
    addPuppiIsolation       = cms.bool(False),
    isoDeposits             = cms.PSet(),
    userIsolation           = cms.PSet(),
    isolationValues         = cms.PSet(),
    isolationValuesNoPFId   = cms.PSet(),
)

makePatLowPtElectronsTask = cms.Task(
    lowPtGsfLinks,
    lowPtElectronMatch,
    patLowPtElectrons,
    )

makePatLowPtElectrons = cms.Sequence(makePatLowPtElectronsTask)


# Modifiers
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(patLowPtElectrons, embedTrack = False)
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17
(run2_miniAOD_80XLegacy | run2_miniAOD_94XFall17).toModify(patLowPtElectrons,
                                                           electronSource = "gedGsfElectrons",
                                                           genParticleMatch = "electronMatch"
                                                           )

# Schedule rekeying of seed BDT ValueMaps by reco::GsfElectron for run2_miniAOD_UL
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronSeedValueMaps_cff import rekeyLowPtGsfElectronSeedValueMaps
makePatLowPtElectronsTask_ = cms.Task(makePatLowPtElectronsTask.copy(),rekeyLowPtGsfElectronSeedValueMaps)
run2_miniAOD_UL.toReplaceWith(makePatLowPtElectronsTask,makePatLowPtElectronsTask_)
