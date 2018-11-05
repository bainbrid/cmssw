import FWCore.ParameterSet.Config as cms

lowPtGsfElectronSeeds = cms.EDProducer(
    "LowPtGsfElectronSeedProducer",
    tracks = cms.InputTag("generalTracks"),
    clusters = cms.InputTag("particleFlowClusterECAL"),
    )
