
import FWCore.ParameterSet.Config as cms

lowPtGsfElectronSuperClusters = cms.EDProducer(
    "LowPtGsfElectronSCProducer",
    gsfPfRecTracks = cms.InputTag("lowPtGsfElePfGsfTracks"),
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
)
