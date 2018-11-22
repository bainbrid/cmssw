import FWCore.ParameterSet.Config as cms

lowPtGsfElectronSeeds = cms.EDProducer(
    "LowPtGsfElectronSeedProducer",
    tracks = cms.InputTag("generalTracks"),
    pfTracks = cms.InputTag("lowPtGsfElePfTracks"),
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    Fitter = cms.string('GsfTrajectoryFitter_forPreId'),
    Smoother = cms.string('GsfTrajectorySmoother_forPreId'),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    UnbiasedWeights = cms.FileInPath('RecoEgamma/EgammaElectronProducers/data/BDT_formatted.xml.gz'),
    PtBiasedWeights = cms.FileInPath('RecoEgamma/EgammaElectronProducers/data/BDT_formatted.xml.gz'),
    UnbiasedThreshold = cms.double(-0.616),
    PtBiasedThreshold = cms.double(-0.616),
    PassThrough = cms.bool(False),
    UsePfTracks = cms.bool(True),
    )
