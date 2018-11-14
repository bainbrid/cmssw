import FWCore.ParameterSet.Config as cms

lowPtGsfElectronSeeds = cms.EDProducer(
    "LowPtGsfElectronSeedProducer",
    tracks = cms.InputTag("generalTracks"),
    clusters = cms.InputTag("particleFlowClusterECAL"),
    Fitter = cms.string('GsfTrajectoryFitter_forPreId'),
    Smoother = cms.string('GsfTrajectorySmoother_forPreId'),
    TTRHBuilder = cms.string('WithAngleAndTemplate'),
    Weights = cms.string('RecoEgamma/EgammaElectronProducers/data/BDT.xml'),
    Threshold = cms.double(-0.616),
    PassThrough = cms.bool(False),
    )
