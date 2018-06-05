import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("PhotonMVANtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

outputFile = "photon_validation_ntuple.root"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/0293A280-B5F3-E711-8303-3417EBE33927.root'
    )
)

useAOD = False

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V0_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2p1_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V0_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V1_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V2_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V2p1_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff',
        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff',
                 ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.ntuplizer = cms.EDAnalyzer('PhotonMVANtuplizer',
        # AOD case
        src                  = cms.InputTag('gedGsfPhotons'),
        vertices             = cms.InputTag('offlinePrimaryVertices'),
        pileup               = cms.InputTag('addPileupInfo'),
        # miniAOD case
        srcMiniAOD           = cms.InputTag('slimmedPhotons'),
        verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
        #
        phoMVAs             = cms.vstring(
                                          ),
        phoMVALabels        = cms.vstring(
                                          ),
        phoMVAValMaps        = cms.vstring(
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV0Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig50nsV0Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig50nsV1Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig50nsV2Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig50nsV2p1Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Values",
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1p1Values",
                                           ),
        phoMVAValMapLabels   = cms.vstring(
                                           "Spring15NonTrig25nsV0",
                                           "Spring15NonTrig25nsV2",
                                           "Spring15NonTrig25nsV2p1",
                                           "Spring15NonTrig50nsV0",
                                           "Spring15NonTrig50nsV1",
                                           "Spring15NonTrig50nsV2",
                                           "Spring15NonTrig50nsV2p1",
                                           "Spring16NonTrigV1",
                                           "Fall17v1",
                                           "Fall17v1p1",
                                           ),
        phoMVACats           = cms.vstring(
                                           "photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1Categories",
                                           ),
        phoMVACatLabels      = cms.vstring(
                                           "PhoMVACats",
                                           ),
        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.p = cms.Path(process.egmPhotonIDSequence * process.ntuplizer)
