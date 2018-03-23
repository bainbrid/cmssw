import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.MVAValueMapProducer_cfi import electronMVAVariableHelper
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("MVATrainingNtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# File with the ID variables to include in the Ntuplizer
mvaVariablesFile = "RecoEgamma/ElectronIdentification/python/Training/electron_variables.txt"

outputFile = "electron_ntuple.root"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/0293A280-B5F3-E711-8303-3417EBE33927.root'
        'file:/data_CMS/cms/rembser/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/50000/022E2036-A2D9-E711-9A8C-0CC47A13D2A4.root',
    )
)

process.ntuplizer = cms.EDAnalyzer('MVATrainingNtuplizer',
        # AOD case
        src                  = cms.InputTag('gedGsfElectrons'),
        vertices             = cms.InputTag('offlinePrimaryVertices'),
        pileup               = cms.InputTag('addPileupInfo'),
        # miniAOD case
        srcMiniAOD           = cms.InputTag('slimmedElectrons'),
        verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
        #
        variableDefinition   = cms.string(mvaVariablesFile),
        isMC                 = cms.bool(True)
        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.electronMVAVariableHelper = electronMVAVariableHelper
process.p = cms.Path(process.electronMVAVariableHelper * process.ntuplizer)
