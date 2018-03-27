import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import electronMVAVariableHelper
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("MVATrainingNtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# File with the ID variables to include in the Ntuplizer
mvaVariablesFile = "RecoEgamma/ElectronIdentification/data/ElectronIDVariables.txt"

outputFile = "electron_ntuple.root"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/0293A280-B5F3-E711-8303-3417EBE33927.root'
    )
)

process.ntuplizer = cms.EDAnalyzer('ElectronMVATrainingNtuplizer',
        # AOD case
        src                  = cms.InputTag('gedGsfElectrons'),
        vertices             = cms.InputTag('offlinePrimaryVertices'),
        pileup               = cms.InputTag('addPileupInfo'),
        genParticles         = cms.InputTag('genParticles'),
        # miniAOD case
        srcMiniAOD           = cms.InputTag('slimmedElectrons'),
        verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
        genParticlesMiniAOD  = cms.InputTag('prunedGenParticles'),
        #
        variableDefinition   = cms.string(mvaVariablesFile),
        isMC                 = cms.bool(True),
        deltaR               = cms.double(0.1)
        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.electronMVAVariableHelper = electronMVAVariableHelper
process.p = cms.Path(process.electronMVAVariableHelper * process.ntuplizer)
