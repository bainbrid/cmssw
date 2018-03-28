import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import electronMVAVariableHelper
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("ElectronMVANtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# File with the ID variables to include in the Ntuplizer
mvaVariablesFile = "RecoEgamma/ElectronIdentification/data/ElectronIDVariables.txt"

outputFile = "electron_ntuple.root"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # /DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM
        '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/0293A280-B5F3-E711-8303-3417EBE33927.root',

        ## /WToENu_M-100_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
        # '/store/mc/RunIIFall17MiniAOD/WToENu_M-100_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/02BF2F02-3EF9-E711-8343-002590E7E010.root',
        # '/store/mc/RunIIFall17MiniAOD/WToENu_M-100_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/AEE8E91B-C5F5-E711-AB79-1866DA7F9160.root',
        # '/store/mc/RunIIFall17MiniAOD/WToENu_M-100_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/DA0CCDF2-91F6-E711-A773-A4BF0112BDCC.root',
        # '/store/mc/RunIIFall17MiniAOD/WToENu_M-100_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/6638350B-3EF9-E711-8C5C-A0369F301924.root',
        # '/store/mc/RunIIFall17MiniAOD/WToENu_M-100_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/60000/729D2D36-8AEB-E711-8A52-FA163E341C2D.root',


        ## /WToENu_M-1000_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
        # '/store/mc/RunIIFall17MiniAOD/WToENu_M-1000_TuneCP5_13TeV-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/18EDE0D2-1CF5-E711-9C5E-0242AC130002.root',

        ## /QCD_Pt_80to120_TuneCP5_13TeV_pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
        # '/store/mc/RunIIFall17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/281B76D0-DCD7-E711-B799-24BE05C6E591.root',

        ## /TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM
        #'/store/mc/RunIIFall17MiniAOD/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/40000/468BB9EA-650D-E811-9122-7CD30AD08E7E.root',
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

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                 ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ntuplizer = cms.EDAnalyzer('ElectronMVANtuplizer',
        # AOD case
        src                  = cms.InputTag('gedGsfElectrons'),
        vertices             = cms.InputTag('offlinePrimaryVertices'),
        pileup               = cms.InputTag('addPileupInfo'),
        genParticles         = cms.InputTag('genParticles'),
        MET                  = cms.InputTag('slimmedMETs'),
        # miniAOD case
        srcMiniAOD           = cms.InputTag('slimmedElectrons'),
        verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
        genParticlesMiniAOD  = cms.InputTag('prunedGenParticles'),
        METMiniAOD           = cms.InputTag('slimmedMETs'),
        #
        eleMVAs             = cms.vstring(
                                          "egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp80",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wpLoose",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpHZZ",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp80",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpLoose",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp90",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80",
                                          "egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose",
                                          ),
        eleMVALabels        = cms.vstring(
                                          "Fall17noIsoV2wp80",
                                          "Fall17noIsoV2wpLoose",
                                          "Fall17noIsoV2wp90",
                                          "Fall17isoV2wpHZZ",
                                          "Fall17isoV2wp80",
                                          "Fall17isoV2wpLoose",
                                          "Fall17isoV2wp90",
                                          "Fall17noIsoV1wp90",
                                          "Fall17noIsoV1wp80",
                                          "Fall17noIsoV1wpLoose",
                                          "Fall17isoV1wp90",
                                          "Fall17isoV1wp80",
                                          "Fall17isoV1wpLoose",
                                          ),
        eleMVAValMaps        = cms.vstring(
                                           "electronMVAValueMapProducer:Fall17NoIsoV2Values",
                                           "electronMVAValueMapProducer:Fall17NoIsoV2RawValues",
                                           "electronMVAValueMapProducer:Fall17IsoV2Values",
                                           "electronMVAValueMapProducer:Fall17IsoV2RawValues",
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values",
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values",
                                           ),
        eleMVAValMapLabels   = cms.vstring(
                                           "Fall17NoIsoV2Vals",
                                           "Fall17NoIsoV2RawVals",
                                           "Fall17IsoV2Vals",
                                           "Fall17IsoV2RawVals",
                                           "Fall17IsoV1Vals",
                                           "Fall17NoIsoV1Vals",
                                           ),
        eleMVACats           = cms.vstring(
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories",
                                           ),
        eleMVACatLabels      = cms.vstring(
                                           "EleMVACats",
                                           ),
        #
        variableDefinition   = cms.string(mvaVariablesFile),
        isMC                 = cms.bool(True),
        deltaR               = cms.double(0.1),
        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.electronMVAVariableHelper = electronMVAVariableHelper
# process.p = cms.Path(process.egmGsfElectronIDSequence * process.electronMVAVariableHelper * process.ntuplizer)
process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntuplizer)

