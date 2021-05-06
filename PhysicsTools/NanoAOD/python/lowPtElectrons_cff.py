import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *

################################################################################
# Modules
################################################################################

from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier106XUL

_lowPtRegressionModifier = regressionModifier106XUL.clone(
    modifierName = 'EGRegressionModifierV3',
    rhoTag = 'fixedGridRhoFastjetAll',
    eleRegs = dict(
        ecalOnlyMean = dict(
            lowEtHighEtBoundary = 20.,
            ebLowEtForestName = ":lowPtElectron_eb_ecalOnly_05To50_mean",
            ebHighEtForestName = ":lowPtElectron_eb_ecalOnly_05To50_mean",
            eeLowEtForestName = ":lowPtElectron_ee_ecalOnly_05To50_mean",
            eeHighEtForestName = ":lowPtElectron_ee_ecalOnly_05To50_mean",
            ),
        ecalOnlySigma = dict(
            lowEtHighEtBoundary = 20.,
            ebLowEtForestName = ":lowPtElectron_eb_ecalOnly_05To50_sigma",
            ebHighEtForestName = ":lowPtElectron_eb_ecalOnly_05To50_sigma",
            eeLowEtForestName = ":lowPtElectron_ee_ecalOnly_05To50_sigma",
            eeHighEtForestName = ":lowPtElectron_ee_ecalOnly_05To50_sigma",
            ),
        epComb = dict(
            ecalTrkRegressionConfig = dict(
                lowEtHighEtBoundary = 20.,
                ebLowEtForestName = ":lowPtElectron_eb_ecalTrk_05To50_mean",
                ebHighEtForestName = ":lowPtElectron_eb_ecalTrk_05To50_mean",
                eeLowEtForestName = ":lowPtElectron_ee_ecalTrk_05To50_mean",
                eeHighEtForestName = ":lowPtElectron_ee_ecalTrk_05To50_mean",
                ),
            ecalTrkRegressionUncertConfig = dict(
                lowEtHighEtBoundary = 20.,
                ebLowEtForestName = ":lowPtElectron_eb_ecalTrk_05To50_sigma",
                ebHighEtForestName = ":lowPtElectron_eb_ecalTrk_05To50_sigma",
                eeLowEtForestName = ":lowPtElectron_ee_ecalTrk_05To50_sigma",
                eeHighEtForestName = ":lowPtElectron_ee_ecalTrk_05To50_sigma",
                ),
        )
    ),
)

import PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi

updatedLowPtElectrons = cms.EDProducer(
    "PATElectronUpdater",
    src = cms.InputTag("slimmedLowPtElectrons"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    computeMiniIso = cms.bool(True),
    fixDxySign = cms.bool(True),
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    miniIsoParamsB = PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi.patElectrons.miniIsoParamsB,
    miniIsoParamsE = PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi.patElectrons.miniIsoParamsE,
    applyEnergyRegression = cms.bool(True),
    regressionConfig = _lowPtRegressionModifier,
)

#lowPtPATElectronID = cms.EDProducer(
#    "LowPtGsfElectronIDProducer",
#    electrons = cms.InputTag("updatedLowPtElectrons"),
#    unbiased = cms.InputTag(""),
#    rho = cms.InputTag("fixedGridRhoFastjetAll"),
#    ModelNames = cms.vstring(['']),
#    ModelWeights = cms.vstring([
#        'RecoEgamma/ElectronIdentification/data/LowPtElectrons/LowPtElectrons_ID_2020Nov28.root',
#    ]),
#    ModelThresholds = cms.vdouble([-99.]),
#    PassThrough = cms.bool(False),
#    MinPtThreshold = cms.double(0.5),
#    MaxPtThreshold = cms.double(15.),
#    Version = cms.string('V1'),
#)

isoForLowPtEle = cms.EDProducer(
    "EleIsoValueMapProducer",
    src = cms.InputTag("updatedLowPtElectrons"),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    rho_PFIso = cms.InputTag("fixedGridRhoFastjetAll"),
    EAFile_MiniIso = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
    EAFile_PFIso = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
)

updatedLowPtElectronsWithUserData = cms.EDProducer(
    "PATElectronUserDataEmbedder",
    src = cms.InputTag("updatedLowPtElectrons"),
    userFloats = cms.PSet(
        #ID = cms.InputTag("lowPtPATElectronID"),
        miniIsoChg = cms.InputTag("isoForLowPtEle:miniIsoChg"),
        miniIsoAll = cms.InputTag("isoForLowPtEle:miniIsoAll"),
    ),
    userIntFromBools = cms.PSet(),
    userInts = cms.PSet(),
    userCands = cms.PSet(),
)

finalLowPtElectrons = cms.EDFilter(
    "PATElectronRefSelector",
    src = cms.InputTag("updatedLowPtElectronsWithUserData"),
    cut = cms.string("pt > 1. && electronID('ID') > -0.25"),
)

################################################################################
# electronTable 
################################################################################

lowPtElectronTable = cms.EDProducer(
    "SimpleCandidateFlatTableProducer",
    src = cms.InputTag("finalLowPtElectrons"),
    cut = cms.string(""),
    name= cms.string("LowPtElectron"),
    doc = cms.string("slimmedLowPtElectrons after basic selection (" + finalLowPtElectrons.cut.value()+")"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the electrons
    variables = cms.PSet(
        # Basic variables
        CandVars,
        # BDT scores and WPs
        ID = Var("electronID('ID')",float,doc="ID, BDT (raw) score"),
        unbiased = Var("electronID('unbiased')",float,doc="ElectronSeed, pT- and dxy- agnostic BDT (raw) score"),
        ptbiased = Var("electronID('ptbiased')",float,doc="ElectronSeed, pT- and dxy- dependent BDT (raw) score"),
        # Isolation
        miniRelIso_chg = Var("miniPFIsolation().chargedHadronIso()/pt",float,
                             doc="mini PF rel. iso. (charged) from miniPFIsolation()"),
        miniRelIso_all = Var("(miniPFIsolation().chargedHadronIso()"+\
                             "+miniPFIsolation().neutralHadronIso()"+\
                             "+miniPFIsolation().photonIso())/pt",float,
                             doc="mini PF rel. iso. (total) from miniPFIsolation()"),
        miniPFRelIso_chg = Var("userFloat('miniIsoChg')/pt",float,
                               doc="mini PF rel. iso., charged component"),
        miniPFRelIso_all = Var("userFloat('miniIsoAll')/pt",float,
                               doc="mini PF rel. iso., total (with scaled rho*EA PU corrections)"),
        # Misc
        convVeto = Var("passConversionVeto()",bool,doc="pass conversion veto"),
        lostHits = Var("gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS')","uint8",doc="number of missing inner hits"),
        # Cluster-related
        energyErr = Var("p4Error('P4_COMBINATION')",float,doc="energy error of the cluster-track combination",precision=6),
        deltaEtaSC = Var("superCluster().eta()-eta()",float,doc="delta eta (SC,ele) with sign",precision=10),
        r9 = Var("full5x5_r9()",float,doc="R9 of the SC, calculated with full 5x5 region",precision=10),
        sieie = Var("full5x5_sigmaIetaIeta()",float,doc="sigma_IetaIeta of the SC, calculated with full 5x5 region",precision=10),
        eInvMinusPInv = Var("(1-eSuperClusterOverP())/ecalEnergy()",float,doc="1/E_SC - 1/p_trk",precision=10),
        scEtOverPt = Var("(superCluster().energy()/(pt*cosh(superCluster().eta())))-1",float,doc="(SC energy)/pt-1",precision=8),
        hoe = Var("hadronicOverEm()",float,doc="H over E",precision=8),
        # Displacement
        dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV, in cm",precision=10),
        # Cross-referencing
        #jetIdx
        #photonIdx
    ),
)

################################################################################
# electronTable (MC)
################################################################################

lowPtElectronsMCMatchForTable = cms.EDProducer(
    "MCMatcher",                                     # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = lowPtElectronTable.src,            # final reco collection
    matched     = cms.InputTag("finalGenParticles"), # final mc-truth particle collection
    mcPdgId     = cms.vint32(11),           # one or more PDG ID (11 = ele); absolute values (see below)
    checkCharge = cms.bool(False),          # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),            # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.3),          # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),          # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True), # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
)

lowPtElectronMCTable = cms.EDProducer(
    "CandMCMatchTableProducer",
    src     = lowPtElectronTable.src,
    mcMap   = cms.InputTag("lowPtElectronsMCMatchForTable"),
    objName = lowPtElectronTable.name,
    objType = lowPtElectronTable.name,
    branchName = cms.string("genPart"),
    docString = cms.string("MC matching to status==1 electrons"),
)

################################################################################
# Sequences
################################################################################

lowPtElectronSequence = cms.Sequence(updatedLowPtElectrons
                                     #+lowPtPATElectronID
                                     +isoForLowPtEle
                                     +updatedLowPtElectronsWithUserData
                                     +finalLowPtElectrons)
lowPtElectronTables = cms.Sequence(lowPtElectronTable)
lowPtElectronMC = cms.Sequence(lowPtElectronsMCMatchForTable
                               +lowPtElectronMCTable)
