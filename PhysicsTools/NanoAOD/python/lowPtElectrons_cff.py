import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *

################################################################################
# slimmedLowPtElectrons_WithUserData
################################################################################

isoForLowPtEle = cms.EDProducer(
    "EleIsoValueMapProducer",
    src = cms.InputTag("slimmedLowPtElectrons"),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    rho_PFIso = cms.InputTag("fixedGridRhoFastjetAll"),
    EAFile_MiniIso = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
    EAFile_PFIso = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
)

slimmedLowPtElectronsWithUserData = cms.EDProducer(
    "PATElectronUserDataEmbedder",
    src = cms.InputTag("slimmedLowPtElectrons"),
    userFloats = cms.PSet(
        miniIsoChg = cms.InputTag("isoForLowPtEle:miniIsoChg"),
        miniIsoAll = cms.InputTag("isoForLowPtEle:miniIsoAll"),
    ),
    userIntFromBools = cms.PSet(),
    userInts = cms.PSet(),
    userCands = cms.PSet(),
)

################################################################################
# finalElectrons
################################################################################

finalLowPtElectrons = cms.EDFilter(
    "PATElectronRefSelector",
    src = cms.InputTag("slimmedLowPtElectronsWithUserData"),
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
        miniPFRelIso_chg = Var("userFloat('miniIsoChg')/pt",float,doc="mini PF rel. iso., charged component"),
        miniPFRelIso_all = Var("userFloat('miniIsoAll')/pt",float,doc="mini PF rel. iso., total (with scaled rho*EA PU corrections)"),
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
# electronTable END
################################################################################

lowPtElectronSequence = cms.Sequence(isoForLowPtEle + slimmedLowPtElectronsWithUserData + finalLowPtElectrons)
lowPtElectronTables = cms.Sequence(lowPtElectronTable)
