import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *

################################################################################
# finalElectrons
################################################################################

finalLowPtElectrons = cms.EDFilter(
    "PATElectronRefSelector",
    src = cms.InputTag("slimmedLowPtElectrons"),
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
        CandVars,
        ID = Var("electronID('ID')",float,doc="ID, BDT (raw) score"),
        unbiased = Var("electronID('unbiased')",float,doc="ElectronSeed, pT- and dxy- agnostic BDT (raw) score"),
        ptbiased = Var("electronID('ptbiased')",float,doc="ElectronSeed, pT- and dxy- dependent BDT (raw) score"),
    ),
)

################################################################################
# electronTable END
################################################################################

lowPtElectronSequence = cms.Sequence(finalLowPtElectrons)
lowPtElectronTables = cms.Sequence(lowPtElectronTable)
