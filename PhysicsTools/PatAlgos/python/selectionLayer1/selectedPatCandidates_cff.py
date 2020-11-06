import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.lowPtElectronSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.ootPhotonSelector_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
#from PhysicsTools.PatAlgos.producersLayer1.hemisphereProducer_cfi import *

# One module to count objects
selectedPatCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    logName = cms.untracked.string("selectedPatCanddiates|PATSummaryTables"),
    candidates = cms.VInputTag(
        cms.InputTag("selectedPatElectrons"),
        cms.InputTag("selectedPatLowPtElectrons"),
        cms.InputTag("selectedPatMuons"),
        cms.InputTag("selectedPatTaus"),
        cms.InputTag("selectedPatPhotons"),
        cms.InputTag("selectedPatOOTPhotons"),
        cms.InputTag("selectedPatJets"),
    )
)

selectedPatCandidatesTask = cms.Task(
    selectedPatElectrons,
    selectedPatLowPtElectrons,
    selectedPatMuons,
    selectedPatTaus,
    selectedPatPhotons,
    selectedPatOOTPhotons,
    selectedPatJets
)

selectedPatCandidates = cms.Sequence(selectedPatCandidateSummary, selectedPatCandidatesTask)

<<<<<<< HEAD
from Configuration.ProcessModifiers.pp_on_AA_cff import pp_on_AA
pp_on_AA.toReplaceWith(selectedPatCandidatesTask, selectedPatCandidatesTask.copyAndExclude([selectedPatOOTPhotons]))
pp_on_AA.toModify(selectedPatCandidateSummary.candidates, func = lambda list: list.remove(cms.InputTag("selectedPatOOTPhotons")) )
=======
from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
from Configuration.Eras.Modifier_pp_on_PbPb_run3_cff import pp_on_PbPb_run3
HI = (pp_on_AA_2018 | pp_on_PbPb_run3)
HI.toReplaceWith(selectedPatCandidatesTask,
                 selectedPatCandidatesTask.copyAndExclude([selectedPatOOTPhotons]))
HI.toModify(selectedPatCandidateSummary.candidates,
            func = lambda list: list.remove(cms.InputTag("selectedPatOOTPhotons")))

from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17
from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
mAOD = (run2_miniAOD_94XFall17 | run2_miniAOD_80XLegacy)
(HI | mAOD).toReplaceWith(selectedPatCandidatesTask,
                          selectedPatCandidatesTask.copyAndExclude([selectedPatLowPtElectrons]))
(HI | mAOD).toModify(selectedPatCandidateSummary.candidates,
                     func = lambda list: list.remove(cms.InputTag("selectedPatLowPtElectrons")) )
>>>>>>> remove low-pT electrons from PAT steps for run2_miniAOD_94XFall17 and run2_miniAOD_80XLegacy
