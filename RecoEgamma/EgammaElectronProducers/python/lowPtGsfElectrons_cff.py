import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronsPreRegression_cfi import lowPtGsfElectronsPreRegression

lowPtGsfElectrons = lowPtGsfElectronsPreRegression.clone()

################################################################################
# LowPtGsfElectronProducer above is run by default in RECO
# LowPtGsfElectronFinalizer below is scheduled for bParking and run2_miniAOD_devel

from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier106XUL
_lowPtRegressionModifier = regressionModifier106XUL.clone(
    modifierName = 'EGRegressionModifierV3',
    rhoTag = 'fixedGridRhoFastjetAll',
    eleRegs = dict(
        ecalOnlyMean = dict(
            lowEtHighEtBoundary = 20.,
            ebLowEtForestName = "lowPtElectron_eb_ecalOnly_05To50_mean",
            ebHighEtForestName = "lowPtElectron_eb_ecalOnly_05To50_mean",
            eeLowEtForestName = "lowPtElectron_ee_ecalOnly_05To50_mean",
            eeHighEtForestName = "lowPtElectron_ee_ecalOnly_05To50_mean",
            ),
        ecalOnlySigma = dict(
            lowEtHighEtBoundary = 20.,
            ebLowEtForestName = "lowPtElectron_eb_ecalOnly_05To50_sigma",
            ebHighEtForestName = "lowPtElectron_eb_ecalOnly_05To50_sigma",
            eeLowEtForestName = "lowPtElectron_ee_ecalOnly_05To50_sigma",
            eeHighEtForestName = "lowPtElectron_ee_ecalOnly_05To50_sigma",
            ),
        epComb = dict(
            ecalTrkRegressionConfig = dict(
                lowEtHighEtBoundary = 20.,
                ebLowEtForestName = "lowPtElectron_eb_ecalTrk_05To50_mean",
                ebHighEtForestName = "lowPtElectron_eb_ecalTrk_05To50_mean",
                eeLowEtForestName = "lowPtElectron_ee_ecalTrk_05To50_mean",
                eeHighEtForestName = "lowPtElectron_ee_ecalTrk_05To50_mean",
                ),
            ecalTrkRegressionUncertConfig = dict(
                lowEtHighEtBoundary = 20.,
                ebLowEtForestName = "lowPtElectron_eb_ecalTrk_05To50_sigma",
                ebHighEtForestName = "lowPtElectron_eb_ecalTrk_05To50_sigma",
                eeLowEtForestName = "lowPtElectron_ee_ecalTrk_05To50_sigma",
                eeHighEtForestName = "lowPtElectron_ee_ecalTrk_05To50_sigma",
                ),
        )
    ),
)

from RecoEgamma.EgammaElectronProducers.lowPtGsfElectronFinalizer_cfi import lowPtGsfElectronFinalizer
_lowPtGsfElectrons = lowPtGsfElectronFinalizer.clone(
    previousGsfElectronsTag = "lowPtGsfElectrons::@skipCurrentProcess",
    regressionConfig = _lowPtRegressionModifier,
)

from Configuration.Eras.Modifier_bParking_cff import bParking
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel
(bParking & ~run2_miniAOD_UL & ~run2_miniAOD_devel).toModify(
    _lowPtGsfElectrons,
    previousGsfElectronsTag = "lowPtGsfElectronsPreRegression",
)

(bParking | run2_miniAOD_UL | run2_miniAOD_devel).toReplaceWith(lowPtGsfElectrons,_lowPtGsfElectrons)
