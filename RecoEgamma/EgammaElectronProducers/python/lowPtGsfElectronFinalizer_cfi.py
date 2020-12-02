import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier106XUL

lowPtRegressionModifier = regressionModifier106XUL.clone(
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

lowPtGsfElectrons = cms.EDProducer("LowPtGsfElectronFinalizer",
                                   previousGsfElectronsTag = cms.InputTag("lowPtGsfElectronsPreRegression"),
                                   regressionConfig = lowPtRegressionModifier,
)
