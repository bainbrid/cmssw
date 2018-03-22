import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools import *

# Documentation of the MVA
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2
# https://rembserj.web.cern.ch/rembserj/notes/Electron_MVA_ID_2017_documentation

#
# In this file we define the locations of the MVA weights, cuts on the MVA values
# for specific working points, and configure those cuts in VID
#

# The tag is an extra string attached to the names of the products
# such as ValueMaps that needs to distinguish cases when the same MVA estimator
# class is used with different tuning/weights
mvaTag = "Fall17NoIsoV2"

# There are 6 categories in this MVA. They have to be configured in this strict order
# (cuts and weight files order):
#   0   EB1 (eta<0.8)  pt 5-10 GeV     |   pt < ptSplit && |eta| < ebSplit
#   1   EB2 (eta>=0.8) pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   2   EE             pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebeeSplit
#   3   EB1 (eta<0.8)  pt 10-inf GeV   |   pt >= ptSplit && |eta| < ebSplit
#   4   EB2 (eta>=0.8) pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   5   EE             pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebeeSplit


mvaFall17WeightFiles_V2 = cms.vstring(
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_noIso_V2_EB1_5.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_noIso_V2_EB2_5.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_noIso_V2_EE_5.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_noIso_V2_EB1_10.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_noIso_V2_EB2_10.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_noIso_V2_EE_10.weights.xml"
    )

## The working point for this MVA that is expected to have about 90% signal
# WP tuned to give about 90 and 80% signal efficiecny for electrons from Drell-Yan with pT > 25 GeV
# The working point for the low pt categories is just taken over from the high pt
idName90 = "mvaEleID-Fall17-noIso-V2-wp90"
MVA_WP90 = EleMVARaw_WP(
    idName = idName90, mvaTag = mvaTag,
    cutCategory0 = "3.3047416314137683 - exp(-pt /  4.611633062427176) *  7.424015786705544 ", # EB1 low pt
    cutCategory1 = "2.0786464662928745 - exp(-pt / 2.0612007367173075) * 13.352125309724341 ", # EB2 low pt
    cutCategory2 = "2.1255209997763536 - exp(-pt / 1.9747852821628373) * 17.706228113364627 ", # EE low pt
    cutCategory3 = " 5.918311416695114 - exp(-pt / 13.582265968662094) *  9.165579811721944 ", # EB1
    cutCategory4 = "4.9030319389144905 - exp(-pt / 13.332700807017755) *  9.033779959079803 ", # EB2
    cutCategory5 = " 4.120553396223688 - exp(-pt / 13.269186969731724) *  8.825818561343114 ", # EE
    )

idName80 = "mvaEleID-Fall17-noIso-V2-wp80"
MVA_WP80 = EleMVARaw_WP(
    idName = idName80, mvaTag = mvaTag,
    cutCategory0 = "3.6453530020414115 - exp(-pt / 3.9082833638220147) * 7.8309551284673375 ", # EB1 low pt
    cutCategory1 = "2.9554004255474244 - exp(-pt /  2.045291042897663) * 12.605893748230844 ", # EB2 low pt
    cutCategory2 = " 3.253550061804974 - exp(-pt / 1.8519992293714957) * 19.307915265477774 ", # EE low pt
    cutCategory3 = " 7.069090518356044 - exp(-pt / 16.259536070266417) *  8.122166134990497 ", # EB1
    cutCategory4 = " 6.059154751478629 - exp(-pt / 15.414188530991165) *  7.746729971800896 ", # EB2
    cutCategory5 = " 5.362818717135066 - exp(-pt / 15.621000722732715) *  7.342612013384582 ", # EE
)

### WP tuned for HZZ analysis with very high efficiency (about 98%)
# The working points were found by requiring the same signal efficiencies in
# each category as for the Spring 16 HZZ ID
# (see RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Spring16_HZZ_V2_cff.py)
idNamewpLoose = "mvaEleID-Fall17-noIso-V2-wpLoose"
MVA_WPLoose = EleMVARaw_WP(
    idName = idNamewpLoose, mvaTag = mvaTag,
    cutCategory0 = "-0.5162221154815269",  # EB1 low pt
    cutCategory1 = "-1.0232856117484725",  # EB2 low pt
    cutCategory2 = "-0.4609807572335405",  # EE low pt
    cutCategory3 = "-0.2697591050874118",  # EB1
    cutCategory4 = "-0.3714291214199583",  # EB2
    cutCategory5 = "-0.1198743839863436"   # EE
    )

#
# Finally, set up VID configuration for all cuts
#

# Create the PSet that will be fed to the MVA value map producer
mvaEleID_Fall17_noIso_V2_producer_config = cms.PSet(
    mvaName             = cms.string(mvaClassName),
    mvaTag              = cms.string(mvaTag),
    # Category parameters
    nCategories         = cms.int32(6),
    categoryCuts        = EleMVA_6CategoriesCuts,
    # Weight files and variable definitions
    weightFileNames     = mvaFall17WeightFiles_V2,
    variableDefinition  = cms.string(mvaVariablesFile)
    )
# Create the VPset's for VID cuts
mvaEleID_Fall17_V2_wpLoose = configureVIDMVAEleID( MVA_WPLoose )
mvaEleID_Fall17_V2_wp90 = configureVIDMVAEleID( MVA_WP90 )
mvaEleID_Fall17_V2_wp80 = configureVIDMVAEleID( MVA_WP80 )

mvaEleID_Fall17_V2_wpLoose.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wp90.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wp80.isPOGApproved = cms.untracked.bool(True)
