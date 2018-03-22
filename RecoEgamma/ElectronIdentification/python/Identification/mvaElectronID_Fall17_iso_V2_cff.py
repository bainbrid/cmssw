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
mvaTag = "Fall17IsoV2"

# There are 6 categories in this MVA. They have to be configured in this strict order
# (cuts and weight files order):
#   0   EB1 (eta<0.8)  pt 5-10 GeV     |   pt < ptSplit && |eta| < ebSplit
#   1   EB2 (eta>=0.8) pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   2   EE             pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebeeSplit
#   3   EB1 (eta<0.8)  pt 10-inf GeV   |   pt >= ptSplit && |eta| < ebSplit
#   4   EB2 (eta>=0.8) pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   5   EE             pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebeeSplit


mvaFall17WeightFiles_V2 = cms.vstring(
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB1_5.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB2_5.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EE_5.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB1_10.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB2_10.weights.xml",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EE_10.weights.xml"
    )

## The working point for this MVA that is expected to have about 90% signal
# WP tuned to give about 90 and 80% signal efficiecny for electrons from Drell-Yan with pT > 25 GeV
# The working point for the low pt categories is just taken over from the high pt
idName90 = "mvaEleID-Fall17-iso-V2-wp90"
MVA_WP90 = EleMVARaw_WP(
    idName = idName90, mvaTag = mvaTag,
    cutCategory0 = "3.3300182400168525  - exp(-pt /  3.998316465442375) *  8.099043200084127", # EB1 low pt
    cutCategory1 = " 2.260175667879445  - exp(-pt / 1.9189641506694057) *   13.9102227144152", # EB2 low pt
    cutCategory2 = " 2.211562592402872  - exp(-pt / 1.8588312432333167) * 21.073158261818474", # EE low pt
    cutCategory3 = "6.1646928728924255  - exp(-pt / 13.215485781314957) *  8.723072659853388", # EB1
    cutCategory4 = " 5.220773396152184  - exp(-pt /  13.76692857786133) *  8.205427471565141", # EB2
    cutCategory5 = " 4.384204383058256  - exp(-pt / 14.690856427853634) *  8.198432664405424", # EE
)

idName80 = "mvaEleID-Fall17-iso-V2-wp80"
MVA_WP80 = EleMVARaw_WP(
    idName = idName80, mvaTag = mvaTag,
    cutCategory0 = "  3.879373677541345 - exp(-pt /   3.56640268963677) *  8.462849664367083", # EB1 low pt
    cutCategory1 = " 3.2060589508492483 - exp(-pt / 1.7930465020050448) * 15.846826583233312", # EB2 low pt
    cutCategory2 = "  3.368655131557081 - exp(-pt /  1.722605342040908) * 24.143220505560134", # EE low pt
    cutCategory3 = " 7.2998806130831975 - exp(-pt / 15.249064444355554) *  7.635961722826079", # EB1
    cutCategory4 = "  6.354710470440623 - exp(-pt /  15.44798606129189) * 7.0022714937821595", # EB2
    cutCategory5 = "  5.640966612632747 - exp(-pt / 17.344950075131518) *  6.860659136176625", # EE
)

### WP tuned for HZZ analysis with very high efficiency (about 98%)
# The working points were found by requiring the same signal efficiencies in
# each category as for the Spring 16 HZZ ID
idNamewpLoose = "mvaEleID-Fall17-iso-V2-wpLoose"
MVA_WPLoose = EleMVARaw_WP(
    idName = idNamewpLoose, mvaTag = mvaTag,
    cutCategory0 = "-0.5036707043540236",    # EB1 low pt
    cutCategory1 = "-0.9419498655189662",    # EB2 low pt
    cutCategory2 = "-0.420452615888308",     # EE low pt
    cutCategory3 = "-0.07112322732699157",   # EB1
    cutCategory4 = "-0.11684180190115502",   # EB2
    cutCategory5 = "-0.030233133764291563"   # EE
    )

### WP tuned for HZZ analysis with very high efficiency (about 98%)
# The working points were found by requiring the same signal efficiencies in
# each category as for the Spring 16 HZZ ID plus isolation cut
idNamewpHZZ = "mvaEleID-Fall17-iso-V2-wpHZZ"
MVA_WPHZZ = EleMVARaw_WP(
    idName = idNamewpHZZ, mvaTag = mvaTag,
    cutCategory0 = "1.2341054581543813",   # EB1 low pt
    cutCategory1 = "1.1512280403260653",   # EB2 low pt
    cutCategory2 = "1.3404835451646275",   # EE low pt
    cutCategory3 = "2.435790605727733",    # EB1
    cutCategory4 = "1.9587613557218568",   # EB2
    cutCategory5 = "1.0805127088478825"    # EE
    )

#
# Finally, set up VID configuration for all cuts
#

# Create the PSet that will be fed to the MVA value map producer
mvaEleID_Fall17_iso_V2_producer_config = cms.PSet(
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
mvaEleID_Fall17_V2_wpHZZ = configureVIDMVAEleID( MVA_WPHZZ )
mvaEleID_Fall17_V2_wp90 = configureVIDMVAEleID( MVA_WP90)
mvaEleID_Fall17_V2_wp80 = configureVIDMVAEleID( MVA_WP80)

mvaEleID_Fall17_V2_wpHZZ.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wpLoose.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wp90.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wp80.isPOGApproved = cms.untracked.bool(True)
