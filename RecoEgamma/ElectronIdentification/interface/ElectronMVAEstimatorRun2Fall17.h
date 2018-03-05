#ifndef RecoEgamma_ElectronIdentification_ElectronMVAEstimatorRun2Fall17_H
#define RecoEgamma_ElectronIdentification_ElectronMVAEstimatorRun2Fall17_H

#include <vector>
#include <string>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "RecoEgamma/EgammaTools/interface/AnyMVAEstimatorRun2Base.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/GBRForestTools.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class ElectronMVAEstimatorRun2Fall17 : public AnyMVAEstimatorRun2Base{

 public:

  // Constructor and destructor
  ElectronMVAEstimatorRun2Fall17(const edm::ParameterSet& conf, bool withIso);
  ~ElectronMVAEstimatorRun2Fall17() override;
  // For use with FWLite/Python
  ElectronMVAEstimatorRun2Fall17(const std::string &mvaTag,
                                 const std::string &mvaName,
                                 bool withIso,
                                 const std::string &conversionsTag = "reducedEgamma:reducedConversions",
                                 const std::string &beamspotTag = "offlineBeamSpot",
                                 const double ptSplit = 10., // The category split parameters are taken over from the python configuration file
                                 const double ebSplit = 0.800,
                                 const double ebeeSplit = 1.479,
                                 const bool debug = false);

  void init(const std::vector<std::string> &weightFileNames);
  void setClips(const std::vector<double> &clipsLowerValues, const std::vector<double> &clipsUpperValues);

  // Calculation of the MVA value (VID accessor)
  float mvaValue( const edm::Ptr<reco::Candidate>& particle, const edm::Event&) const override;
  // Calculation of the MVA value (fwlite-compatible accessor)
  float mvaValue( const edm::Ptr<reco::GsfElectron>& particle, const edm::EventBase & iEvent) const ;
  // Calculation of the MVA value (bare version)
  float mvaValue( const int iCategory, const std::vector<float> & vars) const ;

  // Utility functions
  int getNCategories() const override { return nCategories_; }
  const std::string& getName() const final { return name_; }
  const std::string& getTag() const final { return tag_; }

  // Functions that should work on both pat and reco electrons
  // (use the fact that pat::Electron inherits from reco::GsfElectron)
  std::vector<float> fillMVAVariables(const edm::Ptr<reco::Candidate>& particle, const edm::Event&) const override;

  template<class EventType>
  std::vector<float> fillMVAVariables(const edm::Ptr<reco::GsfElectron>& eleRecoPtr, const EventType& iEvent) const;

  int findCategory( const edm::Ptr<reco::Candidate>& particle) const override;
  int findCategory( const edm::Ptr<reco::GsfElectron>& particle) const ;
  // The function below ensures that the variables passed to MVA are
  // within reasonable bounds
  void constrainMVAVariables(std::vector<float>&) const;

  // Call this function once after the constructor to declare
  // the needed event content pieces to the framework
  void setConsumes(edm::ConsumesCollector&&) const final;

 protected:

  // This is a stuct to store an instruction of clipping one of the input variables with a lower or an upper limit
  struct Clip {
      unsigned int varIdx;
      bool upper;
      float value;
  };

  // Define here the number and the meaning of the categories
  // for this specific MVA
  const int nCategories_ = 6;
  const int nVar_ = 22;
  enum MVACategories_ {
    UNDEFINED = -1,
    CAT_EB1_PTLow = 0,
    CAT_EB2_PTLow = 1,
    CAT_EE_PTLow  = 2,
    CAT_EB1_PTHig = 3,
    CAT_EB2_PTHig = 4,
    CAT_EE_PTHig  = 5
  };

  // MVA tag. This is an additional string variable to distinguish
  // instances of the estimator of this class configured with different
  // weight files.
  const std::string tag_;

  // MVA name. This is a unique name for this MVA implementation.
  // It will be used as part of ValueMap names.
  // For simplicity, keep it set to the class name.

  const std::string name_;

  // Data members
  std::vector< std::unique_ptr<const GBRForest> > gbrForests_;

  const std::string methodName_;

  //edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  //edm::EDGetTokenT<reco::ConversionCollection> conversionsTokenAOD_;
  //edm::EDGetTokenT<reco::ConversionCollection> conversionsTokenMiniAOD_;
  //edm::EDGetTokenT<double> rhoToken_;

  double ptSplit_;   // we have high and low pt categories
  double ebSplit_;    // barrel is split into two regions
  double ebeeSplit_; // division between barrel and endcap

  std::vector<std::string> varNames_;

  // To store the variable clipping operations
  std::vector<Clip> clipsLower_;
  std::vector<Clip> clipsUpper_;

  bool withIso_;

  bool debug_;

  // The functions to obtain the input variables

  /*
   * Defining variables for the MVA evaluation.
   * CAREFUL: It is critical that all the variables that are packed into “vars” are
   * exactly in the order they are found in the weight files
   */
  const std::string gsfEleFuncStrings_ [21] = {
      // Pure ECAL -> shower shapes
      "full5x5_sigmaIetaIeta",
      "full5x5_sigmaIphiIphi",
      "1. - full5x5_e1x5 / full5x5_e5x5",
      "full5x5_r9",
      "superCluster.etaWidth",
      "superCluster.phiWidth",
      "full5x5_hcalOverEcal", //hadronicOverEm
      //Pure tracking variables
      "gsfTrack.normalizedChi2",
      // Energy matching
      "fbrem",
      "gsfTrack.hitPattern.trackerLayersWithMeasurement",
      "gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS')",
      // Energy matching
      "eSuperClusterOverP",
      "eEleClusterOverPout",
      "1.0 / ecalEnergy - 1.0 / trackMomentumAtVtx.R",
      // Geometrical matchings
      "deltaEtaSuperClusterTrackAtVtx",
      "deltaPhiSuperClusterTrackAtVtx",
      "deltaEtaSeedClusterTrackAtCalo",
      // Isolation variables
      "pfIsolationVariables.sumChargedHadronPt", //chargedHadronIso
      "pfIsolationVariables.sumNeutralHadronEt", //neutralHadronIso
      "pfIsolationVariables.sumPhotonEt", //photonIso
      // Endcap only
      "superCluster.preshowerEnergy / superCluster.rawEnergy"
  };
  std::vector<StringObjectFunction<reco::GsfElectron>> gsfEleFunctions_;
};

#endif
