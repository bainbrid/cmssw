#ifndef RecoEgamma_ElectronIdentification_ElectronMVAEstimatorRun2_H
#define RecoEgamma_ElectronIdentification_ElectronMVAEstimatorRun2_H

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

class ElectronMVAEstimatorRun2 : public AnyMVAEstimatorRun2Base{

 public:

  // Constructor and destructor
  ElectronMVAEstimatorRun2(const edm::ParameterSet& conf);
  ~ElectronMVAEstimatorRun2() override;
  // For use with FWLite/Python
  ElectronMVAEstimatorRun2(const std::string &mvaTag,
                           const std::string &mvaName,
                           bool withIso,
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
  std::vector<float> fillMVAVariables(const edm::Ptr<reco::Candidate>& particle, const edm::Event&, const int iCategory) const;

  template<class EventType>
  std::vector<float> fillMVAVariables(const edm::Ptr<reco::GsfElectron>& eleRecoPtr, const EventType& iEvent, const int iCategory) const;

  int findCategory( const edm::Ptr<reco::Candidate>& particle) const override;
  int findCategory( const edm::Ptr<reco::GsfElectron>& particle) const ;

  // Call this function once after the constructor to declare
  // the needed event content pieces to the framework
  void setConsumes(edm::ConsumesCollector&&) const final;

 protected:

  // This class stores all the information that is needed about a variable
  class Variable {
    public:
      Variable(std::string &name, std::string &formula, bool hasLowerClip, bool hasUpperClip, float lowerClipValue, float upperClipValue, bool fromVariableHelper) {
          name_ =  name;
          formula_ = formula;
          hasLowerClip_ = hasLowerClip;
          hasUpperClip_ = hasUpperClip;
          lowerClipValue_ = lowerClipValue;
          upperClipValue_ = upperClipValue;
          fromVariableHelper_ = fromVariableHelper;

          if (!fromVariableHelper) {
              function_ = std::make_shared<StringObjectFunction<reco::GsfElectron>>(formula);
          }
      };

      std::string getName() const {
          return name_;
      };

      template<class EventType>
      float getValue(const edm::Ptr<reco::GsfElectron>& eleRecoPtr, const EventType& iEvent) const {
          float value;
          if (fromVariableHelper_) {
              edm::Handle<edm::ValueMap<float>> vMap;
              iEvent.getByLabel(edm::InputTag(formula_), vMap);
              value = (*vMap)[eleRecoPtr];
          } else {
              value = (*function_)(*eleRecoPtr);
          }
          if (hasLowerClip_ && value < lowerClipValue_) {
              value = lowerClipValue_;
          }
          if (hasUpperClip_ && value > upperClipValue_) {
              value = upperClipValue_;
          }
          return value;
      };

    private:
      std::string name_;
      std::string formula_;
      bool hasLowerClip_;
      bool hasUpperClip_;
      float lowerClipValue_;
      float upperClipValue_;
      bool fromVariableHelper_;
      std::shared_ptr<StringObjectFunction<reco::GsfElectron>> function_;
  };

  // Define here the number and the meaning of the categories
  // for this specific MVA
  //const int nCategories_ = 6;
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

  // The number of categories and number of variables per category
  int nCategories_;
  std::vector<int> nVariables_;

  // Data members
  std::vector< std::unique_ptr<const GBRForest> > gbrForests_;

  const std::string methodName_;

  double ptSplit_;   // we have high and low pt categories
  double ebSplit_;    // barrel is split into two regions
  double ebeeSplit_; // division between barrel and endcap

  // There might be different variables for each category, so the variables
  // names vector is itself a vector of length nCategories
  std::vector<std::vector<Variable>> variables_;

  bool withIso_;

  bool debug_;

  const std::string variableDefinitionFileName_;
};

#endif
