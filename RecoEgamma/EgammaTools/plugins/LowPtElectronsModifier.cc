#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoEgamma/EgammaTools/interface/LowPtConversion.h"

////////////////////////////////////////////////////////////////////////////////
//
class LowPtElectronModifier : public ModifyObjectValueBase {
public:
  LowPtElectronModifier(const edm::ParameterSet& conf, edm::ConsumesCollector&);
  ~LowPtElectronModifier() override{};

  void setEvent(const edm::Event&) final;
  void setEventContent(const edm::EventSetup&) final;

  void modifyObject(pat::Electron& ele) const final;

private:
  const edm::EDGetTokenT<edm::View<reco::Conversion> > convT_;
  edm::Handle<edm::View<reco::Conversion> > convH_;
  const edm::EDGetTokenT<reco::BeamSpot> bsT_;
  edm::Handle<reco::BeamSpot> bsH_;
  bool extra_;
};

////////////////////////////////////////////////////////////////////////////////
//
LowPtElectronModifier::LowPtElectronModifier(const edm::ParameterSet& conf, edm::ConsumesCollector& cc)
    : ModifyObjectValueBase(conf),
      convT_(cc.consumes<edm::View<reco::Conversion> >(conf.getParameter<edm::InputTag>("conversions"))),
      convH_(),
      bsT_(cc.consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpot"))),
      bsH_(),
      extra_(conf.getParameter<bool>("addExtraUserVars")) {
  ;
}

////////////////////////////////////////////////////////////////////////////////
//
void LowPtElectronModifier::setEvent(const edm::Event& iEvent) {
  iEvent.getByToken(convT_, convH_);
  if (!convH_.isValid()) {
    throw cms::Exception("InvalidHandle", "Problem accessing Conversions collection");
  }
  iEvent.getByToken(bsT_, bsH_);
  if (!convH_.isValid()) {
    throw cms::Exception("InvalidHandle", "Problem accessing BeamSpot collection");
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void LowPtElectronModifier::setEventContent(const edm::EventSetup& iSetup) {}

////////////////////////////////////////////////////////////////////////////////
//
void LowPtElectronModifier::modifyObject(pat::Electron& ele) const {
  LowPtConversion conv;
  LowPtConversion::match(bsH_, convH_, ele, conv);
  conv.addUserVars(ele);
  if (extra_) {
    conv.addExtraUserVars(ele);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
DEFINE_EDM_PLUGIN(ModifyObjectValueFactory, LowPtElectronModifier, "LowPtElectronModifier");
