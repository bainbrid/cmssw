#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
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
  const edm::EDGetTokenT<edm::View<reco::Conversion>> convT_;
  edm::Handle<edm::View<reco::Conversion>> convH_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotT_;
  edm::Handle<reco::BeamSpot> beamSpotH_;
  const edm::EDGetTokenT<reco::VertexCollection> verticesT_;
  edm::Handle<reco::VertexCollection> verticesH_;
  bool extra_;
};

////////////////////////////////////////////////////////////////////////////////
//
LowPtElectronModifier::LowPtElectronModifier(const edm::ParameterSet& conf, edm::ConsumesCollector& cc)
    : ModifyObjectValueBase(conf),
      convT_(cc.consumes<edm::View<reco::Conversion>>(conf.getParameter<edm::InputTag>("conversions"))),
      convH_(),
      beamSpotT_(cc.consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpot"))),
      beamSpotH_(),
      verticesT_(cc.consumes<std::vector<reco::Vertex>>(conf.getParameter<edm::InputTag>("vertices"))),
      verticesH_(),
      extra_(conf.getParameter<bool>("addExtraUserVars")) {
  ;
}

////////////////////////////////////////////////////////////////////////////////
//
void LowPtElectronModifier::setEvent(const edm::Event& iEvent) {
  iEvent.getByToken(convT_, convH_);
  if (!convH_.isValid()) {
    throw cms::Exception("InvalidHandle", "Problem accessing reco::Conversions collection");
  }
  iEvent.getByToken(beamSpotT_, beamSpotH_);
  if (!beamSpotH_.isValid()) {
    throw cms::Exception("InvalidHandle", "Problem accessing reco::BeamSpot collection");
  }
  iEvent.getByToken(verticesT_, verticesH_);
  if (!verticesH_.isValid()) {
    throw cms::Exception("InvalidHandle", "Problem accessing reco::Vertex collection");
  }
}

////////////////////////////////////////////////////////////////////////////////
//
void LowPtElectronModifier::setEventContent(const edm::EventSetup& iSetup) {}

////////////////////////////////////////////////////////////////////////////////
//
void LowPtElectronModifier::modifyObject(pat::Electron& ele) const {
  // Embed Conversion info
  LowPtConversion conv;
  LowPtConversion::match(beamSpotH_, convH_, ele, conv);
  conv.addUserVars(ele);
  if (extra_) {
    conv.addExtraUserVars(ele);
  }
  // Set impact parameters
  if (!verticesH_->empty()) {
    const reco::Vertex& pv = verticesH_->front();
    ele.setDB(ele.gsfTrack()->dxy(pv.position()),
              ele.gsfTrack()->dxyError(pv.position(), pv.covariance()),
              pat::Electron::PV2D);  // PV2D
    ele.setDB(ele.gsfTrack()->dz(pv.position()),
              std::hypot(ele.gsfTrack()->dzError(), pv.zError()),
              pat::Electron::PVDZ);  // PVDZ
  }
  ele.setDB(ele.gsfTrack()->dxy(*beamSpotH_), ele.gsfTrack()->dxyError(*beamSpotH_),
            pat::Electron::BS2D);  // BS2D
}

////////////////////////////////////////////////////////////////////////////////
//
DEFINE_EDM_PLUGIN(ModifyObjectValueFactory, LowPtElectronModifier, "LowPtElectronModifier");
