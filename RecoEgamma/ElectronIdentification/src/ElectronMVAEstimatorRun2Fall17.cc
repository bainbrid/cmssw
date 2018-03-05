#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimatorRun2Fall17.h"

ElectronMVAEstimatorRun2Fall17::ElectronMVAEstimatorRun2Fall17(const edm::ParameterSet& conf, bool withIso):
  AnyMVAEstimatorRun2Base(conf),
  tag_(conf.getParameter<std::string>("mvaTag")),
  name_(conf.getParameter<std::string>("mvaName")),
  methodName_("BDTG method"),
  ptSplit_                (conf.getParameter<double>                  ("ptSplit")),
  ebSplit_                (conf.getParameter<double>                  ("ebSplit")),
  ebeeSplit_              (conf.getParameter<double>                  ("ebeeSplit")),
  varNames_               (conf.getParameter<std::vector<std::string>>("varNames"))
{

  const std::vector <std::string> weightFileNames
    = conf.getParameter<std::vector<std::string> >("weightFileNames");

  std::vector<double> clipsLowerValues = conf.getParameter<std::vector<double>>("clipLower");
  std::vector<double> clipsUpperValues = conf.getParameter<std::vector<double>>("clipUpper");

  // Initialize GBRForests and and clipping instructions
  init(weightFileNames);
  setClips(clipsLowerValues, clipsUpperValues);

  withIso_ = withIso;
  debug_ = conf.getUntrackedParameter<bool>("debug", false);

}

ElectronMVAEstimatorRun2Fall17::ElectronMVAEstimatorRun2Fall17(
        const std::string &mvaTag, const std::string &mvaName, bool withIso, const std::string &conversionsTag, const std::string &beamspotTag,
        const double ptSplit, const double ebSplit, const double ebeeSplit, const bool debug):
  AnyMVAEstimatorRun2Base( edm::ParameterSet() ),
  tag_                    (mvaTag),
  name_                   (mvaName),
  methodName_             ("BDTG method"),
  ptSplit_                (ptSplit),
  ebSplit_                (ebSplit),
  ebeeSplit_              (ebeeSplit)
{

  // Set if this is the ID with or without PF isolations
  withIso_ = withIso;

  // Set debug flag
  debug_ = debug;

}

void ElectronMVAEstimatorRun2Fall17::init(const std::vector<std::string> &weightFileNames) {

  // Initialize GBRForests
  if( (int)(weightFileNames.size()) != nCategories_ )
    throw cms::Exception("MVA config failure: ")
      << "wrong number of weightfiles" << std::endl;

  gbrForests_.clear();
  // Create a TMVA reader object for each category
  for(int i=0; i<nCategories_; i++){

    // Use unique_ptr so that all readers are properly cleaned up
    // when the vector clear() is called in the destructor

    edm::FileInPath weightFile( weightFileNames[i] );
    gbrForests_.push_back( GBRForestTools::createGBRForest( weightFile ) );

  }

  // Initialize the functions
  gsfEleFunctions_.clear();
  for (int i = 0; i < 21; ++i) {
      StringObjectFunction<reco::GsfElectron> f(gsfEleFuncStrings_[i]);
      gsfEleFunctions_.push_back(f);
  }

}

void ElectronMVAEstimatorRun2Fall17::setClips(const std::vector<double> &clipsLowerValues, const std::vector<double> &clipsUpperValues) {

  // Set up the variable clipping intructions
  unsigned int i = 0;
  for(auto const& value: clipsLowerValues) {
      if (edm::isFinite(value)) {
          Clip clip = {i, false, (float)value};
          clipsLower_.push_back(clip);
      }
      ++i;
  }

  i = 0;
  for(auto const& value: clipsUpperValues) {
      if (edm::isFinite(value)) {
          Clip clip = {i, true, (float)value};
          clipsUpper_.push_back(clip);
      }
      ++i;
  }

}

ElectronMVAEstimatorRun2Fall17::
~ElectronMVAEstimatorRun2Fall17(){
}

void ElectronMVAEstimatorRun2Fall17::setConsumes(edm::ConsumesCollector&& cc) const {

  // All tokens for event content needed by this MVA

  // Tags from the variable helper
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:convVtxFitProb"));
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:kfchi2"));
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:kfhits"));
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:rho"));
}

float ElectronMVAEstimatorRun2Fall17::
mvaValue( const edm::Ptr<reco::Candidate>& particle, const edm::Event& iEvent) const {
  // Try to cast the particle into a reco particle.
  // This should work for both reco and pat.
  const edm::Ptr<reco::GsfElectron> eleRecoPtr = ( edm::Ptr<reco::GsfElectron> )particle;
  if( eleRecoPtr.get() == nullptr ) {
    throw cms::Exception("MVA failure: ")
      << " given particle is expected to be reco::GsfElectron or pat::Electron," << std::endl
      << " but appears to be neither" << std::endl;
  }

  const int iCategory = findCategory( eleRecoPtr );
  const std::vector<float> vars = fillMVAVariables( particle, iEvent );
  return mvaValue(iCategory, vars);
}

float ElectronMVAEstimatorRun2Fall17::
mvaValue( const edm::Ptr<reco::GsfElectron>& particle, const edm::EventBase & iEvent) const {

  const int iCategory = findCategory( particle );
  const std::vector<float> vars = fillMVAVariables( particle, iEvent );
  return mvaValue(iCategory, vars);
}

float ElectronMVAEstimatorRun2Fall17::
mvaValue( const int iCategory, const std::vector<float> & vars) const  {
  const float result = gbrForests_.at(iCategory)->GetClassifier(vars.data());

  if(debug_) {
    std::cout << " *** Inside the class methodName_ " << methodName_ << std::endl;
    std::cout << " bin "                      << iCategory << std::endl
              << " see "                      << vars[0] << std::endl
              << " spp "                      << vars[1] << std::endl
              << " circularity "              << vars[2] << std::endl
              << " r9 "                       << vars[3] << std::endl
              << " etawidth "                 << vars[4] << std::endl
              << " phiwidth "                 << vars[5] << std::endl
              << " hoe "                      << vars[6] << std::endl
              << " kfhits "                   << vars[7] << std::endl
              << " kfchi2 "                   << vars[8] << std::endl
              << " gsfchi2 "                  << vars[9] << std::endl
              << " fbrem "                    << vars[10] << std::endl
              << " gsfhits "                  << vars[11] << std::endl
              << " expectedMissingInnerHits " << vars[12] << std::endl
              << " convVtxFitProbability "    << vars[13] << std::endl
              << " eop "                      << vars[14] << std::endl
              << " eleeopout "                << vars[15] << std::endl
              << " oneOverEminusOneOverP "    << vars[16] << std::endl
              << " deta "                     << vars[17] << std::endl
              << " dphi "                     << vars[18] << std::endl
              << " detacalo "                 << vars[19] << std::endl;
    if (withIso_) {
      std::cout << " ele_pfPhotonIso "          << vars[20] << std::endl
                << " ele_pfChargedHadIso "      << vars[21] << std::endl
                << " ele_pfNeutralHadIso "      << vars[22] << std::endl
                << " rho "                      << vars[23] << std::endl
                << " preShowerOverRaw "         << vars[24] << std::endl;
    }
    else {
      std::cout << " rho "                      << vars[20] << std::endl
                << " preShowerOverRaw "         << vars[21] << std::endl;
    }
    std::cout << " ### MVA " << result << std::endl << std::endl;
  }

  return result;
}

int ElectronMVAEstimatorRun2Fall17::findCategory( const edm::Ptr<reco::Candidate>& particle) const {

  // Try to cast the particle into a reco particle.
  // This should work for both reco and pat.
  const edm::Ptr<reco::GsfElectron> eleRecoPtr = ( edm::Ptr<reco::GsfElectron> )particle;
  if( eleRecoPtr.get() == nullptr ) {
    throw cms::Exception("MVA failure: ")
      << " given particle is expected to be reco::GsfElectron or pat::Electron," << std::endl
      << " but appears to be neither" << std::endl;
  }
  return findCategory(eleRecoPtr);
}

int ElectronMVAEstimatorRun2Fall17::findCategory( const edm::Ptr<reco::GsfElectron>& eleRecoPtr ) const {
  float pt = eleRecoPtr->pt();
  float eta = eleRecoPtr->superCluster()->eta();
  float absEta = std::abs(eta);

  //
  // Determine the category
  //
  int  iCategory = UNDEFINED;

  if (pt < ptSplit_ && absEta < ebSplit_) {
    iCategory = CAT_EB1_PTLow;
  }

  else if (pt < ptSplit_ && absEta >= ebSplit_ && absEta < ebeeSplit_) {
    iCategory = CAT_EB2_PTLow;
  }

  else if (pt < ptSplit_ && absEta >= ebeeSplit_) {
    iCategory = CAT_EE_PTLow;
  }

  else if (pt >= ptSplit_ && absEta < ebSplit_) {
    iCategory = CAT_EB1_PTHig;
  }

  else if (pt >= ptSplit_ && absEta >= ebSplit_ && absEta < ebeeSplit_) {
    iCategory = CAT_EB2_PTHig;
  }

  else if (pt >= ptSplit_ && absEta >= ebeeSplit_) {
    iCategory = CAT_EE_PTHig;
  }

  return iCategory;
}

// A function that should work on both pat and reco objects
std::vector<float> ElectronMVAEstimatorRun2Fall17::
fillMVAVariables(const edm::Ptr<reco::Candidate>& particle, const edm::Event& iEvent) const {

  // Try to cast the particle into a reco particle.
  // This should work for both reco and pat.
  const edm::Ptr<reco::GsfElectron> eleRecoPtr = ( edm::Ptr<reco::GsfElectron> )particle;
  if( eleRecoPtr.get() == nullptr ) {
    throw cms::Exception("MVA failure: ")
      << " given particle is expected to be reco::GsfElectron or pat::Electron," << std::endl
      << " but appears to be neither" << std::endl;
  }

  return fillMVAVariables(eleRecoPtr, iEvent);
}

// A function that should work on both pat and reco objects
template<class EventType>
std::vector<float> ElectronMVAEstimatorRun2Fall17::
fillMVAVariables(const edm::Ptr<reco::GsfElectron>& eleRecoPtr, const EventType& iEvent) const {

  //
  // Declare all value maps corresponding to the products we defined earlier
  //

  // Get the variables from the helper class
  edm::Handle<edm::ValueMap<float>> kfhits;
  edm::Handle<edm::ValueMap<float>> kfchi2;
  edm::Handle<edm::ValueMap<float>> convVtxFitProb;
  edm::Handle<edm::ValueMap<float>> rho;

  iEvent.getByLabel(edm::InputTag("electronMVAVariableHelper:kfhits"), kfhits);
  iEvent.getByLabel(edm::InputTag("electronMVAVariableHelper:kfchi2"), kfchi2);
  iEvent.getByLabel(edm::InputTag("electronMVAVariableHelper:convVtxFitProb"), convVtxFitProb);
  iEvent.getByLabel(edm::InputTag("electronMVAVariableHelper:rho"), rho);

  if(withIso_)
  {

    std::vector<float> vars = packMVAVariables(
                                  (float)gsfEleFunctions_[0](*eleRecoPtr),
                                  (float)gsfEleFunctions_[1](*eleRecoPtr),
                                  (float)gsfEleFunctions_[2](*eleRecoPtr),
                                  (float)gsfEleFunctions_[3](*eleRecoPtr),
                                  (float)gsfEleFunctions_[4](*eleRecoPtr),
                                  (float)gsfEleFunctions_[5](*eleRecoPtr),
                                  (float)gsfEleFunctions_[6](*eleRecoPtr),
                                  //Pure tracking variables
                                  (float)(*kfhits)[eleRecoPtr],
                                  (float)(*kfchi2)[eleRecoPtr],
                                  (float)gsfEleFunctions_[7](*eleRecoPtr),
                                  (float)gsfEleFunctions_[8](*eleRecoPtr),
                                  (float)gsfEleFunctions_[9](*eleRecoPtr),
                                  (float)gsfEleFunctions_[10](*eleRecoPtr),
                                  (float)(*convVtxFitProb)[eleRecoPtr],
                                  (float)gsfEleFunctions_[11](*eleRecoPtr),
                                  (float)gsfEleFunctions_[12](*eleRecoPtr),
                                  (float)gsfEleFunctions_[13](*eleRecoPtr),
                                  (float)gsfEleFunctions_[14](*eleRecoPtr),
                                  (float)gsfEleFunctions_[15](*eleRecoPtr),
                                  (float)gsfEleFunctions_[16](*eleRecoPtr),
                                  (float)gsfEleFunctions_[17](*eleRecoPtr),
                                  (float)gsfEleFunctions_[18](*eleRecoPtr),
                                  (float)gsfEleFunctions_[19](*eleRecoPtr),
                                  // Pileup
                                  (float)(*rho)[eleRecoPtr],

                                  // Endcap only variables NOTE: we don't need
                                  // to check if we are actually in the endcap
                                  // or not, as it is the last variable in the
                                  // vector and it will be ignored by the
                                  // GBRForest for the barrel.
                                  //
                                  // The GBRForest classification just needs an
                                  // array with the input variables in the
                                  // right order, what comes after doesn't
                                  // matter.
                                  (float)gsfEleFunctions_[20](*eleRecoPtr)
                              );

    constrainMVAVariables(vars);

    return vars;
  }
  else
  {
    std::vector<float> vars = packMVAVariables(
                                  (float)gsfEleFunctions_[0](*eleRecoPtr),
                                  (float)gsfEleFunctions_[1](*eleRecoPtr),
                                  (float)gsfEleFunctions_[2](*eleRecoPtr),
                                  (float)gsfEleFunctions_[3](*eleRecoPtr),
                                  (float)gsfEleFunctions_[4](*eleRecoPtr),
                                  (float)gsfEleFunctions_[5](*eleRecoPtr),
                                  (float)gsfEleFunctions_[6](*eleRecoPtr),
                                  (float)(*kfhits)[eleRecoPtr],
                                  (float)(*kfchi2)[eleRecoPtr],
                                  (float)gsfEleFunctions_[7](*eleRecoPtr),
                                  (float)gsfEleFunctions_[8](*eleRecoPtr),
                                  (float)gsfEleFunctions_[9](*eleRecoPtr),
                                  (float)gsfEleFunctions_[10](*eleRecoPtr),
                                  (float)(*convVtxFitProb)[eleRecoPtr],
                                  (float)gsfEleFunctions_[11](*eleRecoPtr),
                                  (float)gsfEleFunctions_[12](*eleRecoPtr),
                                  (float)gsfEleFunctions_[13](*eleRecoPtr),
                                  (float)gsfEleFunctions_[14](*eleRecoPtr),
                                  (float)gsfEleFunctions_[15](*eleRecoPtr),
                                  (float)gsfEleFunctions_[16](*eleRecoPtr),
                                  (float)(*rho)[eleRecoPtr],
                                  (float)gsfEleFunctions_[20](*eleRecoPtr)
                              );

    constrainMVAVariables(vars);

    return vars;
  }

}

void ElectronMVAEstimatorRun2Fall17::constrainMVAVariables(std::vector<float>& vars) const {

  // Check that variables do not have crazy values

  for(auto const& clip: clipsLower_) {
      if ( vars[clip.varIdx] < clip.value  ) {
          vars[clip.varIdx] =   clip.value;
      }
  }

  for(auto const& clip: clipsUpper_) {
      if ( vars[clip.varIdx] > clip.value  ) {
          vars[clip.varIdx] =   clip.value;
      }
  }

}
