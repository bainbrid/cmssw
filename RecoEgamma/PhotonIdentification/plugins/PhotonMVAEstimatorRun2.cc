#include "RecoEgamma/PhotonIdentification/plugins/PhotonMVAEstimatorRun2.h"

PhotonMVAEstimatorRun2::PhotonMVAEstimatorRun2(const edm::ParameterSet& conf):
  AnyMVAEstimatorRun2Base(conf),
  name_(conf.getParameter<std::string>("mvaName")),
  tag_(conf.getParameter<std::string>("mvaTag")),
  methodName_("BDTG method"),
  mvaVarMngr_(conf.getParameter<std::string>("variableDefinition")),
  debug_(conf.getUntrackedParameter<bool>("debug", false))
{
  //
  // Construct the MVA estimators
  //

  const std::vector <std::string> weightFileNames
    = conf.getParameter<std::vector<std::string> >("weightFileNames");

  // Initialize GBRForests
  if( (int)(weightFileNames.size()) != nCategories_ )
    throw cms::Exception("MVA config failure: ")
      << "wrong number of weightfiles" << std::endl;

  gbrForests_.clear();
  // Create a TMVA reader object for each category
  for(int i=0; i<nCategories_; i++){

    std::vector<std::string> variableNamesInCategory;
    std::vector<int> variablesInCategory;

    // Use unique_ptr so that all readers are properly cleaned up
    // when the vector clear() is called in the destructor

    gbrForests_.push_back( GBRForestTools::createGBRForest( weightFileNames[i], variableNamesInCategory ) );

    nVariables_.push_back(variableNamesInCategory.size());

    variables_.push_back(variablesInCategory);

    for (int j=0; j<nVariables_[i];++j) {
        int index = mvaVarMngr_.getVarIndex(variableNamesInCategory[j]);
        if(index == -1) {
            throw cms::Exception("MVA config failure: ")
               << "Concerning " << name_ << tag_ << std::endl
               << "Variable " << variableNamesInCategory[j]
               << " not found in variable definition file!" << std::endl;
        }
        variables_[i].push_back(index);

    }
  }
}

PhotonMVAEstimatorRun2::
~PhotonMVAEstimatorRun2(){
}

float PhotonMVAEstimatorRun2::
mvaValue(const edm::Ptr<reco::Candidate>& particle, const edm::Event& iEvent) const {  

  const int iCategory = findCategory( particle );
  //const std::vector<float> vars = std::move( fillMVAVariables( particle, iEvent ) );  
  const edm::Ptr<reco::Photon> phoRecoPtr = ( edm::Ptr<reco::Photon> )particle;
  std::vector<float> vars;

  for (int i = 0; i < nVariables_[iCategory]; ++i) {
      vars.push_back(mvaVarMngr_.getValue(variables_[iCategory][i], phoRecoPtr, iEvent));
  }
  
  const float response = gbrForests_.at(iCategory)->GetResponse(vars.data());

  if(debug_) {
    std::cout << " *** Inside " << name_ << tag_ << std::endl;
    std::cout << " category " << iCategory << std::endl;
    for (int i = 0; i < nVariables_[iCategory]; ++i) {
        std::cout << " " << mvaVarMngr_.getName(variables_[iCategory][i]) << " " << vars[i] << std::endl;
    }
    std::cout << " ### MVA " << response << std::endl << std::endl;
  }
  
  return response;
}

int PhotonMVAEstimatorRun2::findCategory( const edm::Ptr<reco::Candidate>& particle) const {
  
  // Try to cast the particle into a reco particle.
  // This should work for both reco and pat.
  const edm::Ptr<reco::Photon> phoRecoPtr = ( edm::Ptr<reco::Photon> )particle;
  if( phoRecoPtr.isNull() )
    throw cms::Exception("MVA failure: ")
      << " given particle is expected to be reco::Photon or pat::Photon," << std::endl
      << " but appears to be neither" << std::endl;

  float eta = phoRecoPtr->superCluster()->eta();

  //
  // Determine the category
  //
  int  iCategory = UNDEFINED;
  const float ebeeSplit = 1.479; // division between barrel and endcap

  if ( std::abs(eta) < ebeeSplit)  
    iCategory = CAT_EB;

  if (std::abs(eta) >= ebeeSplit) 
    iCategory = CAT_EE;

  return iCategory;
}

void PhotonMVAEstimatorRun2::setConsumes(edm::ConsumesCollector&& cc) const {
  // All tokens for event content needed by this MVA
  // Tags from the variable helper
  for (auto &tag : mvaVarMngr_.getHelperInputTags()) {
      cc.consumes<edm::ValueMap<float>>(tag);
  }
  for (auto &tag : mvaVarMngr_.getGlobalInputTags()) {
      cc.consumes<double>(tag);
  }
}

// Dummy fonction just to make the template happy
std::vector<float> PhotonMVAEstimatorRun2::
fillMVAVariables(const edm::Ptr<reco::Candidate>& particle, const edm::Event& iEvent) const {
  std::vector<float> vars;
  return vars;
}

DEFINE_EDM_PLUGIN(AnyMVAEstimatorRun2Factory,
          PhotonMVAEstimatorRun2,
          "PhotonMVAEstimatorRun2");
