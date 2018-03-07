#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimatorRun2.h"

ElectronMVAEstimatorRun2::ElectronMVAEstimatorRun2(const edm::ParameterSet& conf):
  AnyMVAEstimatorRun2Base(conf),
  tag_(conf.getParameter<std::string>("mvaTag")),
  name_(conf.getParameter<std::string>("mvaName")),
  nCategories_            (conf.getParameter<int>                     ("nCategories")),
  methodName_             ("BDTG method"),
  variableDefinitionFileName_(conf.getParameter<std::string>("variableDefinition"))
{

  const std::vector <std::string> weightFileNames
    = conf.getParameter<std::vector<std::string> >("weightFileNames");

  const std::vector <std::string> categoryCutStrings
    = conf.getParameter<std::vector<std::string> >("categoryCuts");

  if( (int)(categoryCutStrings.size()) != nCategories_ )
    throw cms::Exception("MVA config failure: ")
      << "wrong number of category cuts" << std::endl;

  for (int i = 0; i < nCategories_; ++i) {
      StringCutObjectSelector<reco::GsfElectron> select(categoryCutStrings[i]);
      categoryFunctions_.push_back(select);
  }

  // Initialize GBRForests from weight files
  init(weightFileNames);

  debug_ = conf.getUntrackedParameter<bool>("debug", false);

}

ElectronMVAEstimatorRun2::ElectronMVAEstimatorRun2(
        const std::string &mvaTag, const std::string &mvaName, const bool debug):
  AnyMVAEstimatorRun2Base( edm::ParameterSet() ),
  tag_                    (mvaTag),
  name_                   (mvaName),
  methodName_             ("BDTG method"),
  debug_                  (debug) { }

void ElectronMVAEstimatorRun2::init(const std::vector<std::string> &weightFileNames) {

  edm::FileInPath variableDefinitionFileEdm(variableDefinitionFileName_);

  // Initialize GBRForests
  if( (int)(weightFileNames.size()) != nCategories_ )
    throw cms::Exception("MVA config failure: ")
      << "wrong number of weightfiles" << std::endl;

  gbrForests_.clear();
  // Create a TMVA reader object for each category
  for(int i=0; i<nCategories_; i++){

    std::vector<std::string> variableNamesInCategory;
    std::vector<Variable> variablesInCategory;

    // Use unique_ptr so that all readers are properly cleaned up
    // when the vector clear() is called in the destructor

    edm::FileInPath weightFile( weightFileNames[i] );
    gbrForests_.push_back( GBRForestTools::createGBRForest( weightFile, variableNamesInCategory ) );

    nVariables_.push_back(variableNamesInCategory.size());

    variables_.push_back(variablesInCategory);

    for (int j=0; j<nVariables_[i];++j) {
        std::string formula;
        bool hasLowerClip = false;
        bool hasUpperClip = false;
        float lowerClipValue = 0.;
        float upperClipValue = 0.;
        bool fromVariableHelper = false;

        std::ifstream file(variableDefinitionFileEdm.fullPath());
        std::string line;
        int lineNumberFromVarName = 0;
        bool variableFound = false;
        while (std::getline(file, line))
        {
            if (line.find("#") != std::string::npos) {
                continue;
            }
            if (line.find(variableNamesInCategory[j]) != std::string::npos) {
                variableFound = true;
            }
            if (lineNumberFromVarName == 1) {
                if (line.find("electronMVAVariableHelper") != std::string::npos) {
                    fromVariableHelper = true;
                }
                formula = line;
            }
            if (lineNumberFromVarName == 2) {
                if (line.find("None") == std::string::npos) {
                    hasLowerClip = true;
                    lowerClipValue = ::atof(line.c_str());
                }
            }
            if (lineNumberFromVarName == 3) {
                if (line.find("None") == std::string::npos) {
                    hasUpperClip = true;
                    upperClipValue = ::atof(line.c_str());
                }
                break;
            }
            if (variableFound) {
                ++lineNumberFromVarName;
            }
        }

        if(!variableFound) {
            throw cms::Exception("MVA config failure: ")
               << "Concerning " << name_ << tag_ << std::endl
               << "Variable " << variableNamesInCategory[j]
               << " not found in variable definition file!" << std::endl << "Check "
               << variableDefinitionFileName_ << std::endl;
        }

        variables_[i].push_back(Variable(variableNamesInCategory[j], formula, hasLowerClip, hasUpperClip, lowerClipValue, upperClipValue, fromVariableHelper));
    }
  }
}

ElectronMVAEstimatorRun2::
~ElectronMVAEstimatorRun2(){
}

void ElectronMVAEstimatorRun2::setConsumes(edm::ConsumesCollector&& cc) const {

  // All tokens for event content needed by this MVA

  // Tags from the variable helper
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:convVtxFitProb"));
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:kfchi2"));
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:kfhits"));
  cc.consumes<edm::ValueMap<float>>(edm::InputTag("electronMVAVariableHelper:rho"));
}

float ElectronMVAEstimatorRun2::
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
  const std::vector<float> vars = fillMVAVariables( particle, iEvent, iCategory );
  return mvaValue(iCategory, vars);
}

float ElectronMVAEstimatorRun2::
mvaValue( const edm::Ptr<reco::GsfElectron>& particle, const edm::EventBase & iEvent) const {

  const int iCategory = findCategory( particle );
  const std::vector<float> vars = fillMVAVariables( particle, iEvent, iCategory );
  return mvaValue(iCategory, vars);
}

float ElectronMVAEstimatorRun2::
mvaValue( const int iCategory, const std::vector<float> & vars) const  {
  const float result = gbrForests_.at(iCategory)->GetClassifier(vars.data());

  if(debug_) {
    std::cout << " *** Inside " << name_ << tag_ << std::endl;
    std::cout << " category " << iCategory << std::endl;
    for (int i = 0; i < nVariables_[iCategory]; ++i) {
        std::cout << " " << variables_[iCategory][i].getName() << " " << vars[i] << std::endl;
    }
    std::cout << " ### MVA " << result << std::endl << std::endl;
  }

  return result;
}

int ElectronMVAEstimatorRun2::findCategory( const edm::Ptr<reco::Candidate>& particle) const {

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

int ElectronMVAEstimatorRun2::findCategory( const edm::Ptr<reco::GsfElectron>& eleRecoPtr ) const {
  for (int i = 0; i < nCategories_; ++i) {
      if (categoryFunctions_[i](*eleRecoPtr)) return i;
  }
  return -1;
}

// Dummy fonction just to make the template happy
std::vector<float> ElectronMVAEstimatorRun2::
fillMVAVariables(const edm::Ptr<reco::Candidate>& particle, const edm::Event& iEvent) const {
  std::vector<float> vars;
  return vars;
}

// A function that should work on both pat and reco objects
std::vector<float> ElectronMVAEstimatorRun2::
fillMVAVariables(const edm::Ptr<reco::Candidate>& particle, const edm::Event& iEvent, const int iCategory) const {

  // Try to cast the particle into a reco particle.
  // This should work for both reco and pat.
  const edm::Ptr<reco::GsfElectron> eleRecoPtr = ( edm::Ptr<reco::GsfElectron> )particle;
  if( eleRecoPtr.get() == nullptr ) {
    throw cms::Exception("MVA failure: ")
      << " given particle is expected to be reco::GsfElectron or pat::Electron," << std::endl
      << " but appears to be neither" << std::endl;
  }

  return fillMVAVariables(eleRecoPtr, iEvent, iCategory);
}

// A function that should work on both pat and reco objects
// The EventType will be edm::Event for the VID accessor and edm::EventBase for for the fwlite accessor
template<class EventType>
std::vector<float> ElectronMVAEstimatorRun2::
fillMVAVariables(const edm::Ptr<reco::GsfElectron>& eleRecoPtr, const EventType& iEvent, const int iCategory) const {

  std::vector<float> vars;

  for (int i = 0; i < nVariables_[iCategory]; ++i) {
      vars.push_back(variables_[iCategory][i].getValue(eleRecoPtr, iEvent));
  }

  return vars;

}
