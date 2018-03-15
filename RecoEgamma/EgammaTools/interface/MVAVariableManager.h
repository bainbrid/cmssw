#ifndef RecoEgamma_EgammaTools_MVAVariableManager_H
#define RecoEgamma_EgammaTools_MVAVariableManager_H

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

using namespace std;

template <class ParticleType>
class MVAVariableManager {

  public:
    MVAVariableManager() {
        nVars = 0;
    };

    MVAVariableManager(const string variableDefinitionFileName) {
        init(variableDefinitionFileName);
    };

    int init(const string variableDefinitionFileName) {
        nVars = 0;

        variableInfos_.clear();
        functions_.clear();
        formulas_.clear();
        names_.clear();
        helperInputTags_.clear();
        globalInputTags_.clear();

        edm::FileInPath variableDefinitionFileEdm(variableDefinitionFileName);
        ifstream file(variableDefinitionFileEdm.fullPath());

        string name, formula, upper, lower;
        while( true ) {
            file >> name;
            if (name.find("#") != string::npos) {
                file.ignore(numeric_limits<streamsize>::max(), '\n');
                continue;
            }
            file >> formula >> lower >> upper;
            if (file.eof()) {
                break;
            }
            addVariable_(name, formula, lower, upper);
        }
        return nVars;
    };

    int getVarIndex(string &name) {
        map<string,int>::iterator it = indexMap_.find(name);
        if (it == indexMap_.end()) {
            return -1;
        } else {
            return it->second;
        }
    }

    const string getName(int index) const {
        return names_[index];
    }

    vector<edm::InputTag> getHelperInputTags() const {
        return helperInputTags_;
    }

    vector<edm::InputTag> getGlobalInputTags() const {
        return globalInputTags_;
    }

    template <class EventType>
    float getValue(int index, const edm::Ptr<ParticleType>& ptclPtr, const EventType& iEvent) const {
        float value;
        MVAVariableInfo varInfo = variableInfos_[index];
        if (varInfo.fromVariableHelper) {
            edm::Handle<edm::ValueMap<float>> vMap;
            iEvent.getByLabel(edm::InputTag(formulas_[index]), vMap);
            value = (*vMap)[ptclPtr];
        } else if (varInfo.isGlobalVariable) {
            edm::Handle<double> valueHandle;
            iEvent.getByLabel(edm::InputTag(formulas_[index]), valueHandle);
            value = *valueHandle;
        } else {
            value = functions_[index](*ptclPtr);
        }
        if (varInfo.hasLowerClip && value < varInfo.lowerClipValue) {
            value = varInfo.lowerClipValue;
        }
        if (varInfo.hasUpperClip && value > varInfo.upperClipValue) {
            value = varInfo.upperClipValue;
        }
        return value;
    }

  private:

    struct MVAVariableInfo {
        bool hasLowerClip;
        bool hasUpperClip;
        float lowerClipValue;
        float upperClipValue;
        bool fromVariableHelper;
        bool isGlobalVariable;
    };

    void addVariable_(string &name, string &formula, string &lowerClip, string &upperClip) {
        bool hasLowerClip = lowerClip.find("None") == string::npos;
        bool hasUpperClip = upperClip.find("None") == string::npos;
        bool fromVariableHelper = formula.find("MVAVariableHelper") != string::npos ||
                                  formula.find("IDValueMapProducer") != string::npos ||
                                  formula.find("egmPhotonIsolation") != string::npos;
        float lowerClipValue = hasLowerClip ? (float)::atof(lowerClip.c_str()) : 0.;
        float upperClipValue = hasUpperClip ? (float)::atof(upperClip.c_str()) : 0.;

        // fixedGridRhoFastjetAll is the only global variable used ever, so its hardcoded...
        bool isGlobalVariable = formula.find("Rho") != string::npos;

        functions_.push_back(!(fromVariableHelper || isGlobalVariable) ? StringObjectFunction<ParticleType>(formula) : StringObjectFunction<ParticleType>("pt"));
        formulas_.push_back(formula);
        if (fromVariableHelper) {
            helperInputTags_.push_back(edm::InputTag(formula));
        }
        if (isGlobalVariable) {
            globalInputTags_.push_back(edm::InputTag(formula));
        }
        MVAVariableInfo varInfo = {
            .hasLowerClip       = hasLowerClip,
            .hasUpperClip       = hasUpperClip,
            .lowerClipValue     = lowerClipValue,
            .upperClipValue     = upperClipValue,
            .fromVariableHelper = fromVariableHelper,
            .isGlobalVariable   = isGlobalVariable};
        variableInfos_.push_back(varInfo);
        names_.push_back(name);
        indexMap_[name] = nVars;
        nVars++;
    };


    int nVars;
    vector<MVAVariableInfo> variableInfos_;
    vector<StringObjectFunction<ParticleType>> functions_;
    vector<string> formulas_;
    vector<string> names_;
    map<string, int> indexMap_;

    // To store the MVAVariableHelper input tags needed for the variables in this container
    vector<edm::InputTag> helperInputTags_;

    vector<edm::InputTag> globalInputTags_;
};

#endif
