// -*- C++ -*-
//
// Package:    RecoEgamma/ElectronIdentification
// Class:      ElectronMVATrainingNtuplizer
//
/**\class ElectronMVATrainingNtuplizer ElectronMVATrainingNtuplizer.cc RecoEgamma/ElectronIdentification/plugins/ElectronMVATrainingNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jonas REMBSER
//         Created:  Thu, 22 Mar 2018 14:54:24 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "RecoEgamma/EgammaTools/interface/MVAVariableManager.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//

class ElectronMVATrainingNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronMVATrainingNtuplizer(const edm::ParameterSet&);
      ~ElectronMVATrainingNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void findFirstNonElectronMother2(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);

      template<class T, class V>
      int matchToTruth(const T &el, const V &genParticles, int &genIdx);

      // ----------member data ---------------------------

      // for AOD case
      edm::EDGetToken src_;
      edm::EDGetToken vertices_;
      edm::EDGetToken pileup_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticles_;

      // for miniAOD case
      edm::EDGetToken srcMiniAOD_;
      edm::EDGetToken verticesMiniAOD_;
      edm::EDGetToken pileupMiniAOD_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAOD_;

      // other
      TTree* tree_;

      MVAVariableManager<reco::GsfElectron> mvaVarMngr_;
      float vars_[200];
      int nVars_;

      //global variables
      int nEvent_, nRun_, nLumi_;
      int genNpu_;
      int vtxN_;
      int matchedToGenEle_;
      int matchedGenIdx_;

      // gap variables
      bool eleIsEB_;
      bool eleIsEE_;
      bool eleIsEBEtaGap_;
      bool eleIsEBPhiGap_;
      bool eleIsEBEEGap_;
      bool eleIsEEDeeGap_;
      bool eleIsEERingGap_;

      // config
      const bool isMC_;
      const double deltaR_;
      const double ptThreshold_;
};

//
// constants, enums and typedefs
//

enum ElectronMatchType {
                        UNMATCHED,
                        TRUE_PROMPT_ELECTRON,
                        TRUE_ELECTRON_FROM_TAU,
                        TRUE_NON_PROMPT_ELECTRON,
                       }; // The last does not include tau parents

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronMVATrainingNtuplizer::ElectronMVATrainingNtuplizer(const edm::ParameterSet& iConfig)
 :
  src_                   (consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("src"))),
  vertices_              (consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
  pileup_                (consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileup"))),
  genParticles_          (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  srcMiniAOD_            (consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("srcMiniAOD"))),
  verticesMiniAOD_       (consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("verticesMiniAOD"))),
  pileupMiniAOD_         (consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileupMiniAOD"))),
  genParticlesMiniAOD_   (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticlesMiniAOD"))),
  mvaVarMngr_            (iConfig.getParameter<std::string>("variableDefinition")),
  isMC_                  (iConfig.getParameter<bool>("isMC")),
  deltaR_                (iConfig.existsAs<double>("deltaR")        ? iConfig.getParameter<double>("deltaR"): 0.1),
  ptThreshold_           (iConfig.existsAs<double>("ptThreshold")   ? iConfig.getParameter<double>("ptThreshold"): 5)

{
   //now do what ever initialization is needed

   // Book tree
   usesResource(TFileService::kSharedResource);
   edm::Service<TFileService> fs ;
   tree_  = fs->make<TTree>("tree","tree");

   nVars_ = mvaVarMngr_.getNVars();

   tree_->Branch("nEvent",  &nEvent_);
   tree_->Branch("nRun",    &nRun_);
   tree_->Branch("nLumi",   &nLumi_);
   tree_->Branch("genNpu", &genNpu_);
   tree_->Branch("vtxN",   &vtxN_);

   if (isMC_) {
       tree_->Branch("matchedToGenEle",   &matchedToGenEle_);
   }

   for (int i = 0; i < nVars_; ++i) {
       tree_->Branch(mvaVarMngr_.getName(i).c_str(), &(vars_[i]));
   }

   tree_->Branch("ele_isEB",&eleIsEB_);
   tree_->Branch("ele_isEE",&eleIsEE_);
   tree_->Branch("ele_isEBEtaGap",&eleIsEBEtaGap_);
   tree_->Branch("ele_isEBPhiGap",&eleIsEBPhiGap_);
   tree_->Branch("ele_isEBEEGap", &eleIsEBEEGap_);
   tree_->Branch("ele_isEEDeeGap",&eleIsEEDeeGap_);
   tree_->Branch("ele_isEERingGap",&eleIsEERingGap_);

   // All tokens for event content needed by this MVA
   // Tags from the variable helper
   for (auto &tag : mvaVarMngr_.getHelperInputTags()) {
       consumes<edm::ValueMap<float>>(tag);
   }
   for (auto &tag : mvaVarMngr_.getGlobalInputTags()) {
       consumes<double>(tag);
   }
}


ElectronMVATrainingNtuplizer::~ElectronMVATrainingNtuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronMVATrainingNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Fill global event info
    nEvent_ = iEvent.id().event();
    nRun_   = iEvent.id().run();
    nLumi_  = iEvent.luminosityBlock();


    // Retrieve Vertecies
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertices_, vertices);
    if( !vertices.isValid() ){
      iEvent.getByToken(verticesMiniAOD_,vertices);
      if( !vertices.isValid() )
        throw cms::Exception(" Collection not found: ")
          << " failed to find a standard AOD or miniAOD vertex collection " << std::endl;
    }

    vtxN_ = vertices->size();

    // Retrieve Pileup Info
    edm::Handle<std::vector< PileupSummaryInfo > >  pileup;
    iEvent.getByToken(pileup_, pileup);
    if( !pileup.isValid() ){
      iEvent.getByToken(pileupMiniAOD_,pileup);
      if( !pileup.isValid() )
        throw cms::Exception(" Collection not found: ")
          << " failed to find a standard AOD or miniAOD pileup collection " << std::endl;
    }

    // Fill with true number of pileup
    if(isMC_) {
       for(const auto& pu : *pileup)
       {
           int bx = pu.getBunchCrossing();
           if(bx == 0)
           {
               genNpu_ = pu.getPU_NumInteractions();
               break;
           }
       }
    }

    // Retrieve genParticles
    edm::Handle<edm::View<reco::GenParticle> >  genParticles;
    iEvent.getByToken(genParticles_, genParticles);
    if( !genParticles.isValid() ){
      iEvent.getByToken(genParticlesMiniAOD_, genParticles);
      if( !genParticles.isValid() )
        throw cms::Exception(" Collection not found: ")
          << " failed to find a standard AOD or miniAOD genParticle collection " << std::endl;
    }


    edm::Handle<edm::View<reco::GsfElectron> > src;

    // Retrieve the collection of particles from the event.
    // If we fail to retrieve the collection with the standard AOD
    // name, we next look for the one with the stndard miniAOD name.
    iEvent.getByToken(src_, src);
    if( !src.isValid() ){
      iEvent.getByToken(srcMiniAOD_,src);
      if( !src.isValid() )
        throw cms::Exception(" Collection not found: ")
          << " failed to find a standard AOD or miniAOD particle collection " << std::endl;
    }

    int nEle = src->size();

    for(int iEle = 0; iEle < nEle; ++iEle) {

        const auto ele =  src->ptrAt(iEle);

        if (ele->pt() < ptThreshold_) {
            continue;
        }

        for (int iVar = 0; iVar < nVars_; ++iVar) {
            vars_[iVar] = mvaVarMngr_.getValue(iVar, ele, iEvent);
        }

        if (isMC_) {
            matchedToGenEle_ = matchToTruth( ele, genParticles, matchedGenIdx_);
        }

        // gap variables
        eleIsEB_ = ele->isEB();
        eleIsEE_ = ele->isEE();
        eleIsEBEEGap_ = ele->isEBEEGap();
        eleIsEBEtaGap_ = ele->isEBEtaGap();
        eleIsEBPhiGap_ = ele->isEBPhiGap();
        eleIsEEDeeGap_ = ele->isEEDeeGap();
        eleIsEERingGap_ = ele->isEERingGap();

        tree_->Fill();
    }

}

void ElectronMVATrainingNtuplizer::findFirstNonElectronMother2(const reco::Candidate *particle,
                         int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother2(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

template<class T, class V>
int ElectronMVATrainingNtuplizer::matchToTruth(const T &el, const V &prunedGenParticles, int &genIdx){

  //
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
      genIdx = i;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < deltaR_) ) {
    return UNMATCHED;
  }

  //
  int ancestorPID = -999;
  int ancestorStatus = -999;
  findFirstNonElectronMother2(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }

  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

// ------------ method called once each job just before starting event loop  ------------
void
ElectronMVATrainingNtuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronMVATrainingNtuplizer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronMVATrainingNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMVATrainingNtuplizer);
