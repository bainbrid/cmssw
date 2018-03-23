// -*- C++ -*-
//
// Package:    RecoEgamma/EgammaTools
// Class:      MVATrainingNtuplizer
//
/**\class MVATrainingNtuplizer MVATrainingNtuplizer.cc RecoEgamma/EgammaTools/plugins/MVATrainingNtuplizer.cc

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

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class MVATrainingNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MVATrainingNtuplizer(const edm::ParameterSet&);
      ~MVATrainingNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      // for AOD case
      edm::EDGetToken src_;
      edm::EDGetToken vertices_;
      edm::EDGetToken pileup_;

      // for miniAOD case
      edm::EDGetToken srcMiniAOD_;
      edm::EDGetToken verticesMiniAOD_;
      edm::EDGetToken pileupMiniAOD_;

      // other
      TTree* tree_;

      MVAVariableManager<reco::GsfElectron> mvaVarMngr_;
      float vars_[200];
      int nVars_;

      //global variables
      int nEvent_, nRun_, nLumi_;
      int genNpu_;
      int vtxN_;

      // config
      const bool isMC_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MVATrainingNtuplizer::MVATrainingNtuplizer(const edm::ParameterSet& iConfig)
 :
  src_            (consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("src"))),
  vertices_       (consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
  pileup_         (consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileup"))),
  srcMiniAOD_     (consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("srcMiniAOD"))),
  verticesMiniAOD_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("verticesMiniAOD"))),
  pileupMiniAOD_  (consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileupMiniAOD"))),
  mvaVarMngr_     (iConfig.getParameter<std::string>("variableDefinition")),
  isMC_           (iConfig.getParameter<bool>("isMC"))

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
   tree_->Branch("gen_npu", &genNpu_);
   tree_->Branch("vtx_N",   &vtxN_);

   for (int i = 0; i < nVars_; ++i) {
       tree_->Branch(mvaVarMngr_.getName(i).c_str(), &(vars_[i]));
   }

   // All tokens for event content needed by this MVA
   // Tags from the variable helper
   for (auto &tag : mvaVarMngr_.getHelperInputTags()) {
       consumes<edm::ValueMap<float>>(tag);
   }
   for (auto &tag : mvaVarMngr_.getGlobalInputTags()) {
       consumes<double>(tag);
   }
}


MVATrainingNtuplizer::~MVATrainingNtuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MVATrainingNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

        for (int iVar = 0; iVar < nVars_; ++iVar) {
            vars_[iVar] = mvaVarMngr_.getValue(iVar, ele, iEvent);
        }

        tree_->Fill();
    }

}


// ------------ method called once each job just before starting event loop  ------------
void
MVATrainingNtuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MVATrainingNtuplizer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MVATrainingNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(MVATrainingNtuplizer);
