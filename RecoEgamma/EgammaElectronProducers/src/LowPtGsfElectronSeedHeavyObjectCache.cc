#include "RecoEgamma/EgammaElectronProducers/interface/LowPtGsfElectronSeedHeavyObjectCache.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Reader.h"
#include <string>

namespace lowptgsfeleseed {
  HeavyObjectCache::HeavyObjectCache( const edm::ParameterSet& conf ) {
    
    edm::FileInPath weights( conf.getParameter<edm::FileInPath>("Weights") );
    
    TMVA::Reader reader("!Color:Silent");
    reader.AddVariable("trk_pt",&trk_pt);
    reader.AddVariable("trk_eta",&trk_eta);
    reader.AddVariable("trk_phi",&trk_phi);
    reader.AddVariable("trk_p",&trk_p);
    reader.AddVariable("trk_charge",&trk_charge);
    reader.AddVariable("trk_nhits",&trk_nhits);
    reader.AddVariable("trk_high_purity",&trk_high_purity);
    reader.AddVariable("trk_inp",&trk_inp);
    reader.AddVariable("trk_outp",&trk_outp);
    reader.AddVariable("trk_chi2red",&trk_chi2red);
    reader.AddVariable("preid_trk_ecal_Deta",&preid_trk_ecal_Deta);
    reader.AddVariable("preid_trk_ecal_Dphi",&preid_trk_ecal_Dphi);
    reader.AddVariable("preid_e_over_p",&preid_e_over_p);
    reader.AddVariable("preid_gsf_dpt",&preid_gsf_dpt);
    reader.AddVariable("preid_trk_gsf_chiratio",&preid_trk_gsf_chiratio);
    reader.AddVariable("preid_gsf_chi2red",&preid_gsf_chi2red);
    reader.BookMVA("BDT", weights.fullPath().c_str());
    
    gbr.reset( new GBRForest(dynamic_cast<TMVA::MethodBDT*>(reader.FindMVA("BDT"))) );
    
  }
}
