#ifndef LowPtGsfElectronProducer_h
#define LowPtGsfElectronProducer_h

#include "RecoEgamma/EgammaElectronProducers/plugins/GsfElectronBaseProducer.h"

class LowPtGsfElectronProducer : public GsfElectronBaseProducer {
  
 public:
  
  explicit LowPtGsfElectronProducer( const edm::ParameterSet&, 
				     const gsfAlgoHelpers::HeavyObjectCache* );

  ~LowPtGsfElectronProducer();

  void produce( edm::Event&, const edm::EventSetup& );

 private:

};

#endif // LowPtGsfElectronProducer_h
