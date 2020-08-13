import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.slimming.slimmedLowPtElectrons_cfi import *

from Configuration.Eras.Modifier_bParking_cff import bParking
bParking.toModify(lowPtGsfElectronID, 
                  electrons='slimmedLowPtElectrons'
                  rho='fixedGridRhoFastjetAll')
