#ifndef BACONPROD_NTUPLER_FILLERGENJET_HH
#define BACONPROD_NTUPLER_FILLERGENJET_HH

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include <vector>
#include <string>

// forward class declarations
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TClonesArray;

namespace baconhep
{
  class FillerGenJets
  {
    public:
       FillerGenJets(const edm::ParameterSet &iConfig);
      ~FillerGenJets();
    void    fill(TClonesArray *array,TClonesArray *fatJetArray,const edm::Event &iEvent);
    protected:
       double* genCone(const reco::GenJet *iJet,const reco::GenParticleCollection &iGenParticles,double iDRMin,double iDRMax,int iType);
       int     flavor (const reco::GenJet *iJet,const reco::GenParticleCollection &iGenParticles);            
       void    trim(const reco::GenJet *iJet,float &iMTrim,float &iTau1,float &iTau2);
      // EDM object collection names
      std::string fGenParName;
      std::string fGenJetName;
      std::string fGenFatJetName;

      fastjet::Filter* fTrimmer1;
      fastjet::JetDefinition*       fCAJetDef;
      fastjet::ActiveAreaSpec*      fActiveArea;
      fastjet::AreaDefinition*      fAreaDefinition;
      fastjet::ClusterSequenceArea* fClustering;

  };
}
#endif
