#include "BaconProd/Ntupler/interface/FillerJet.hh"
#include "BaconProd/Utils/interface/EnergyCorrelator.hh"
#include "BaconProd/Utils/interface/TriggerTools.hh"
#include "BaconProd/Utils/interface/JetTools.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
//#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h" 
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "RecoBTag/SecondaryVertex/interface/TrackSelector.h"
#include "RecoBTag/SecondaryVertex/interface/V0Filter.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <iostream>

using namespace baconhep;

class tauAxes
{
public:
  tauAxes (const fastjet::PseudoJet tauAxis, float dRSV)
  {
    dR = dRSV;
    tau = tauAxis;
  }
  ~tauAxes() { }

  float dR;
  fastjet::PseudoJet tau;
};

//--------------------------------------------------------------------------------------------------
bool compAxesSV(tauAxes tauAxis1, tauAxes tauAxis2) {
  return tauAxis1.dR<tauAxis2.dR;
}

//--------------------------------------------------------------------------------------------------
FillerJet::FillerJet(const edm::ParameterSet &iConfig, const bool useAOD,edm::ConsumesCollector && iC):
  fMinPt              (iConfig.getUntrackedParameter<double>("minPt",20)),
  fUseGen             (iConfig.getUntrackedParameter<bool>("doGenJet",true)),
  fApplyJEC           (iConfig.getUntrackedParameter<bool>("applyJEC",true)),
  fPVName             (iConfig.getUntrackedParameter<std::string>("edmPVName","offlinePrimaryVertices")),
  fRhoName            (iConfig.getUntrackedParameter<std::string>("edmRhoName","fixedGridRhoFastjetAll")),
  fJetName            (iConfig.getUntrackedParameter<std::string>("jetName","ak4PFJetsCHS")),
  fGenJetName         (iConfig.getUntrackedParameter<std::string>("genJetName","AK4GenJetsCHS")),
  fJetFlavorName      (iConfig.getUntrackedParameter<std::string>("jetFlavorName","AK4byValAlgoCHS")),
  fPrunedJetName      (iConfig.getUntrackedParameter<std::string>("prunedJetName","AK4caPFJetsPrunedCHS")),
  fTrimmedJetName     (iConfig.getUntrackedParameter<std::string>("trimmedJetName","AK4caPFJetsTrimmedCHS")),
  fSoftDropJetName    (iConfig.getUntrackedParameter<std::string>("softdropJetName","AK4caPFJetsSoftDropCHS")),
  fSubJetName         (iConfig.getUntrackedParameter<std::string>("subJetName","AK4caPFJetsTrimmedCHS__SubJets")),
  fCSVbtagName        (iConfig.getUntrackedParameter<std::string>("csvBTagName","combinedInclusiveSecondaryVertexV2BJetTags")),
  fCSVbtagSubJetName  (iConfig.getUntrackedParameter<std::string>("csvBTagSubJetName","AK4CombinedInclusiveSecondaryVertexV2BJetTagsSJCHS")),
  fCSVDoubleBtagName  (iConfig.getUntrackedParameter<std::string>("csvDoubleBTagName","AK4PFBoostedDoubleSecondaryVertexBJetTagsCHS")),
  fJettinessName      (iConfig.getUntrackedParameter<std::string>("jettiness","AK4NjettinessCHS")),
  fIPTagInfos         (iConfig.getUntrackedParameter<std::string>("ipTagInfos","AK4ImpactParameterTagInfosCHS")),
  fSVTagInfos         (iConfig.getUntrackedParameter<std::string>("svTagInfos","AK4InclusiveSecondaryVertexFinderTagInfosCHS")),
  fQGLikelihood       (iConfig.getUntrackedParameter<std::string>("qgLikelihood","QGLikelihood")),
  fQGLikelihoodSubJets(iConfig.getUntrackedParameter<std::string>("qgLikelihoodSubjet","QGLikelihood")),
  fTopTaggerName      (iConfig.getUntrackedParameter<std::string>("topTaggerName","")),
  fShowerDecoConf     (iConfig.getUntrackedParameter<std::string>("showerDecoConf","")),
  fConeSize           (iConfig.getUntrackedParameter<double>("coneSize",0.4)),
  fComputeFullJetInfo (iConfig.getUntrackedParameter<bool>("doComputeFullJetInfo",false)),  
  fShowerDeco         (0),
  fJetCorr            (0),
  fJetUnc             (0),
  fUseAOD             (useAOD),
  fbeta               (iConfig.getUntrackedParameter<double>("beta_")),
  fR0                 (iConfig.getUntrackedParameter<double>("R0_")),
  fnjettiness         (fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(fbeta,fR0)),
  ftrackPairV0Filter  (iConfig.getUntrackedParameter<edm::ParameterSet>("trackPairV0Filter")),
  ftrackSelector      (iConfig.getParameter<edm::ParameterSet>("trackSelection")),
  fmaxSVDeltaRToJet   (iConfig.getUntrackedParameter<double>("maxSVDeltaRToJet"))
{
    std::vector<std::string> empty_vstring;
    initJetCorr(iConfig.getUntrackedParameter< std::vector<std::string> >("jecFiles",empty_vstring),
                iConfig.getUntrackedParameter< std::vector<std::string> >("jecUncFiles",empty_vstring));

    std::string cmssw_base_src = getenv("CMSSW_BASE");
    cmssw_base_src += "/src/";

    if(fUseAOD) {
      std::vector<std::string> puIDFiles = iConfig.getUntrackedParameter< std::vector<std::string> >("jetPUIDFiles",empty_vstring);
      assert(puIDFiles.size()==2);
      fLowPtWeightFile  = (puIDFiles[0].length()>0) ? (cmssw_base_src + puIDFiles[0]) : "";
      fHighPtWeightFile = (puIDFiles[1].length()>0) ? (cmssw_base_src + puIDFiles[1]) : "";
      initPUJetId();
    }

  fRand = new TRandom2();
  if(fShowerDecoConf.size() > 0) { 
    fShowerDeco = new ShowerDeco(cmssw_base_src+fShowerDecoConf);
  }
  fTokRhoTag        = iC.consumes<double>               (fRhoName);
  if(fUseAOD)  fTokJetName       = iC.consumes<reco::PFJetCollection>(fJetName);
  if(!fUseAOD) fTokPatJetName    = iC.consumes<pat::JetCollection>   (fJetName);
  fTokGenJetName    = iC.consumes<reco::GenJetCollection>(fGenJetName);
  fTokJetFlavorName = iC.consumes<reco::JetFlavourInfoMatchingCollection>(fJetFlavorName);
  fTokPVName        = iC.consumes<reco::VertexCollection>(fPVName);
  fTokCSVbtagName   = iC.consumes<reco::JetTagCollection>(fCSVbtagName);
  edm::InputTag lQGLikelihood(fQGLikelihood,"qgLikelihood");
  edm::InputTag lQGLAxis2    (fQGLikelihood,"axis2");
  edm::InputTag lQGLPtD      (fQGLikelihood,"ptD");
  edm::InputTag lQGLMult     (fQGLikelihood,"mult");
  fTokQGLikelihood        = iC.consumes<edm::ValueMap<float> >   (lQGLikelihood);
  fTokQGLAxis2            = iC.consumes<edm::ValueMap<float> >   (lQGLAxis2);
  fTokQGLPtD              = iC.consumes<edm::ValueMap<float> >   (lQGLPtD);
  fTokQGLMult             = iC.consumes<edm::ValueMap<int> >     (lQGLMult);
  if(fComputeFullJetInfo) { 
    fTokPrunedJetName     = iC.consumes<reco::BasicJetCollection>                  (fPrunedJetName);
    fTokTrimmedJetName    = iC.consumes<reco::BasicJetCollection>                  (fTrimmedJetName);
    fTokSoftDropJetName   = iC.consumes<reco::BasicJetCollection>                  (fSoftDropJetName);
    fTokCSVbtagSubJetName = iC.consumes<reco::JetTagCollection>                    (fCSVbtagSubJetName);
    fTokCSVDoubleBtagName = iC.consumes<reco::JetTagCollection>                    (fCSVDoubleBtagName);
    fTokIPTagInfos        = iC.consumes<reco::CandIPTagInfoCollection>             (fIPTagInfos);
    fTokSVTagInfos        = iC.consumes<reco::CandSecondaryVertexTagInfoCollection>(fSVTagInfos);
    edm::InputTag lTau1(fJettinessName,"tau1");
    edm::InputTag lTau2(fJettinessName,"tau2");
    edm::InputTag lTau3(fJettinessName,"tau3");
    edm::InputTag lTau4(fJettinessName,"tau4");
    edm::InputTag lQGSubJets(fQGLikelihoodSubJets,"qgLikelihood");
    edm::InputTag lSubJets  (fSubJetName,"SubJets");
    edm::InputTag lTopTagSubJet(fTopTaggerName,"caTopSubJets"); //"cmsTopTagPFJetsCHS"
    edm::InputTag lTopTag      (fTopTaggerName);

    fTokTau1Name           = iC.consumes<edm::ValueMap<float> >   (lTau1);
    fTokTau2Name           = iC.consumes<edm::ValueMap<float> >   (lTau2);
    fTokTau3Name           = iC.consumes<edm::ValueMap<float> >   (lTau3);
    fTokTau4Name           = iC.consumes<edm::ValueMap<float> >   (lTau4);
    fTokQGLSubJets         = iC.consumes<edm::ValueMap<float> >   (lQGSubJets);
    fTokSubJets            = iC.consumes<reco::PFJetCollection>   (lSubJets);
    fTokCMSTTJetProduct    = iC.consumes<reco::BasicJetCollection>(lTopTag);
    fTokCMSTTSubJetProduct = iC.consumes<reco::PFJetCollection>   (lTopTagSubJet);
  }
}

//--------------------------------------------------------------------------------------------------
FillerJet::~FillerJet()
{
  delete fJetCorr;
  delete fJetUnc;
}
void FillerJet::initPUJetId() { 
  if(!fUseAOD) return;
  std::cout << "===> Re-initializing" << std::endl;
  fJetPUIDMVACalc.initialize(baconhep::JetPUIDMVACalculator::k53,
			     "BDT",fLowPtWeightFile,
			     "BDT",fHighPtWeightFile);
  
}
//--------------------------------------------------------------------------------------------------
void FillerJet::initJetCorr(const std::vector<std::string> &jecFiles,
                            const std::vector<std::string> &jecUncFiles)
{
  assert(jecFiles.size()>0);
  assert(jecUncFiles.size()>0);

  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
 
  std::vector<JetCorrectorParameters> corrParams;
  for(unsigned int icorr=0; icorr<jecFiles.size(); icorr++) {
    corrParams.push_back(JetCorrectorParameters(cmssw_base_src + jecFiles[icorr]));
  }
  fJetCorr = new FactorizedJetCorrector(corrParams);
  
  JetCorrectorParameters param(cmssw_base_src + jecUncFiles[0]);
  fJetUnc = new JetCorrectionUncertainty(param);
}

//--------------------------------------------------------------------------------------------------
// === filler for AOD ===
void FillerJet::fill(TClonesArray *array, TClonesArray *iExtraArray,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup, 
		     const reco::Vertex	&pv,
		     const std::vector<TriggerRecord> &triggerRecords,
		     const trigger::TriggerEvent *triggerEvent,
                     const pat::TriggerObjectStandAloneCollection *patTriggerObjects)
{
  assert(array);
  assert(!fComputeFullJetInfo || iExtraArray);
  if(fUseAOD) { assert(fJetPUIDMVACalc.isInitialized()); }
  fRand->SetSeed(iEvent.id().event());
  
  // Get jet collection
  edm::Handle<reco::PFJetCollection> hJetProduct;
  iEvent.getByToken(fTokJetName,hJetProduct);
  assert(hJetProduct.isValid());
  const reco::PFJetCollection *jetCol = hJetProduct.product();
  
  // Get gen jet collection
  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJetCol = 0;
  if(fUseGen) { 
    iEvent.getByToken(fTokGenJetName,hGenJetProduct);
    assert(hGenJetProduct.isValid());
    genJetCol = hGenJetProduct.product();
  }

  // Get Jet Flavor Match
  edm::Handle<reco::JetFlavourInfoMatchingCollection> hJetFlavourMatch;
  if(fUseGen) {
    iEvent.getByToken(fTokJetFlavorName, hJetFlavourMatch);
    assert(hJetFlavourMatch.isValid());
  }

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(fTokPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid()); 
 
  // Get b-jet tagger
  edm::Handle<reco::JetTagCollection> hCSVbtags;
  iEvent.getByToken(fTokCSVbtagName, hCSVbtags);
  assert(hCSVbtags.isValid());

  //Get Quark Gluon Likelihood
  edm::Handle<edm::ValueMap<float> > hQGLikelihood; 
  iEvent.getByToken(fTokQGLikelihood,hQGLikelihood); 
  assert(hQGLikelihood.isValid());

  edm::Handle<edm::ValueMap<float> > hQGLaxis2;
  iEvent.getByToken(fTokQGLAxis2,hQGLaxis2);
  assert(hQGLaxis2.isValid());

  edm::Handle<edm::ValueMap<float> > hQGLptD;
  iEvent.getByToken(fTokQGLPtD,hQGLptD);
  assert(hQGLptD.isValid());

  edm::Handle<edm::ValueMap<int> > hQGLmult;
  iEvent.getByToken(fTokQGLMult,hQGLmult);
  assert(hQGLmult.isValid());

  // Get track builder setup                                                                                                                                                                                                     
  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);
  theTTBuilder = trackBuilder.product();

  TClonesArray &rArray      = *array;
  TClonesArray &rExtraArray = *iExtraArray;
  for(reco::PFJetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {
    const double ptRaw = itJet->pt();

    // input to jet corrections
    double jetcorr = 1;
    if(fabs(itJet->eta()) < 5.191 && fApplyJEC) {
      fJetCorr->setJetPt(ptRaw);
      fJetCorr->setJetEta(itJet->eta());
      fJetCorr->setJetPhi(itJet->phi());
      fJetCorr->setJetE(itJet->energy());
      fJetCorr->setRho(*hRho);
      fJetCorr->setJetA(itJet->jetArea());
      fJetCorr->setJetEMF(-99.0);
      jetcorr = fJetCorr->getCorrection();
    }

    // jet pT cut (BOTH raw AND corrected pT must exceed threshold)
    if(ptRaw*jetcorr < fMinPt || ptRaw < fMinPt) continue;
    
    fJetUnc->setJetPt ( ptRaw  );
    fJetUnc->setJetEta( itJet->eta() );
    double jetunc = fJetUnc->getUncertainty(true);
    bool passLoose = JetTools::passPFLooseID(*itJet);
    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TJet();
    baconhep::TJet *pJet = (baconhep::TJet*)rArray[index];
 
    //
    // Kinematics
    //==============================    
    pJet->pt    = ptRaw * jetcorr;
    pJet->eta   = itJet->eta();
    pJet->phi   = itJet->phi();
    pJet->mass  = itJet->mass() * jetcorr;
    pJet->ptRaw = ptRaw;
    pJet->area  = itJet->jetArea();
    pJet->unc   = jetunc;
    //
    // Impact Parameter
    //==============================
    pJet->d0 = JetTools::jetD0(*itJet, pv);
    pJet->dz = JetTools::jetDz(*itJet, pv);

    //
    // Identification
    //==============================
    reco::PFJetRef jetRef(hJetProduct, itJet - jetCol->begin());
    reco::JetBaseRef jetBaseRef(jetRef);
    pJet->csv  = (*(hCSVbtags.product()))[jetBaseRef];

    pJet->qgid  = (*(hQGLikelihood.product()))[jetBaseRef];
    pJet->axis2 = (*(hQGLaxis2.product()))[jetBaseRef];
    pJet->ptD   = (*(hQGLptD.product()))[jetBaseRef];
    pJet->mult  = (*(hQGLmult.product()))[jetBaseRef];

    pJet->q = JetTools::jetCharge(*itJet);
    
    pJet->beta     = JetTools::beta(*itJet, pv);
    pJet->betaStar = JetTools::betaStar(*itJet, pv, pvCol);
    pJet->dR2Mean  = JetTools::dR2Mean(*itJet);
    pJet->mva = -2;
    if(passLoose) {
      double dRMean = JetTools::dRMean(*itJet);
      double frac01 = JetTools::frac(*itJet,0.1);
      double frac02 = JetTools::frac(*itJet,0.2);
      double frac03 = JetTools::frac(*itJet,0.3);
      double frac04 = JetTools::frac(*itJet,0.4);
      double frac05 = JetTools::frac(*itJet,0.5);
      pJet->mva = fJetPUIDMVACalc.mvaValue((float)pvCol->size(), ptRaw*jetcorr, itJet->eta(), itJet->phi(),
                                           pJet->d0, pJet->dz, pJet->beta, pJet->betaStar, itJet->chargedMultiplicity(), itJet->neutralMultiplicity(),
					   dRMean, pJet->dR2Mean, pJet->ptD, frac01, frac02, frac03, frac04, frac05);
    }

    TVector2 lPull = JetTools::jetPull(*itJet,0);
    pJet->pullY      = lPull.X();
    pJet->pullPhi    = lPull.Y();
    TVector2 lChPull = JetTools::jetPull(*itJet,1);
    pJet->chPullY    = lChPull.X();
    pJet->chPullPhi  = lChPull.Y();
    TVector2 lNeuPull = JetTools::jetPull(*itJet,2);
    pJet->neuPullY   = lNeuPull.X();
    pJet->neuPullPhi = lNeuPull.Y();

    // Basic Noise Variables
    pJet->chEmFrac   = itJet->chargedEmEnergy() / itJet->energy();
    pJet->neuEmFrac  = itJet->neutralEmEnergy() / itJet->energy();
    pJet->chHadFrac  = itJet->chargedHadronEnergy() / itJet->energy();
    pJet->neuHadFrac = itJet->neutralHadronEnergy() / itJet->energy();
    pJet->muonFrac   = itJet->muonEnergy() / itJet->energy();

    //
    // Generator matching
    //==============================
    const reco::GenJet * matchGenJet = 0; 
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    if(matchGenJet != 0) { 
      pJet->partonFlavor = (*hJetFlavourMatch)[jetBaseRef].getPartonFlavour();
      pJet->hadronFlavor = (*hJetFlavourMatch)[jetBaseRef].getHadronFlavour();
      pJet->genpt        = matchGenJet->pt();
      pJet->geneta       = matchGenJet->eta();
      pJet->genphi       = matchGenJet->phi();
      pJet->genm         = matchGenJet->mass();
    }

    pJet->nCharged   = itJet->chargedMultiplicity();
    pJet->nNeutrals  = itJet->neutralMultiplicity();
    pJet->nParticles = itJet->nConstituents ();

    if(triggerEvent      != 0) {pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, *triggerEvent); } 
    else                       {pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, *patTriggerObjects); }

    ////Add Extras
    baconhep::TAddJet *pAddJet = 0; 
    if(fComputeFullJetInfo) {
      assert(rExtraArray.GetEntries() < rExtraArray.GetSize());
      const int extraIndex = rExtraArray.GetEntries();
      new(rExtraArray[extraIndex]) baconhep::TAddJet();
      pAddJet = (baconhep::TAddJet*)rExtraArray[extraIndex];
      pAddJet->index = index;
      addJet(pAddJet, iEvent, *itJet, jetBaseRef);
    }
  } 
}

// === filler for MINIAOD ===
void FillerJet::fill(TClonesArray *array, TClonesArray *iExtraArray,
                     const edm::Event &iEvent, const edm::EventSetup &iSetup,
                     const reco::Vertex &pv,
                     const std::vector<TriggerRecord> &triggerRecords,
                     const pat::TriggerObjectStandAloneCollection &triggerObjects)
{
  assert(array);
  assert(!fComputeFullJetInfo || iExtraArray);
  fRand->SetSeed(iEvent.id().event());

  // Get jet collection
  edm::Handle<pat::JetCollection> hJetProduct;
  iEvent.getByToken(fTokPatJetName,hJetProduct);
  assert(hJetProduct.isValid());
  const pat::JetCollection *jetCol = hJetProduct.product();

  // Get vertex collection
  edm::Handle<reco::VertexCollection> hVertexProduct;
  iEvent.getByToken(fTokPVName,hVertexProduct);
  assert(hVertexProduct.isValid());
  const reco::VertexCollection *pvCol = hVertexProduct.product();
  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid()); 

  // Get Quark Gluon Likelihood
  edm::Handle<edm::ValueMap<float> > hQGLikelihood;
  edm::Handle<edm::ValueMap<float> > hQGLaxis2;
  edm::Handle<edm::ValueMap<float> > hQGLptD;
  edm::Handle<edm::ValueMap<int> >   hQGLmult;
  if(fQGLikelihood.size() > 0) { 
    iEvent.getByToken(fTokQGLikelihood,hQGLikelihood);
    assert(hQGLikelihood.isValid());

    iEvent.getByToken(fTokQGLAxis2,hQGLaxis2);
    assert(hQGLaxis2.isValid());
    
    iEvent.getByToken(fTokQGLPtD,hQGLptD);
    assert(hQGLptD.isValid());
    
    iEvent.getByToken(fTokQGLMult,hQGLmult);
    assert(hQGLmult.isValid());
  }
  // Get gen jet collection
  edm::Handle<reco::GenJetCollection> hGenJetProduct;
  const reco::GenJetCollection *genJetCol = 0;
  if(fUseGen) {
    iEvent.getByToken(fTokGenJetName,hGenJetProduct);
    assert(hGenJetProduct.isValid());
    genJetCol = hGenJetProduct.product();
  }
  TClonesArray &rArray      = *array;
  TClonesArray &rExtraArray = *iExtraArray;

  for(pat::JetCollection::const_iterator itJet = jetCol->begin(); itJet!=jetCol->end(); ++itJet) {

    edm::RefToBase<pat::Jet> jetBaseRef( edm::Ref<pat::JetCollection>(hJetProduct, itJet - jetCol->begin()) );

    double ptRaw = itJet->pt()*itJet->jecFactor("Uncorrected");

    // input to jet corrections
    double jetcorr = 1;
    if(fabs(itJet->eta()) < 5.191 && fApplyJEC) {
      fJetCorr->setJetPt(ptRaw);
      fJetCorr->setJetEta(itJet->eta());
      fJetCorr->setJetPhi(itJet->phi());
      fJetCorr->setJetE(itJet->energy());
      fJetCorr->setRho(*hRho);
      fJetCorr->setJetA(itJet->jetArea());
      fJetCorr->setJetEMF(-99.0);
      jetcorr = fJetCorr->getCorrection();
    }

    // jet pT cut (BOTH raw AND corrected pT must exceed threshold)
    if(ptRaw*jetcorr < fMinPt || ptRaw < fMinPt) continue;
    fJetUnc->setJetPt ( ptRaw  );
    fJetUnc->setJetEta( itJet->eta() );
    double jetunc = fJetUnc->getUncertainty(true);

    // construct object and place in array
    assert(rArray.GetEntries() < rArray.GetSize());
    const int index = rArray.GetEntries();
    new(rArray[index]) baconhep::TJet();
    baconhep::TJet *pJet = (baconhep::TJet*)rArray[index];

    //
    // Kinematics
    //==============================
    pJet->pt    = itJet->pt();
    if(fApplyJEC) pJet->pt    = ptRaw*jetcorr;
    pJet->eta   = itJet->eta();
    pJet->phi   = itJet->phi();
    pJet->mass  = itJet->mass();
    //if(fApplyJEC) pJet->mass    = itJet->mass()*jetcorr;
    pJet->ptRaw = ptRaw;
    pJet->area  = itJet->jetArea();
    pJet->unc   = jetunc;
    //
    // Impact Parameter
    //==============================
    pJet->d0 = JetTools::jetD0(*itJet, pv);
    pJet->dz = JetTools::jetDz(*itJet, pv);

    //
    // Identification
    //==============================
    pJet->csv = itJet->bDiscriminator(fCSVbtagName);
    if(fQGLikelihood.size() > 0) { 
      pJet->qgid  = (*hQGLikelihood)[jetBaseRef];
      pJet->axis2 = (*hQGLaxis2)[jetBaseRef];
      pJet->ptD   = (*hQGLptD)[jetBaseRef];
      pJet->mult  = (*hQGLmult)[jetBaseRef];
    }
    pJet->q = JetTools::jetCharge(*itJet);

    pJet->beta     = JetTools::beta(*itJet, pv);
    pJet->betaStar = JetTools::betaStar(*itJet, pv, pvCol);
    pJet->dR2Mean  = JetTools::dR2Mean(*itJet);

    pJet->mva = itJet->userFloat("pileupJetId:fullDiscriminant");

    TVector2 lPull    = JetTools::jetPull(*itJet,0);
    pJet->pullY       = lPull.X();
    pJet->pullPhi     = lPull.Y();
    TVector2 lChPull  = JetTools::jetPull(*itJet,1);
    pJet->chPullY     = lChPull.X();
    pJet->chPullPhi   = lChPull.Y();
    TVector2 lNeuPull = JetTools::jetPull(*itJet,2);
    pJet->neuPullY    = lNeuPull.X();
    pJet->neuPullPhi  = lNeuPull.Y();

    // Basic Noise Variables
    pJet->chEmFrac   = itJet->chargedEmEnergyFraction();
    pJet->neuEmFrac  = itJet->neutralEmEnergyFraction();
    pJet->chHadFrac  = itJet->chargedHadronEnergyFraction();
    pJet->neuHadFrac = itJet->neutralHadronEnergyFraction();
    pJet->muonFrac   = itJet->muonEnergyFraction();


    //
    // Generator matching
    //==============================
    const reco::GenJet * matchGenJet = 0;
    if(fUseGen) matchGenJet = match(&(*itJet),genJetCol);
    if(matchGenJet != 0) {
      pJet->partonFlavor = itJet->partonFlavour();
      pJet->hadronFlavor = itJet->hadronFlavour();
      pJet->genpt        = matchGenJet->pt();
      pJet->geneta       = matchGenJet->eta();
      pJet->genphi       = matchGenJet->phi();
      pJet->genm         = matchGenJet->mass();
    }

    pJet->nCharged   = itJet->chargedMultiplicity();
    pJet->nNeutrals  = itJet->neutralMultiplicity();
    pJet->nParticles = itJet->numberOfDaughters();

    pJet->hltMatchBits = TriggerTools::matchHLT(pJet->eta, pJet->phi, triggerRecords, triggerObjects);

    ////Add Extras
    baconhep::TAddJet *pAddJet = 0;
    if(fComputeFullJetInfo) {
      assert(rExtraArray.GetEntries() < rExtraArray.GetSize());
      const int extraIndex = rExtraArray.GetEntries();
      new(rExtraArray[extraIndex]) baconhep::TAddJet();
      pAddJet = (baconhep::TAddJet*)rExtraArray[extraIndex];
      pAddJet->index = index;
      addJet(pAddJet, iEvent, *itJet);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void FillerJet::addJet(baconhep::TAddJet *pAddJet, const edm::Event &iEvent, 
                       const reco::PFJet &itJet, const reco::JetBaseRef &jetBaseRef)
{ 
  // Get pruned jet collection
  edm::Handle<reco::BasicJetCollection> hPrunedJetProduct;
  iEvent.getByToken(fTokPrunedJetName,hPrunedJetProduct);
  assert(hPrunedJetProduct.isValid());
  const reco::BasicJetCollection *prunedJetCol = hPrunedJetProduct.product();

  // Get trimmed jet collection
  edm::Handle<reco::BasicJetCollection> hTrimmedJetProduct;
  iEvent.getByToken(fTokTrimmedJetName,hTrimmedJetProduct);
  assert(hTrimmedJetProduct.isValid());
  const reco::BasicJetCollection *trimmedJetCol = hTrimmedJetProduct.product();

  // Get soft drop jet collection
  edm::Handle<reco::BasicJetCollection> hSoftDropJetProduct;
  iEvent.getByToken(fTokSoftDropJetName,hSoftDropJetProduct);
  assert(hSoftDropJetProduct.isValid());
  const reco::BasicJetCollection *softdropJetCol = hSoftDropJetProduct.product();

  // Get sub-jet collection
  edm::Handle<reco::PFJetCollection> hSubJetProduct;
  iEvent.getByToken(fTokSubJets,hSubJetProduct);
  assert(hSubJetProduct.isValid());
  const reco::PFJetCollection *subJetCol = hSubJetProduct.product();

  // Get b sub-jets 
  edm::Handle<reco::JetTagCollection> hCSVbtagsSubJets;
  iEvent.getByToken(fTokCSVbtagSubJetName, hCSVbtagsSubJets);
  assert(hCSVbtagsSubJets.isValid());

  // Get double b tag
  edm::Handle<reco::JetTagCollection> hCSVDoubleBtag;
  iEvent.getByToken(fTokCSVDoubleBtagName, hCSVDoubleBtag);
  assert(hCSVDoubleBtag.isValid());

  //Get Quark Gluon Likelihood on subjets
  edm::Handle<edm::ValueMap<float> > hQGLikelihoodSubJets;
  iEvent.getByToken(fTokQGLSubJets,hQGLikelihoodSubJets);
  assert(hQGLikelihoodSubJets.isValid());
  
  // Get N-subjettiness moments
  edm::Handle<edm::ValueMap<float> > hTau1;
  iEvent.getByToken(fTokTau1Name,hTau1);                                                                                                                                                      
  assert(hTau1.isValid());
  edm::Handle<edm::ValueMap<float> > hTau2;
  iEvent.getByToken(fTokTau2Name,hTau2);
  assert(hTau2.isValid());
  edm::Handle<edm::ValueMap<float> > hTau3;
  iEvent.getByToken(fTokTau3Name,hTau3); 
  assert(hTau3.isValid());
  edm::Handle<edm::ValueMap<float> > hTau4;
  iEvent.getByToken(fTokTau4Name,hTau4); 
  assert(hTau4.isValid());

  // Get event energy density for jet correction
  edm::Handle<double> hRho;
  iEvent.getByToken(fTokRhoTag,hRho);
  assert(hRho.isValid());

  // Get Track and Secondary Vertex tag infos
  edm::Handle<reco::CandIPTagInfoCollection> hIPTagInfo;
  iEvent.getByToken(fTokIPTagInfos,hIPTagInfo);
  assert(hIPTagInfo.isValid());
  const reco::CandIPTagInfoCollection & ipTagInfoCollection = *(hIPTagInfo.product());

  edm::Handle<reco::CandSecondaryVertexTagInfoCollection> hSVTagInfo;
  iEvent.getByToken(fTokSVTagInfos, hSVTagInfo);
  assert(hSVTagInfo.isValid());
  const reco::CandSecondaryVertexTagInfoCollection & svTagInfoCollection = *(hSVTagInfo.product());

  pAddJet->pullAngle = JetTools::jetPullAngle(itJet,hSubJetProduct,fConeSize);
  pAddJet->tau1 = (*(hTau1.product()))[jetBaseRef];
  pAddJet->tau2 = (*(hTau2.product()))[jetBaseRef];
  pAddJet->tau3 = (*(hTau3.product()))[jetBaseRef];
  pAddJet->tau4 = (*(hTau4.product()))[jetBaseRef];
  pAddJet->doublecsv = (*(hCSVDoubleBtag.product()))[jetBaseRef];
  if(fShowerDeco != 0) { 
    std::vector<reco::CandidatePtr> pfConstituents = itJet.getJetConstituents();                                                                                                                     
    std::vector<fastjet::PseudoJet>   lClusterParticles;                                                                                                                                     
    for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {                                                                                                                                         
      reco::CandidatePtr pfcand = pfConstituents[ic];                                                                                                                                  
      fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());                                                                                                      
      lClusterParticles.push_back(pPart);                                                                                                                                                       
    }                                                                           
    pAddJet->topchi2 = fShowerDeco->chi(itJet.pt(),lClusterParticles);
    fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 2.0);
    fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
    std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
    fastjet::EnergyCorrelatorDoubleRatio C2beta0 (2,0. ,fastjet::EnergyCorrelator::pt_R);
    fastjet::EnergyCorrelatorDoubleRatio C2beta02(2,0.2,fastjet::EnergyCorrelator::pt_R);
    fastjet::EnergyCorrelatorDoubleRatio C2beta05(2,0.5,fastjet::EnergyCorrelator::pt_R);
    fastjet::EnergyCorrelatorDoubleRatio C2beta10(2,1.0,fastjet::EnergyCorrelator::pt_R);
    fastjet::EnergyCorrelatorDoubleRatio C2beta20(2,2.0,fastjet::EnergyCorrelator::pt_R);
    pAddJet->c2_0   = C2beta0 (inclusive_jets[0]);
    pAddJet->c2_0P2 = C2beta02(inclusive_jets[0]);
    pAddJet->c2_0P5 = C2beta05(inclusive_jets[0]);
    pAddJet->c2_1P0 = C2beta10(inclusive_jets[0]);
    pAddJet->c2_2P0 = C2beta20(inclusive_jets[0]);
  }

  double pCorr=1;
  const reco::BasicJet* matchJet = 0;

  // Pruning
  matchJet = match(&itJet,prunedJetCol);
  if(matchJet) {
    //pCorr = correction(pP1Jet,*hRho);
    pAddJet->mass_prun  = matchJet->mass()*pCorr;
  }

  // Trimming
  matchJet = match(&itJet,trimmedJetCol);
  if(matchJet) {
    //pCorr = correction(pP1Jet,*hRho);
    pAddJet->mass_trim  = matchJet->mass()*pCorr;
  }

  // Soft drop
  matchJet = match(&itJet,softdropJetCol);
  if(matchJet) {
    //pCorr = correction(pP1Jet,*hRho);
    pAddJet->mass_sd0  = matchJet->mass()*pCorr;
    pAddJet->pt_sd0 = matchJet->pt()*pCorr;
    pAddJet->eta_sd0 = matchJet->eta()*pCorr;
    pAddJet->phi_sd0 = matchJet->phi()*pCorr;
  }
  /*
  // Q-Jets
  pAddJet->qjet = 0;
  if(itJet.pt() > 100) pAddJet->qjet = JetTools::qJetVolatility(lClusterParticles,25.,fRand->Rndm());  // (!) why pT > 100 cut? computation time?
  //
  // Subjets
  //
  */
  // find/sort up to 4 hardest subjets
  const reco::PFJet *subjet1=0, *subjet2=0, *subjet3=0, *subjet4=0;
  double csv1=-2, csv2=-2, csv3=-2, csv4=-2;
  double qgid1=-2, qgid2=-2, qgid3=-2, qgid4=-2;
  double q1=-100, q2=-100, q3=-100, q4=-100;
  for(reco::PFJetCollection::const_iterator itSubJet = subJetCol->begin(); itSubJet!=subJetCol->end(); ++itSubJet) {
    if(reco::deltaR(itJet.eta(),itJet.phi(),itSubJet->eta(),itSubJet->phi())>fConeSize) continue;  // (!) get associated subjets by dR...is there a better way???

    reco::PFJetRef subjetRef(hSubJetProduct, itSubJet - subJetCol->begin());
    reco::JetBaseRef subjetBaseRef(subjetRef);

    if(!subjet1 || itSubJet->pt() > subjet1->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;
      
      subjet3 = subjet2;
      csv3    = csv2;
      qgid3   = qgid2;
      q3      = q2;
      
      subjet2 = subjet1;
      csv2    = csv1;
      qgid2   = qgid1;
      q2      = q1;
      
      subjet1 = &(*itSubJet);      
      csv1    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];      
      qgid1   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q1      = JetTools::jetCharge(*itSubJet);
      
    } else if(!subjet2 || itSubJet->pt() > subjet2->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;

      subjet3 = subjet2;
      csv3    = csv2;
      qgid3   = qgid2;
      q3      = q2;

      subjet2 = &(*itSubJet);      
      csv2    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef]; 
      qgid2   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
      q2      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet3 || itSubJet->pt() > subjet3->pt()) {
      subjet4 = subjet3;
      csv4    = csv3;
      qgid4   = qgid3;
      q4      = q3;

      subjet3 = &(*itSubJet);
      csv3    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      qgid3   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef]; 
      q3      = JetTools::jetCharge(*itSubJet);
      
    } else if(!subjet4 || itSubJet->pt() > subjet4->pt()) {
      subjet4 = &(*itSubJet);
      csv4    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      qgid4   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q4      = JetTools::jetCharge(*itSubJet);
    }
  }
  if(subjet1) {
    pAddJet->sj1_pt   = subjet1->pt();
    pAddJet->sj1_eta  = subjet1->eta();
    pAddJet->sj1_phi  = subjet1->phi();
    pAddJet->sj1_m    = subjet1->mass();
    pAddJet->sj1_csv  = csv1;
    pAddJet->sj1_qgid = qgid1;
    pAddJet->sj1_q    = q1;
  }
  if(subjet2) {
    pAddJet->sj2_pt   = subjet2->pt();
    pAddJet->sj2_eta  = subjet2->eta();
    pAddJet->sj2_phi  = subjet2->phi();
    pAddJet->sj2_m    = subjet2->mass();
    pAddJet->sj2_csv  = csv2;
    pAddJet->sj2_qgid = qgid2;
    pAddJet->sj2_q    = q2;
  }
  if(subjet3) {
    pAddJet->sj3_pt   = subjet3->pt();
    pAddJet->sj3_eta  = subjet3->eta();
    pAddJet->sj3_phi  = subjet3->phi();
    pAddJet->sj3_m    = subjet3->mass();
    pAddJet->sj3_csv  = csv3;
    pAddJet->sj3_qgid = qgid3;
    pAddJet->sj3_q    = q3;
  }
  if(subjet4) {
    pAddJet->sj4_pt   = subjet4->pt();
    pAddJet->sj4_eta  = subjet4->eta();
    pAddJet->sj4_phi  = subjet4->phi();
    pAddJet->sj4_m    = subjet4->mass();
    pAddJet->sj4_csv  = csv4;
    pAddJet->sj4_qgid = qgid4;
    pAddJet->sj4_q    = q4;
  }

  //
  // Top Tagging
  //
  pAddJet->topTagType=0;

  if(fTopTaggerName.compare("CMS")==0) {  // CMS Top Tagger

    edm::Handle<reco::BasicJetCollection> hCMSTTJetProduct;
    iEvent.getByToken(fTokCMSTTJetProduct,hCMSTTJetProduct);  // (!) hard-code
    assert(hCMSTTJetProduct.isValid());
    const reco::BasicJetCollection *cmsttJetCol = hCMSTTJetProduct.product();

    edm::Handle<reco::PFJetCollection> hCMSTTSubJetProduct;
    iEvent.getByToken(fTokCMSTTSubJetProduct,hCMSTTSubJetProduct);  // (!) hard-code
    assert(hCMSTTSubJetProduct.isValid());
    const reco::PFJetCollection *cmsttSubJetCol = hCMSTTSubJetProduct.product();

    matchJet = match(&itJet,cmsttJetCol);
    if(matchJet) {
      pAddJet->topTagType |= kCMSTT;

      unsigned int nsubjets=0;
      const reco::PFJet *sub1=0, *sub2=0, *sub3=0;
      for(reco::PFJetCollection::const_iterator itSub = cmsttSubJetCol->begin(); itSub!=cmsttSubJetCol->end(); ++itSub) {
        if(reco::deltaR(itJet.eta(),itJet.phi(),itSub->eta(),itSub->phi())>fConeSize) continue;  // (!) get associated subjets by dR...is there a better way???

        nsubjets++;

        if(!sub1 || itSub->pt() > sub1->pt()) {
          sub3 = sub2;
          sub2 = sub1;
          sub1 = &(*itSub);

        } else if(!sub2 || itSub->pt() > sub2->pt()) {
          sub3 = sub2;
          sub2 = &(*itSub);

        } else if(!sub3 || itSub->pt() > sub3->pt()) {
          sub3 = &(*itSub);
        }
      }
      pAddJet->top_n_subjets = nsubjets;

      TLorentzVector vSub1; if(sub1) { vSub1.SetPtEtaPhiM(sub1->pt(), sub1->eta(), sub1->phi(), sub1->mass()); }
      TLorentzVector vSub2; if(sub2) { vSub2.SetPtEtaPhiM(sub2->pt(), sub2->eta(), sub2->phi(), sub2->mass()); }
      TLorentzVector vSub3; if(sub3) { vSub3.SetPtEtaPhiM(sub3->pt(), sub3->eta(), sub3->phi(), sub3->mass()); }
      double m12 = (sub1 && sub2) ? (vSub1+vSub2).M() : 0;
      double m23 = (sub2 && sub3) ? (vSub2+vSub3).M() : 0;
      double m31 = (sub3 && sub1) ? (vSub3+vSub1).M() : 0;
      pAddJet->top_m_min=0;
      if     (m12 < m23 && m12 < m31) { pAddJet->top_m_min = m12; }
      else if(m23 < m12 && m23 < m31) { pAddJet->top_m_min = m23; }
      else if(m31 < m12 && m31 < m23) { pAddJet->top_m_min = m31; }

      pAddJet->top_m_123 = 0;
      pAddJet->top_fRec  = 0;
    }
  }
  /*
    } else if(fTopTaggerName.compare("HEP")==0) {  // HEP Top Tagger
    
    edm::Handle<reco::HTTTopJetTagInfoCollection> hHTTTagInfos;
    //iEvent.getByLabel("CA15HTTCHSAOD",hHTTTagInfos);  // (!) hard-code
    assert(hHTTTagInfos.isValid());
    const reco::HTTTopJetTagInfoCollection *httTagInfoCol = hHTTTagInfos.product();

    double dRmin = fConeSize;
    for(reco::HTTTopJetTagInfoCollection::const_iterator itTagInfo = httTagInfoCol->begin(); itTagInfo!=httTagInfoCol->end(); ++itTagInfo) {
      double dR = reco::deltaR(itJet.eta(), itJet.phi(),
      itTagInfo->properties().fjEta, itTagInfo->properties().fjPhi);
      
      if(dR < dRmin) {
        dRmin = dR;
        pAddJet->topTagType |= kHEPTT;
        pAddJet->top_n_subjets = 3;
        pAddJet->top_m_min     = 0;
        pAddJet->top_m_123     = itTagInfo->properties().topMass;
        pAddJet->top_fRec      = itTagInfo->properties().fRec;
      }
  }*/

  //
  // Single b tagger
  //
  for (std::size_t svIndex = 0; svIndex < svTagInfoCollection.size(); ++svIndex) {
    for (unsigned int itt=0; itt < ipTagInfoCollection.size(); ++itt) {
      const reco::JetBaseRef jetSv = svTagInfoCollection[svIndex].jet();
      const reco::VertexRef & vertexRef = ipTagInfoCollection[itt].primaryVertex();
      GlobalPoint pvtx(0.,0.,0.);
      if ( ipTagInfoCollection[itt].primaryVertex().isNonnull() )
	pvtx = GlobalPoint(vertexRef->x(),vertexRef->y(),vertexRef->z());

      // Re-calculate N-subjettiness using IVF vertices as composite b candidates - fill currentAxes
      std::vector<fastjet::PseudoJet> currentAxes;
      recalcNsubjettiness(jetSv,pAddJet->tau1,pAddJet->tau2,pAddJet->tau3,pAddJet->tau4,currentAxes);

      const std::vector<reco::CandidatePtr>  & selectedTracks( ipTagInfoCollection[itt].selectedTracks() );
      const std::vector<reco::btag::TrackIPData> & ipData = ipTagInfoCollection[itt].impactParameterData();

      reco::TrackKinematics allKinematics;

      float trackSip3dSig_3, trackSip3dSig_2, trackSip3dSig_1, trackSip3dSig_0;
      float tau1_trackSip3dSig_0, tau1_trackSip3dSig_1;
      float jet_NTracks = 0;

      std::vector<float> IP3Ds, IP3Ds_1, etaRels;
      int contTrk=0;

      for  (unsigned int iTrk=0; iTrk < selectedTracks.size(); ++iTrk) {
	const reco::CandidatePtr ptrackRef = selectedTracks[iTrk];
	const reco::Track & ptrack = *(reco::btag::toTrack(selectedTracks[iTrk]));

	float track_PVweight = 0.;
	setTracksPV(ptrackRef, vertexRef, track_PVweight);
	if (track_PVweight>0.) { allKinematics.add(ptrack, track_PVweight); }

	const reco::btag::TrackIPData &data = ipData[iTrk];
	bool isSelected = false;
	if (ftrackSelector(ptrack, data, itJet, pvtx)) isSelected = true;

	bool isfromV0 = false, isfromV0Tight = false;
	const reco::Track * trackPairV0Test[2];
	trackPairV0Test[0] = reco::btag::toTrack(selectedTracks[iTrk]);

	for (unsigned int jTrk=0; jTrk < selectedTracks.size(); ++jTrk){
	  if (iTrk == jTrk) continue;

	  const reco::btag::TrackIPData & pairTrackData = ipData[jTrk];
	  const reco::CandidatePtr pairTrackRef = selectedTracks[jTrk];
	  const reco::Track * pairTrackPtr = reco::btag::toTrack(pairTrackRef);
	  const reco::Track & pairTrack = *pairTrackPtr;

	  trackPairV0Test[1] = pairTrackPtr;

	  if (!ftrackPairV0Filter(trackPairV0Test, 2))
	    {
	      isfromV0 = true;
	      if (ftrackSelector(pairTrack, pairTrackData, itJet, pvtx) )
		isfromV0Tight = true;
	    }

	  if (isfromV0 && isfromV0Tight)
	    break;
	}

	if( isSelected && !isfromV0Tight ) jet_NTracks += 1.;

	reco::TransientTrack transientTrack = theTTBuilder->build(ptrack);
	GlobalVector direction(itJet.px(), itJet.py(), itJet.pz());

	// Find closest tau axis to track and calculate direction of that axis                                                                                                                                              
	double smallest = 100;
	fastjet::PseudoJet closestTauAxisToTrack;
	for (size_t i=0; i<currentAxes.size(); ++i){
	  if (reco::deltaR2(ptrack,currentAxes[i]) < smallest)
	    {
	      smallest = reco::deltaR2(ptrack,currentAxes[i]);
	      closestTauAxisToTrack = currentAxes[i];
	    }
	  else
	    continue;
	}
	direction = GlobalVector(closestTauAxisToTrack.px(), closestTauAxisToTrack.py(), closestTauAxisToTrack.pz());

	float Track_decayLengthTau = -1;
	float Track_distTauAxis    = -1;
	float Track_IPsig          = -1;

	TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), *vertexRef, direction, transientTrack.field());
	if (closest.isValid())
	  Track_decayLengthTau =  (closest.globalPosition() - RecoVertex::convertPos(vertexRef->position())).mag();
	Track_distTauAxis = std::abs(IPTools::jetTrackDistance(transientTrack, direction, *vertexRef).second.value());
	Track_IPsig = ipTagInfoCollection[itt].impactParameterData()[iTrk].ip3d.significance();

	if(!isfromV0 && Track_decayLengthTau<5. && Track_distTauAxis <0.07) {
	  IP3Ds.push_back( Track_IPsig <-50. ? -50. : Track_IPsig);
	  contTrk++;
	  // store 3DIPsig for the track with the closest tau axis
	  for (size_t i=0; i<currentAxes.size(); ++i){
	    if (reco::deltaR2(ptrack,currentAxes[i]) == reco::deltaR2(ptrack,closestTauAxisToTrack)) {
	      IP3Ds_1.push_back(Track_IPsig <-50. ? -50. : Track_IPsig);
	    }
	    else
	      continue;
	  }
	}
      } // end loop on selectedTracks
      pAddJet->nTracks = jet_NTracks;

      // charm and bottom cuts
      std::vector<size_t> indices = ipTagInfoCollection[itt].sortedIndexes(reco::btag::IP2DSig);
      bool charmThreshSet = false;

      pAddJet->trackSip2dSigAboveCharm = -19;
      pAddJet->trackSip2dSigAboveBottom = -19;

      reco::TrackKinematics kin;
      for (size_t i =0; i<indices.size(); ++i) {
	size_t idx = indices[i];
	const reco::btag::TrackIPData & data = ipData[idx];
	const reco::Track & track = (*reco::btag::toTrack(selectedTracks[idx]));
	kin.add(track);
	if ( kin.vectorSum().M() > 1.5 && !charmThreshSet ) {
	  pAddJet->trackSip2dSigAboveCharm  = data.ip2d.significance();
	  charmThreshSet = true;
	}
	if ( kin.vectorSum().M() > 5.2 ) {
	  pAddJet->trackSip2dSigAboveBottom  = data.ip2d.significance();
	  break;
	}
      }

      float dummyTrack = -50.;
      std::sort( IP3Ds.begin(),IP3Ds.end(),std::greater<float>() );
      std::sort( IP3Ds_1.begin(),IP3Ds_1.end(),std::greater<float>() );
      int num_1 = IP3Ds_1.size();
      switch(contTrk){
      case 0:
	trackSip3dSig_0 = dummyTrack;
	trackSip3dSig_1 = dummyTrack;
	trackSip3dSig_2 = dummyTrack;
	trackSip3dSig_3 = dummyTrack;
	break;
      case 1:
	trackSip3dSig_0 = IP3Ds.at(0);
	trackSip3dSig_1 = dummyTrack;
	trackSip3dSig_2 = dummyTrack;
	trackSip3dSig_3 = dummyTrack;
	break;
      case 2:
	trackSip3dSig_0 = IP3Ds.at(0);
	trackSip3dSig_1 = IP3Ds.at(1);
	trackSip3dSig_2 = dummyTrack;
	trackSip3dSig_3 = dummyTrack;
	break;
      case 3:
	trackSip3dSig_0 = IP3Ds.at(0);
	trackSip3dSig_1 = IP3Ds.at(1);
	trackSip3dSig_2 = IP3Ds.at(2);
	trackSip3dSig_3 = dummyTrack;
	break;
      default:
	trackSip3dSig_0 = IP3Ds.at(0);
	trackSip3dSig_1 = IP3Ds.at(1);
	trackSip3dSig_2 = IP3Ds.at(2);
	trackSip3dSig_3 = IP3Ds.at(3);
      }
      switch(num_1){
      case 0:
	tau1_trackSip3dSig_0 = dummyTrack;
	tau1_trackSip3dSig_1 = dummyTrack;
	break;
      case 1:
	tau1_trackSip3dSig_0 = IP3Ds_1.at(0);
	tau1_trackSip3dSig_1 = dummyTrack;
	break;
      default:
	tau1_trackSip3dSig_0 = IP3Ds_1.at(0);
	tau1_trackSip3dSig_1 = IP3Ds_1.at(1);
      }
      pAddJet->trackSip3dSig_3  = trackSip3dSig_3;
      pAddJet->trackSip3dSig_2  = trackSip3dSig_2;
      pAddJet->trackSip3dSig_1  = trackSip3dSig_1;
      pAddJet->trackSip3dSig_0  = trackSip3dSig_0;

      pAddJet->tau1_trackSip3dSig_0  = tau1_trackSip3dSig_0;
      pAddJet->tau1_trackSip3dSig_1  = tau1_trackSip3dSig_1;

      math::XYZTLorentzVector allSum = allKinematics.weightedVectorSum() ;

      // Secondary Vertices                                                                                                                                                                                                   
      int contSV = 0, vertexNTracks = 0;
      math::XYZVector jetDir = jetSv->momentum().Unit();
      std::map<double, size_t> VTXmass;
      std::map<double, size_t> VTXfd;

      // Fill VTXfd and VTXmass                                                                                                                                                                                               
      for (size_t vtx = 0; vtx < (size_t)svTagInfoCollection[svIndex].nVertices(); ++vtx) {
	vertexNTracks += (svTagInfoCollection[svIndex].secondaryVertex(vtx)).numberOfSourceCandidatePtrs();
	GlobalVector flightDir = svTagInfoCollection[svIndex].flightDirection(vtx);
	if (reco::deltaR2(flightDir, jetDir)<(fmaxSVDeltaRToJet*fmaxSVDeltaRToJet)) {
	  ++contSV;
	  VTXmass[svTagInfoCollection[svIndex].secondaryVertex(vtx).p4().mass()]=vtx;
	  VTXfd[svTagInfoCollection[svIndex].flightDistance(vtx).error()]=vtx;
	}
      }
      pAddJet->nSV = VTXfd.size();
      
      //Find leading SV in mass (higher mass) and flight distance error (lowest error) and tau axis closest to them                                                                                                           
      if (VTXmass.size()>0) {
	std::vector<tauAxes> tausSVmass; //to order closestTauAxistoSVMass and closestTauAxistoSV                                                                                                                           
	const reco::VertexCompositePtrCandidate &vertex_SVmass = svTagInfoCollection[svIndex].secondaryVertex(VTXmass.rbegin()->second);
	GlobalVector flightDir_SVmass = svTagInfoCollection[svIndex].flightDirection(VTXmass.rbegin()->second);
	
	for (size_t i=0; i<3; ++i) {
	  tausSVmass.push_back(tauAxes(currentAxes[i],reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(VTXmass.rbegin()->second),currentAxes[i])));
	}
	std::sort(tausSVmass.begin(),tausSVmass.end(),compAxesSV);
	
	//loop over the vertices that belong to the tau_SVmass                                                                                                                                                              
	pAddJet->tau_SVmass_nSecondaryVertices = 0;
	pAddJet->tau_SVmass_flightDistance2dSig = -1;
	pAddJet->tau_SVmass_vertexDeltaR = -1;
	pAddJet->tau_SVmass_vertexNTracks = 0;
	pAddJet->tau_SVmass_vertexEnergyRatio = -1;
	pAddJet->tau_SVmass_vertexMass = -1;
	pAddJet->tau_SVmass_vertexMass_corrected = -1;
	pAddJet->tau_SVmass_zratio = -5;
	float tau_SVmass_nSecondaryVertices = 0;
	
	reco::TrackKinematics tau_SVmass_Kinematics;
	std::vector<float> tau_SVmass_trackEtaRels;
	for ( std::map<double, size_t>::reverse_iterator iVtx=VTXmass.rbegin(); iVtx!=VTXmass.rend(); ++iVtx) {
	  const reco::VertexCompositePtrCandidate &vertex =  svTagInfoCollection[svIndex].secondaryVertex(iVtx->second);
	  reco::TrackKinematics vtxKinematics;
	  vertexKinematicsAndChange(vertex, vtxKinematics);
	  if ( (reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVmass[0].tau) < reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVmass[1].tau))
	       &&  (reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVmass[0].tau) < reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVmass[2].tau))
	       ) {
	    tau_SVmass_Kinematics  = tau_SVmass_Kinematics + vtxKinematics;
	    pAddJet->tau_SVmass_vertexNTracks  += svTagInfoCollection[svIndex].secondaryVertex(iVtx->second).numberOfSourceCandidatePtrs();
	    if( pAddJet->tau_SVmass_flightDistance2dSig  < 0 ) {
	      pAddJet->tau_SVmass_flightDistance2dSig  = svTagInfoCollection[svIndex].flightDistance(iVtx->second,true).significance();
	      pAddJet->tau_SVmass_vertexDeltaR  = reco::deltaR(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVmass[0].tau);
	    }
	    
	    etaRelToTauAxis(vertex, tausSVmass[0].tau, tau_SVmass_trackEtaRels);
	    tau_SVmass_nSecondaryVertices += 1;
	  }
	}
	pAddJet->tau_SVmass_nSecondaryVertices = tau_SVmass_nSecondaryVertices;
	float tau_SVmass_trackEtaRel_0, tau_SVmass_trackEtaRel_1, tau_SVmass_trackEtaRel_2;
	std::sort( tau_SVmass_trackEtaRels.begin(),tau_SVmass_trackEtaRels.end() );
	float tau_SVmass_dummyEtaRel = -1.;
	switch(tau_SVmass_trackEtaRels.size()) {
	case 0:
	  tau_SVmass_trackEtaRel_0 = tau_SVmass_dummyEtaRel;
	  tau_SVmass_trackEtaRel_1 = tau_SVmass_dummyEtaRel;
	  tau_SVmass_trackEtaRel_2 = tau_SVmass_dummyEtaRel;
	  break;
	case 1:
	  tau_SVmass_trackEtaRel_0 = tau_SVmass_trackEtaRels.at(0);
	  tau_SVmass_trackEtaRel_1 = tau_SVmass_dummyEtaRel;
	  tau_SVmass_trackEtaRel_2 = tau_SVmass_dummyEtaRel;
	  break;
	case 2:
	  tau_SVmass_trackEtaRel_0 = tau_SVmass_trackEtaRels.at(0);
	  tau_SVmass_trackEtaRel_1 = tau_SVmass_trackEtaRels.at(1);
	  tau_SVmass_trackEtaRel_2 = tau_SVmass_dummyEtaRel;
	  break;
	default:
	  tau_SVmass_trackEtaRel_0 = tau_SVmass_trackEtaRels.at(0);
	  tau_SVmass_trackEtaRel_1 = tau_SVmass_trackEtaRels.at(1);
	  tau_SVmass_trackEtaRel_2 = tau_SVmass_trackEtaRels.at(2);
	}
	pAddJet->tau_SVmass_trackEtaRel_2  = tau_SVmass_trackEtaRel_2;
	pAddJet->tau_SVmass_trackEtaRel_1  = tau_SVmass_trackEtaRel_1;
	pAddJet->tau_SVmass_trackEtaRel_0  = tau_SVmass_trackEtaRel_0;

	if ( pAddJet->tau_SVmass_nSecondaryVertices  > 0. ) {
	  // vertexEnergyRatio                                                                                                                                                                                            
	  math::XYZTLorentzVector tau_SVmass_vertexSum = tau_SVmass_Kinematics.weightedVectorSum();
	  pAddJet->tau_SVmass_vertexEnergyRatio  = tau_SVmass_vertexSum.E() / allSum.E();
	  // vertexMass                                                                                                                                                                                                   
	  pAddJet->tau_SVmass_vertexMass  = tau_SVmass_vertexSum.M();
	  double tau_SVmass_vertexPt2 = math::XYZVector(flightDir_SVmass.x(), flightDir_SVmass.y(), flightDir_SVmass.z()).Cross(tau_SVmass_vertexSum).Mag2() / flightDir_SVmass.mag2();
	  pAddJet->tau_SVmass_vertexMass_corrected  =std::sqrt( tau_SVmass_vertexSum.M() *  tau_SVmass_vertexSum.M() + tau_SVmass_vertexPt2) + std::sqrt(tau_SVmass_vertexPt2);
	  // z ratio                                                                                                                                                                                                      
	  pAddJet->tau_SVmass_zratio =  reco::deltaR2(tausSVmass[1].tau+tausSVmass[2].tau,tausSVmass[0].tau)*(vertex_SVmass.p4().pt()/vertex_SVmass.p4().mass());
	}
      }

      if (VTXfd.size()>0) {
	std::vector<tauAxes> tausSVfd;
	const reco::VertexCompositePtrCandidate &vertex_SVfd = svTagInfoCollection[svIndex].secondaryVertex(VTXfd.begin()->second);
	GlobalVector flightDir_SVfd = svTagInfoCollection[svIndex].flightDirection(VTXfd.begin()->second);
	
	for (size_t i=0; i<3; ++i) {
	  tausSVfd.push_back(tauAxes(currentAxes[i],reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(VTXfd.begin()->second),currentAxes[i])));
	}
	// order taus according to dR to SV with lowest flighDistanceError                                                                                                                                                  
	std::sort(tausSVfd.begin(),tausSVfd.end(),compAxesSV);
	
	pAddJet->tau_SVfd_nSecondaryVertices = 0;
	pAddJet->tau_SVfd_flightDistance2dSig = -1;
	pAddJet->tau_SVfd_vertexDeltaR = -1;
	pAddJet->tau_SVfd_vertexNTracks = 0;
	pAddJet->tau_SVfd_vertexEnergyRatio = -1;
	pAddJet->tau_SVfd_vertexMass = -1;
	pAddJet->tau_SVfd_vertexMass_corrected = -1;
	pAddJet->tau_SVfd_zratio = -5;
	float tau_SVfd_nSecondaryVertices = 0;
	
	reco::TrackKinematics tau_SVfd_Kinematics;
	std::vector<float> tau_SVfd_trackEtaRels;
	for ( std::map<double, size_t>::iterator iVtx=VTXfd.begin(); iVtx!=VTXfd.end(); ++iVtx) {
	  const reco::VertexCompositePtrCandidate &vertex =  svTagInfoCollection[svIndex].secondaryVertex(iVtx->second);
	  reco::TrackKinematics vtxKinematics;
	  vertexKinematicsAndChange(vertex, vtxKinematics);
	  if ( (reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVfd[0].tau) < reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVfd[1].tau) )
	       &&  (reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVfd[0].tau) < reco::deltaR2(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVfd[2].tau))
	       ) {
	    tau_SVfd_Kinematics  = tau_SVfd_Kinematics + vtxKinematics;
	    pAddJet->tau_SVfd_vertexNTracks  += svTagInfoCollection[svIndex].secondaryVertex(iVtx->second).numberOfSourceCandidatePtrs();
	    if( pAddJet->tau_SVfd_flightDistance2dSig  < 0 ) {
	      pAddJet->tau_SVfd_flightDistance2dSig  = svTagInfoCollection[svIndex].flightDistance(iVtx->second,true).significance();
	      pAddJet->tau_SVfd_vertexDeltaR  = reco::deltaR(svTagInfoCollection[svIndex].flightDirection(iVtx->second),tausSVfd[0].tau);
	    }
	   
	    etaRelToTauAxis(vertex, tausSVfd[0].tau, tau_SVfd_trackEtaRels);
	    tau_SVfd_nSecondaryVertices += 1;
	  }
	}
	pAddJet->tau_SVfd_nSecondaryVertices = tau_SVfd_nSecondaryVertices;
	float tau_SVfd_trackEtaRel_0, tau_SVfd_trackEtaRel_1, tau_SVfd_trackEtaRel_2;
	std::sort( tau_SVfd_trackEtaRels.begin(),tau_SVfd_trackEtaRels.end() );
	float tau_SVfd_dummyEtaRel = -1.;
	switch(tau_SVfd_trackEtaRels.size()){
	case 0:
	  tau_SVfd_trackEtaRel_0 = tau_SVfd_dummyEtaRel;
	  tau_SVfd_trackEtaRel_1 = tau_SVfd_dummyEtaRel;
	  tau_SVfd_trackEtaRel_2 = tau_SVfd_dummyEtaRel;
	  break;
	case 1:
	  tau_SVfd_trackEtaRel_0 = tau_SVfd_trackEtaRels.at(0);
	  tau_SVfd_trackEtaRel_1 = tau_SVfd_dummyEtaRel;
	  tau_SVfd_trackEtaRel_2 = tau_SVfd_dummyEtaRel;
	  break;
	case 2:
	  tau_SVfd_trackEtaRel_0 = tau_SVfd_trackEtaRels.at(0);
	  tau_SVfd_trackEtaRel_1 = tau_SVfd_trackEtaRels.at(1);
	  tau_SVfd_trackEtaRel_2 = tau_SVfd_dummyEtaRel;
	  break;
	default:
	  tau_SVfd_trackEtaRel_0 = tau_SVfd_trackEtaRels.at(0);
	  tau_SVfd_trackEtaRel_1 = tau_SVfd_trackEtaRels.at(1);
	  tau_SVfd_trackEtaRel_2 = tau_SVfd_trackEtaRels.at(2);
	}
	pAddJet->tau_SVfd_trackEtaRel_2  = tau_SVfd_trackEtaRel_2;
	pAddJet->tau_SVfd_trackEtaRel_1  = tau_SVfd_trackEtaRel_1;
	pAddJet->tau_SVfd_trackEtaRel_0  = tau_SVfd_trackEtaRel_0;
	
	if ( pAddJet->tau_SVfd_nSecondaryVertices  > 0. ) {
	  math::XYZTLorentzVector tau_SVfd_vertexSum = tau_SVfd_Kinematics.weightedVectorSum();
	  pAddJet->tau_SVfd_vertexEnergyRatio  = tau_SVfd_vertexSum.E() / allSum.E();
	  pAddJet->tau_SVfd_vertexMass  = tau_SVfd_vertexSum.M();
	  double tau_SVfd_vertexPt2 = math::XYZVector(flightDir_SVfd.x(), flightDir_SVfd.y(), flightDir_SVfd.z()).Cross(tau_SVfd_vertexSum).Mag2() / flightDir_SVfd.mag2();
	  pAddJet->tau_SVfd_vertexMass_corrected  =std::sqrt(tau_SVfd_vertexSum.M() * tau_SVfd_vertexSum.M() + tau_SVfd_vertexPt2) + std::sqrt(tau_SVfd_vertexPt2);
	  pAddJet->tau_SVfd_zratio =  reco::deltaR2(tausSVfd[1].tau+tausSVfd[2].tau,tausSVfd[0].tau)*(vertex_SVfd.p4().pt()/vertex_SVfd.p4().mass());
	}
      }
    }
  }
}

void FillerJet::addJet(baconhep::TAddJet *pAddJet, const edm::Event &iEvent, const pat::Jet &itJet)
{
  //pAddJet->pullAngle = JetTools::jetPullAngle(itJet,hSubJetProduct,fConeSize);
  pAddJet->tau1 = itJet.userFloat(fJettinessName + std::string(":tau1"));
  pAddJet->tau2 = itJet.userFloat(fJettinessName + std::string(":tau2"));
  pAddJet->tau3 = itJet.userFloat(fJettinessName + std::string(":tau3"));
  pAddJet->tau4 = -1;  //(!) Not in MINIAOD

  // Pruning
  pAddJet->mass_prun = itJet.userFloat(fPrunedJetName + std::string("Mass"));

  // Trimming
  pAddJet->mass_trim = itJet.userFloat(fTrimmedJetName + std::string("Mass"));

  // Soft drop
  pAddJet->mass_sd0 = itJet.userFloat(fSoftDropJetName + std::string("Mass"));

  // Jet Shape Correlation observables
//  fastjet::JetDefinition lCJet_def(fastjet::cambridge_algorithm, 2.0);
//  fastjet::ClusterSequence lCClust_seq(lClusterParticles, lCJet_def);
//  std::vector<fastjet::PseudoJet> inclusive_jets = lCClust_seq.inclusive_jets(0);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta0 (2,0. ,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta02(2,0.2,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta05(2,0.5,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta10(2,1.0,fastjet::EnergyCorrelator::pt_R);
//  fastjet::EnergyCorrelatorDoubleRatio C2beta20(2,2.0,fastjet::EnergyCorrelator::pt_R);
//  pAddJet->c2_0   = C2beta0 (inclusive_jets[0]);
//  pAddJet->c2_0P2 = C2beta02(inclusive_jets[0]);
//  pAddJet->c2_0P5 = C2beta05(inclusive_jets[0]);
//  pAddJet->c2_1P0 = C2beta10(inclusive_jets[0]);
//  pAddJet->c2_2P0 = C2beta20(inclusive_jets[0]);

  // Q-Jets
//  pAddJet->qjet = 0;
//  if(itJet->pt() > 100) pAddJet->qjet = JetTools::qJetVolatility(lClusterParticles,25.,fRand->Rndm());  // (!) why pT > 100 cut? computation time?

  //
  // Subjets
  //
/*
  // find/sort up to 4 hardest subjets
  const pat::Jet *subjet1=0, *subjet2=0, *subjet3=0, *subjet4=0;
  double csv1=-2, csv2=-2, csv3=-2, csv4=-2;
  double qgid1=-2, qgid2=-2, qgid3=-2, qgid4=-2;
  double q1=-100, q2=-100, q3=-100, q4=-100;
  pat::JetPtrCollection const &subjetCol = itJet.subjets(fSubJetName);
  for(edm::Ptr<pat::Jet> const &itSubJet : subjetCol) {

    if(!subjet1 || itSubJet->pt() > subjet1->pt()) {
      subjet4 = subjet3;
      //csv4    = csv3;
      //qgid4   = qgid3;
      q4      = q3;

      subjet3 = subjet2;
      //csv3    = csv2;
      //qgid3   = qgid2;
      q3      = q2;

      subjet2 = subjet1;
      //csv2    = csv1;
      //qgid2   = qgid1;
      q2      = q1;

      subjet1 = itSubJet.get();
      //csv1    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid1   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q1      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet2 || itSubJet->pt() > subjet2->pt()) {
      subjet4 = subjet3;
      //csv4    = csv3;
      //qgid4   = qgid3;
      q4      = q3;

      subjet3 = subjet2;
      //csv3    = csv2;
      //qgid3   = qgid2;
      q3      = q2;

      subjet2 = itSubJet.get();
      //csv2    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid2   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q2      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet3 || itSubJet->pt() > subjet3->pt()) {
      subjet4 = subjet3;
      //csv4    = csv3;
      //qgid4   = qgid3;
      q4      = q3;

      subjet3 = itSubJet.get();
      //csv3    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid3   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q3      = JetTools::jetCharge(*itSubJet);

    } else if(!subjet4 || itSubJet->pt() > subjet4->pt()) {
      subjet4 = itSubJet.get();
      //csv4    = (*(hCSVbtagsSubJets.product()))[subjetBaseRef];
      //qgid4   = (*(hQGLikelihoodSubJets.product()))[subjetBaseRef];
      q4      = JetTools::jetCharge(*itSubJet);
    }
  }
  if(subjet1) {
    pAddJet->sj1_pt   = subjet1->pt();
    pAddJet->sj1_eta  = subjet1->eta();
    pAddJet->sj1_phi  = subjet1->phi();
    pAddJet->sj1_m    = subjet1->mass();
    pAddJet->sj1_csv  = csv1;
    pAddJet->sj1_qgid = qgid1;
    pAddJet->sj1_q    = q1;
  }
  if(subjet2) {
    pAddJet->sj2_pt   = subjet2->pt();
    pAddJet->sj2_eta  = subjet2->eta();
    pAddJet->sj2_phi  = subjet2->phi();
    pAddJet->sj2_m    = subjet2->mass();
    pAddJet->sj2_csv  = csv2;
    pAddJet->sj2_qgid = qgid2;
    pAddJet->sj2_q    = q2;
  }
  if(subjet3) {
    pAddJet->sj3_pt   = subjet3->pt();
    pAddJet->sj3_eta  = subjet3->eta();
    pAddJet->sj3_phi  = subjet3->phi();
    pAddJet->sj3_m    = subjet3->mass();
    pAddJet->sj3_csv  = csv3;
    pAddJet->sj3_qgid = qgid3;
    pAddJet->sj3_q    = q3;
  }
  if(subjet4) {
    pAddJet->sj4_pt   = subjet4->pt();
    pAddJet->sj4_eta  = subjet4->eta();
    pAddJet->sj4_phi  = subjet4->phi();
    pAddJet->sj4_m    = subjet4->mass();
    pAddJet->sj4_csv  = csv4;
    pAddJet->sj4_qgid = qgid4;
    pAddJet->sj4_q    = q4;
  }
*/
  //
  // Top Tagging
  //
  if(fTopTaggerName.compare("CMS")==0) {  // CMS Top Tagger
    reco::CATopJetTagInfo const *tagInfo = dynamic_cast<reco::CATopJetTagInfo const *>(itJet.tagInfo("caTop"));  // (!) hard-code
    if(tagInfo) {
      pAddJet->topTagType |= kCMSTT;
      pAddJet->top_n_subjets = tagInfo->properties().nSubJets;
      pAddJet->top_m_min     = tagInfo->properties().minMass;
      pAddJet->top_m_123     = 0;
      pAddJet->top_fRec      = 0;
    }
/*  } else if(fTopTaggerName.compare("HEP")==0) {  // HEP Top Tagger
    reco::HTTTopJetTagInfo const *tagInfo = dynamic_cast<reco::HTTTopJetTagInfo const *>(itJet.tagInfo("CA15HTTCHSMINIAOD"));  // (!) hard-code
    if(tagInfo) {
      pAddJet->topTagType |= kHEPTT;
      pAddJet->top_n_subjets = 3;
      pAddJet->top_m_min     = 0;
      pAddJet->top_m_123     = tagInfo->properties().topMass;
      pAddJet->top_fRec      = tagInfo->properties().fRec;
    }*/
  }
}

//--------------------------------------------------------------------------------------------------
const reco::BasicJet* FillerJet::match( const reco::PFJet *iJet,const reco::BasicJetCollection *jets ) { 
  int lId = -1;
  double dRmin = 999.;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::BasicJet* jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if(dR > fConeSize) continue;
    if ( dR < dRmin ) {
      dRmin = dR;
      lId = i;
    }
  }
  const reco::BasicJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

const reco::BasicJet* FillerJet::match(const pat::Jet *iJet, const reco::BasicJetCollection *jets) {
  int lId = -1;
  double dRmin = 999.;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::BasicJet* jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if(dR > fConeSize) continue;
    if ( dR < dRmin ) {
      dRmin = dR;
      lId = i;
    }
  }
  const reco::BasicJet* lJet = 0;
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

const reco::GenJet* FillerJet::match( const reco::PFJet *iJet,const reco::GenJetCollection *jets ) { 
  int lId = -1;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::GenJet *jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if ( dR < 0.25 ) {
      lId = i;
      break;
    }
  }
  const reco::GenJet* lJet = 0; 
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

const reco::GenJet* FillerJet::match(const pat::Jet *iJet, const reco::GenJetCollection *jets) {
  int lId = -1;
  for ( unsigned int i=0; i< jets->size(); ++i ) {
    const reco::GenJet *jet = &(*jets)[i];
    float dR = deltaR( iJet->eta(), iJet->phi(), jet->eta(), jet->phi() );
    if ( dR < 0.25 ) {
      lId = i;
      break;
    }
  }
  const reco::GenJet* lJet = 0;
  if(lId != -1) lJet = &((*jets)[lId]);
  return lJet;
}

void FillerJet::recalcNsubjettiness(const reco::JetBaseRef &jet, float & tau1, float & tau2, float & tau3, float & tau4, std::vector<fastjet::PseudoJet> & currentAxes)
{
  std::vector<fastjet::PseudoJet> fjParticles;

  for(const reco::CandidatePtr & daughter : jet->daughterPtrVector())
    {
      if ( daughter.isNonnull() && daughter.isAvailable() )
        fjParticles.push_back( fastjet::PseudoJet( daughter->px(), daughter->py(), daughter->pz(), daughter->energy() ) );
      else
	edm::LogWarning("MissingJetConstituent") << "Jet constituent required for N-subjettiness computation is missing!";
    }
  // calculate N-subjetiness                                                                                                                                                                                                      
  tau1 = fnjettiness.getTau(1, fjParticles);
  tau2 = fnjettiness.getTau(2, fjParticles);
  tau3 = fnjettiness.getTau(3, fjParticles);
  tau4 = fnjettiness.getTau(4, fjParticles);
  currentAxes = fnjettiness.currentAxes();
}

void FillerJet::setTracksPVBase(const reco::TrackRef & trackRef, const reco::VertexRef & vertexRef, float & PVweight) const
{
  PVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::Vertex::trackRef_iterator IT;

  const reco::Vertex & vtx = *(vertexRef);
  for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
    {
      const reco::TrackBaseRef & baseRef = *it;
      if( baseRef == trackBaseRef )
        {
          PVweight = vtx.trackWeight(baseRef);
          break;
	}
    }
}

void FillerJet::setTracksPV(const reco::CandidatePtr & trackRef, const reco::VertexRef & vertexRef, float & PVweight) const
{
  PVweight = 0.;
  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(trackRef.get());
  if(pcand) // MiniAOD case                                                                                                                                                                                                       
    {
      if( pcand->fromPV() == pat::PackedCandidate::PVUsedInFit )
        {
          PVweight = 1.;
        }
    }
  else
    {
      const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(trackRef.get());
      setTracksPVBase(pfcand->trackRef(), vertexRef, PVweight);
    }
}

void FillerJet::vertexKinematicsAndChange(const reco::VertexCompositePtrCandidate & vertex, reco::TrackKinematics & vertexKinematics)
{
  const std::vector<reco::CandidatePtr> & tracks = vertex.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track) {
    const reco::Track& mytrack = *(*track)->bestTrack();
    vertexKinematics.add(mytrack, 1.0);
  }
}

void FillerJet::etaRelToTauAxis(const reco::VertexCompositePtrCandidate & vertex, const fastjet::PseudoJet & tauAxis, std::vector<float> & tau_trackEtaRel)
{
  math::XYZVector direction(tauAxis.px(), tauAxis.py(), tauAxis.pz());
  const std::vector<reco::CandidatePtr> & tracks = vertex.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track)
    tau_trackEtaRel.push_back(std::abs(reco::btau::etaRel(direction.Unit(), (*track)->momentum())));
}
