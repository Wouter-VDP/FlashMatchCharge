////////////////////////////////////////////////////////////////////////
// Class:       FlashChargeAna
// Plugin Type: analyzer (art v2_06_03)
// File:        FlashChargeAna_module.cc
//
// by Wouter Van De Pontseele 
////////////////////////////////////////////////////////////////////////


#ifndef ELECTRON_EVENT_SELECTION_ALG_H
#define ELECTRON_EVENT_SELECTION_ALG_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "SpaceChargeMicroBooNE.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "TTree.h"
#include "TVector3.h"

#include <algorithm>

namespace flashcharge {
  class FlashChargeAna;
}
class flashcharge::FlashChargeAna : public art::EDAnalyzer {

public:

  explicit FlashChargeAna(fhicl::ParameterSet const & pset);

  void reconfigure(fhicl::ParameterSet const & pset);
  // Plugins should not be copied or assigned.
  FlashChargeAna(FlashChargeAna const &) = delete;
  FlashChargeAna(FlashChargeAna &&) = delete;
  FlashChargeAna & operator = (FlashChargeAna const &) = delete;
  FlashChargeAna & operator = (FlashChargeAna &&) = delete;
  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // functions
  void fillTree(art::Event const & e);
  void fillTruthTree(art::Event const & e);
  void fillPandoraTree(art::Event const & e);
  void fillOticalTree(art::Event const & e);
  bool is_fiducial(const std::vector<double> & x) const;
  TVector3 calculateChargeCenter(
                        size_t top_particle_index,
                        const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                        const art::Event & evt);

  double trackEnergy(const art::Ptr<recob::Track>& track, const art::Event & e);
  void traversePFParticleTree(
                        const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                        size_t top_index,
                        std::vector<size_t> & unordered_daugthers );

  // variables
  SpaceChargeMicroBooNE SCE = SpaceChargeMicroBooNE("SCEoffsets_MicroBooNE_E273.root");
  TTree*      fTree;
  

  bool bool_truth   = true;
  bool bool_pandora = true;
  bool bool_optical = true;

  /* FCL VARIABLES */

  double m_fidvolXstart;
  double m_fidvolXend;

  double m_fidvolYstart;
  double m_fidvolYend;

  double m_fidvolZstart;
  double m_fidvolZend;


  double m_fractionsigmaflashwidth;
  double m_absoluteflashdist;

  int    m_nrPMT;
  double m_chargetolight;


  /* TREE VARIABLES*/

  //Run Subrun Event
  Short_t    run;
  Short_t    subrun;
  Int_t    event;

  //Truth information
  Bool_t    true_fid;                                  ///< Is there a neutrino in the fiducial volume? 
  Short_t   mcevts_truth;                              ///< number of neutrino Interactions in the spill
  std::vector<Short_t>   nuPDG_truth;                  ///< neutrino PDG code
  std::vector<Short_t>   ccnc_truth;                   ///< 0=CC 1=NC
  std::vector<Short_t>   mode_truth;                   ///< 0=QE/El, 1=RES, 2=DIS, 3=Coherent production
  std::vector<double>   enu_truth;                     ///< true neutrino energy
  std::vector<double>   nuvtxx_truth;                  ///< neutrino vertex x
  std::vector<double>   nuvtxy_truth;                  ///< neutrino vertex y
  std::vector<double>   nuvtxz_truth;                  ///< neutrino vertex z
  std::vector<double>   nuvtxx_sc;                     ///< neutrino vertex x
  std::vector<double>   nuvtxy_sc;                     ///< neutrino vertex y
  std::vector<double>   nuvtxz_sc;                     ///< neutrino vertex z

  //PandoraNu information
  Short_t             nnuvtx;                          ///< Number of PandoraNu neutrino candidate vertices
  std::vector<double> nuvtxx;                          ///< x coordinate
  std::vector<double> nuvtxy;                          ///< y coordinate
  std::vector<double> nuvtxz;                          ///< z coordinate
  std::vector<Short_t> nuvtxpdg;                        ///< PDG code assigned by PandoraNu
  std::vector<double> center_of_charge_x;              ///< x Center of deposited charge
  std::vector<double> center_of_charge_y;              ///< y Center of deposited charge
  std::vector<double> center_of_charge_z;              ///< z Center of deposited charge

  std::vector<std::vector< TVector3 >> shwr_dir;        ///< The direction of shower for every shower connected to the passed pandoraNu candidates
  std::vector<std::vector< double >> shwr_en;          ///< Shower energy
  std::vector<std::vector< double >> shwr_angle;       ///< Shower opening angle

  std::vector<std::vector< TVector3 >> trck_dir;        ///< The start direction of track for every track connected to the passed pandoraNu candidates
  std::vector<std::vector< double >> trck_len;         ///< Length of the track
  std::vector<std::vector< double >> trck_dedxavg;     ///< Average dedx of the track


  //Optical information
  Short_t nfls;                                         ///< Number of reconstructed flashes
  std::vector<double> flsTime;                         ///< Flash time (us)
  std::vector<double> flsPe;                           ///< Flash total PE
  std::vector<double> flsPePMT;                        ///< Flash PE per PMT, needs a reshape, nrPMT
  std::vector<double> flsYcenter;                      ///< Flash Y center (cm)
  std::vector<double> flsZcenter;                      ///< Flash Z center (cm)
  std::vector<double> flsYwidth;                       ///< Flash Y width (cm)
  std::vector<double> flsZwidth;                       ///< Flash Z width (cm)
};




void flashcharge::FlashChargeAna::reconfigure(fhicl::ParameterSet const & pset)
{
  // Implementation of optional member function here.
  m_fidvolXstart             = pset.get<double>("fidvolXstart", 10);
  m_fidvolXend               = pset.get<double>("fidvolXstart", 10);

  m_fidvolYstart             = pset.get<double>("fidvolYstart", 20);
  m_fidvolYend               = pset.get<double>("fidvolYend", 20);

  m_fidvolZstart             = pset.get<double>("fidvolZstart", 10);
  m_fidvolZend               = pset.get<double>("fidvolZend", 50);

  m_fractionsigmaflashwidth  = pset.get<double>("fractionsigmaflashwidth", 2.0);
  m_absoluteflashdist        = pset.get<double>("absoluteflashdist", 50.0);

  m_chargetolight            = pset.get<double>("chargetolight", 10000.0);
  m_nrPMT                    = pset.get<int>   ("nrPMT", 32);
}

flashcharge::FlashChargeAna::FlashChargeAna(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{
  //initialize output tree
  std::cout << "Initializing output tree..." << std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fTree  = tfs->make<TTree>("flashtree","FlashAnalysis Tree");

  //Set branches for (Run Subrun Event)
  fTree->Branch("run",     &run,     "run/S"       );
  fTree->Branch("subrun",  &subrun,  "subrun/S"    );
  fTree->Branch("event",   &event,   "event/I"     );

  //Set branches for truth information
  fTree->Branch("true_fid",     &true_fid,                 "true_fid/O"            );
  fTree->Branch("mcevts_truth", &mcevts_truth,             "mcevts_truth/S"        );
  fTree->Branch("nuPDG_truth",  "std::vector<Short_t>",    &nuPDG_truth            );
  fTree->Branch("ccnc_truth",   "std::vector<Short_t>",    &ccnc_truth             );
  fTree->Branch("mode_truth",   "std::vector<Short_t>",    &mode_truth             );
  fTree->Branch("enu_truth",    "std::vector<double>",    &enu_truth               );
  fTree->Branch("nuvtxx_truth", "std::vector<double>",    &nuvtxx_truth            );
  fTree->Branch("nuvtxy_truth", "std::vector<double>",    &nuvtxy_truth            );
  fTree->Branch("nuvtxz_truth", "std::vector<double>",    &nuvtxz_truth            );
  fTree->Branch("nuvtxx_sc",    "std::vector<double>",    &nuvtxx_sc               );
  fTree->Branch("nuvtxy_sc",    "std::vector<double>",    &nuvtxy_sc               );
  fTree->Branch("nuvtxz_sc",    "std::vector<double>",    &nuvtxz_sc               );

  //Set branches for PandoraNU information
  fTree->Branch("nnuvtx",                   &nnuvtx,                  "nnuvtx/S"             );
  fTree->Branch("nuvtxx",                   "std::vector<double>",    &nuvtxx                );
  fTree->Branch("nuvtxy",                   "std::vector<double>",    &nuvtxy                );
  fTree->Branch("nuvtxz",                   "std::vector<double>",    &nuvtxz                );
  fTree->Branch("nuvtxpdg",                 "std::vector<Short_t>",   &nuvtxpdg              );
  fTree->Branch("center_of_charge_x",       "std::vector<double>",    &center_of_charge_x    );
  fTree->Branch("center_of_charge_y",       "std::vector<double>",    &center_of_charge_y    );
  fTree->Branch("center_of_charge_z",       "std::vector<double>",    &center_of_charge_z    );

  fTree->Branch("shwr_dir",   "std::vector<std::vector<TVector3>>",   &shwr_dir              );
  fTree->Branch("shwr_en",    "std::vector<std::vector<double>>",     &shwr_en               );
  fTree->Branch("shwr_angle", "std::vector<std::vector<double>>",     &shwr_angle            );

  fTree->Branch("trck_dir",   "std::vector<std::vector<TVector3>>",   &trck_dir              );
  fTree->Branch("trck_len",   "std::vector<std::vector<double>>",     &shwr_en               );
  fTree->Branch("trck_dedxavg","std::vector<std::vector<double>>",    &shwr_angle            );


  //Set branches for optical information
  fTree->Branch("nfls",       &nfls,       "nfls/S"                );
  fTree->Branch("flsTime",    "std::vector<double>",   &flsTime    );
  fTree->Branch("flsPe",      "std::vector<double>",   &flsPe      );
  fTree->Branch("flsPePMT",   "std::vector<double>",   &flsPePMT   );
  fTree->Branch("flsYcenter", "std::vector<double>",   &flsYcenter );
  fTree->Branch("flsZcenter", "std::vector<double>",   &flsZcenter );
  fTree->Branch("flsYwidth",  "std::vector<double>",   &flsYwidth  );
  fTree->Branch("flsZwidth",  "std::vector<double>",   &flsZwidth  );

  this->reconfigure(pset);
}



#endif // ELECTRON_EVENT_SELECTION_ALG_H