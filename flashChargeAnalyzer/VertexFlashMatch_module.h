////////////////////////////////////////////////////////////////////////
// Class:       FitSimPhotons
// File:        FitSimPhotons_module.cc
// Date:        May 4, 2017
// Author:      Wouter Van De Pontseele
////////////////////////////////////////////////////////////////////////

#ifndef VERTEX_FLASH_MATCH_H
#define VERTEX_FLASH_MATCH_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "fhiclcpp/ParameterSet.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightCharge.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"

#include "Spacecharge/SpaceChargeMicroBooNE.h"

#include "TTree.h"
#include "TVector3.h"
#include <numeric>

class VertexFlashMatch;

class VertexFlashMatch : public art::EDAnalyzer {
public:
  explicit VertexFlashMatch(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  ~VertexFlashMatch();

  // Plugins should not be copied or assigned.
  VertexFlashMatch(VertexFlashMatch const &) = delete;
  VertexFlashMatch(VertexFlashMatch &&) = delete;
  VertexFlashMatch &operator=(VertexFlashMatch const &) = delete;
  VertexFlashMatch &operator=(VertexFlashMatch &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  /* LEVEL 0 FUCNTIONS */
  // L0: Clears all the variables to store in the TTree
  void resetTreeVar();

  // L0: Main function will call other functions to actually do the filling
  void fillTree(art::Event const &e);

  /* LEVEL 1 FUCNTIONS */

  // L1: Fills the true information returns.
  void fillTrueTree(art::Event const &e);

  // L1: Fills the pandora information returns a vector of QClusters that can be
  // used in flashmatching
  std::vector<flashana::QCluster_t> fillPandoraTree(art::Event const &e);

  // L1: Fills the optical tree and returns the flash that we want to use for
  // matching
  flashana::Flash_t fillOticalTree(art::Event const &e);

  // L1: Matches a flash to a cluster, given a list of clusters, the
  // hypothetical flash is also included in the returned result.
  std::vector<flashana::FlashMatch_t>
  makeMatch(std::vector<flashana::QCluster_t> const &cluster,
            flashana::Flash_t const &flashReco);

  // L1: Store the match information to
  void fillMatchTree(std::vector<flashana::FlashMatch_t> const &matchvec);

  /* LEVEL 2 FUCNTIONS */
  // L2: Calculates the center of charge given a list of pfpindices (ony using
  // collection plane currently)
  std::vector<double> calculateChargeCenter(art::Event const &e,
                                            std::vector<size_t> const &pfplist);

  // L2: Returns a QCluster given a list of pfparticles (ony using collection
  // plane currently)
  flashana::QCluster_t collect3DHits(art::Event const &e,
                                     std::vector<size_t> const &pfplist);

  /* LEVEL 3 FUCNTIONS */
  // L3: Returns an unordered list of the daughters of a PFParticle, recursively
  // and not including itself
  void traversePFParticleTree(
      size_t top_index, std::vector<size_t> &unordered_daugthers,
      const art::ValidHandle<std::vector<recob::PFParticle>> pfparticles);

  // Variables
  ::flashana::FlashMatchManager m_mgr;
  art::ServiceHandle<geo::Geometry> m_geo;
  TTree *m_tree;
  SpaceChargeMicroBooNE *m_sce;

  /* TREE VARIABLES*/
  // Run Subrun Event
  unsigned short run;
  unsigned short subrun;
  unsigned int event;

  // MC information
  unsigned short
      true_pdg; ///< The true particle pdgcode, for single particle generation
  double
      true_energy; ///< The true particle energy, for single particle generation
  double true_time; ///< The true particle interaction time, for single particle
                    ///generation
  double true_x;    ///< The true particle interaction positon
  double true_y;
  double true_z;
  double true_px; ///< The true particle start momentum (direction)
  double true_py;
  double true_pz;
  double true_sce_x; ///< True vertex position, corrected for spacecharge
  double true_sce_y;
  double true_sce_z;
  std::vector<int> simphot_spectrum; ///< Array with PMT of the simphotons

  // PandoraNu information (vector like fields with a star have nr_pfp+1 size,
  // the plus one is when the event is considered as a whole.)
  unsigned short nr_pfp; ///< Number of pfparticles in the event

  unsigned short nr_nupfp; ///< Number of PandoraNu neutrino candidate vertices
  std::vector<double> nuvtxx;     ///< x coordinate
  std::vector<double> nuvtxy;     ///< y coordinate
  std::vector<double> nuvtxz;     ///< z coordinate
  std::vector<Short_t> nupfp_pdg; ///< PDG code assigned by PandoraNu

  std::vector<double> q_Y_sps; ///< The charge on the collection plane from
                               ///spacepoints associated to pandora hierachy
  std::vector<double> center_of_charge_x; ///< x Center of deposited charge
  std::vector<double> center_of_charge_y; ///< y Center of deposited charge
  std::vector<double> center_of_charge_z; ///< z Center of deposited charge
  std::vector<double> min_x_sps;          ///< minimal spacepoint x value

  std::vector<Short_t>
      tpc_id; ///< the xth neutrino candidate corresponds to this match.
  std::vector<double> width_of_flash_x;  ///< x width of opFlash
  std::vector<double> center_of_flash_x; ///< x Center of opFlash
  std::vector<double> matchscore; ///< Matchscore of the single opflash with the
                                  ///qcluster containing all the pfparticles
  std::vector<std::vector<double>> hypo_spectrum; ///< Array with PMT of the
                                                  ///hypothetical flash made by
                                                  ///the flashmatcher

  // Reconstructed photon information
  short nr_flash;                 ///< Number of flashes in window
  std::vector<int> reco_spectrum; ///< Array with PMT of the flash
  double recphot_time;            ///< time of the simplebeamflash
  double center_of_flash_y;       ///< y Center of opFlash
  double center_of_flash_z;       ///< z center of opFlash
  double width_of_flash_y;        ///< y width of opFlash
  double width_of_flash_z;        ///< z width of opFlash

  /* FCL VARIABLES */
  double m_startbeamtime;
  double m_endbeamtime;

  double m_ly_proton; // coeficient of charge light for particles
  double m_ly_electron;
  double m_ly_muon;
  double m_ly_gamma;
  std::map<unsigned short, double> m_ly_map;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

VertexFlashMatch::~VertexFlashMatch(){
    if (m_sce){
        delete m_sce;
    }
}

VertexFlashMatch::VertexFlashMatch(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) {
  // Prepare fixed length arrays
  simphot_spectrum.resize(m_geo->NOpDets(), 0);
  reco_spectrum.resize(m_geo->NOpDets(), 0);

  // initialize fcl parameters
  m_startbeamtime = p.get<double>("FlashVetoTimeStart", 3.2);
  m_endbeamtime = p.get<double>("FlashVetoTimeEnd", 4.8);

  m_ly_proton = p.get<double>("LightYieldProton", 1.0);
  m_ly_electron = p.get<double>("LightYieldElectron", 1.0);
  m_ly_muon = p.get<double>("LightYieldMuon", 1.0);
  m_ly_gamma = p.get<double>("LightYieldGamma", 1.0);

  m_ly_map = {{2212, m_ly_proton},
              {11, m_ly_electron},
              {13, m_ly_muon},
              {22, m_ly_gamma}};

  m_mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));

  // Initialize sce:
  std::string _env = std::getenv("UBOONE_DATA_DIR");
  _env = _env + "/SpaceCharge/SCEoffsets_MicroBooNE_E273.root";
  m_sce = new SpaceChargeMicroBooNE(_env);

  // initialize output tree
  art::ServiceHandle<art::TFileService> tfs;
  m_tree = tfs->make<TTree>("VertexFlashMatch", "VertexFlashMatch Tree");

  // Set branches for (Run Subrun Event)
  m_tree->Branch("run", &run, "run/s");
  m_tree->Branch("subrun", &subrun, "subrun/s");
  m_tree->Branch("event", &event, "event/i");

  // Set branches for MC information
  m_tree->Branch("true_pdg", &true_pdg, "true_pdg/s");
  m_tree->Branch("true_energy", &true_energy, "true_energy/D");
  m_tree->Branch("true_time", &true_time, "true_time/D");
  m_tree->Branch("true_x", &true_x, "true_x/D");
  m_tree->Branch("true_y", &true_y, "true_y/D");
  m_tree->Branch("true_z", &true_z, "true_z/D");
  m_tree->Branch("true_px", &true_px, "true_px/D");
  m_tree->Branch("true_py", &true_py, "true_py/D");
  m_tree->Branch("true_pz", &true_pz, "true_pz/D");
  m_tree->Branch("true_sce_x", &true_sce_x, "true_sce_x/D");
  m_tree->Branch("true_sce_y", &true_sce_y, "true_sce_y/D");
  m_tree->Branch("true_sce_z", &true_sce_z, "true_sce_z/D");
  m_tree->Branch("simphot_spectrum", "std::vector<int>", &simphot_spectrum);

  // Reconstructed photon information
  m_tree->Branch("nr_flash", &nr_flash, "nr_flash/s");
  m_tree->Branch("reco_spectrum", "std::vector<int>", &reco_spectrum);
  m_tree->Branch("recphot_time", &recphot_time, "recphot_time/D");
  m_tree->Branch("center_of_flash_y", &center_of_flash_y,
                 "center_of_flash_y/D");
  m_tree->Branch("center_of_flash_z", &center_of_flash_z,
                 "center_of_flash_z/D");
  m_tree->Branch("width_of_flash_y", &width_of_flash_y, "width_of_flash_y/D");
  m_tree->Branch("width_of_flash_z", &width_of_flash_z, "width_of_flash_z/D");
  m_tree->Branch("tpc_id", "std::vector<Short_t>", &tpc_id);

  // Set branches for PandoraNU information
  m_tree->Branch("nr_pfp", &nr_pfp, "nr_pfp/s");

  m_tree->Branch("q_Y_sps", "std::vector<double>", &q_Y_sps);

  m_tree->Branch("nr_nupfp", &nr_nupfp, "nr_nupfp/s");
  m_tree->Branch("nuvtxx", "std::vector<double>", &nuvtxx);
  m_tree->Branch("nuvtxy", "std::vector<double>", &nuvtxy);
  m_tree->Branch("nuvtxz", "std::vector<double>", &nuvtxz);
  m_tree->Branch("nupfp_pdg", "std::vector<Short_t>", &nupfp_pdg);

  m_tree->Branch("center_of_charge_x", "std::vector<double>",
                 &center_of_charge_x);
  m_tree->Branch("center_of_charge_y", "std::vector<double>",
                 &center_of_charge_y);
  m_tree->Branch("center_of_charge_z", "std::vector<double>",
                 &center_of_charge_z);
  m_tree->Branch("min_x_sps", "std::vector<double>", &min_x_sps);

  m_tree->Branch("width_of_flash_x", "std::vector<double>", &width_of_flash_x);
  m_tree->Branch("center_of_flash_x", "std::vector<double>",
                 &center_of_flash_x);
  m_tree->Branch("matchscore", "std::vector<double>", &matchscore);
  m_tree->Branch("hypo_spectrum", "std::vector<std::vector<double>>",
                 &hypo_spectrum);
}

void VertexFlashMatch::resetTreeVar() {
  /* TREE VARIABLES*/
  // Run Subrun Event
  run = 0;
  subrun = 0;
  event = 0;

  // MC information
  true_pdg = 0;
  true_energy = 0;
  true_time = 0;
  true_x = 0;
  true_y = 0;
  true_z = 0;
  true_px = 0;
  true_py = 0;
  true_pz = 0;
  true_sce_x = 0;
  true_sce_y = 0;
  true_sce_z = 0;

  // PandoraNu information
  nr_pfp = 0;
  nr_nupfp = 0;
  nr_flash = 0;
  recphot_time = 0;
  center_of_flash_y = 0;
  center_of_flash_z = 0;
  width_of_flash_y = 0;
  width_of_flash_z = 0;

  // Arrays of fixed length
  std::fill(simphot_spectrum.begin(), simphot_spectrum.end(), 0);
  std::fill(reco_spectrum.begin(), reco_spectrum.end(), 0);

  // Arrays with variable length
  hypo_spectrum.clear(); // double array, contains arrays of fixed length
  q_Y_sps.clear();
  nuvtxx.clear();
  nuvtxy.clear();
  nuvtxz.clear();
  nupfp_pdg.clear();
  center_of_charge_x.clear();
  center_of_charge_y.clear();
  center_of_charge_z.clear();
  min_x_sps.clear();

  width_of_flash_x.clear();
  center_of_flash_x.clear();
  tpc_id.clear();
  matchscore.clear();
}
#endif // VERTEX_FLASH_MATCH_H