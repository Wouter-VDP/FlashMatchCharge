////////////////////////////////////////////////////////////////////////
// Class:       PandoraAnalyzer
// Module Type: analyzer
// File:        PandoraAnalyzer_module.cc
//
// Generated at Thu Jun 23 00:24:52 2016 by Lorena Escudero Sanchez using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// TODO
// - Put fidvol, electron_energy_threshold and proton_energy_threshold as fcl parameters
// - Use Geometry service for TPC size

#include <fstream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//uncomment the lines below as you use these objects

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "PandoraAnalysis/PandoraAnalysis.hh"

using namespace lar_pandora;

namespace geo { class Geometry; }

namespace test {
  class PandoraAnalyzer;
}

double x_start = 0;
double x_end = 256.35;
double y_start = -116.5;
double y_end = 116.5;
double z_start = 0;
double z_end = 1036.8;
double fidvol = 10;

bool is_fiducial(double x[3], double d) {
  bool is_x = x[0] > (x_start+d) && x[0] < (x_end-d);
  bool is_y = x[1] > (y_start+d) && x[1] < (y_end-d);
  bool is_z = x[2] > (z_start+d) && x[2] < (z_end-d);
  return is_x && is_y && is_z;
}

double distance(double a[3], double b[3]) {
  double d = 0;

  for (int i = 0; i < 3; i++) {
    d += pow((a[i]-b[i]),2);
  }

  return sqrt(d);
}

class test::PandoraAnalyzer : public art::EDAnalyzer {
public:
  explicit PandoraAnalyzer(fhicl::ParameterSet const & pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  virtual ~PandoraAnalyzer();

  // Plugins should not be copied or assigned.
  PandoraAnalyzer(PandoraAnalyzer const &) = delete;
  PandoraAnalyzer(PandoraAnalyzer &&) = delete;
  PandoraAnalyzer & operator = (PandoraAnalyzer const &) = delete;
  PandoraAnalyzer & operator = (PandoraAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  void reconfigure(fhicl::ParameterSet const &pset) override;

private:

  test::PandoraAnalysis fMyAnalysisObj;
  TFile * myTFile;
  TTree * myTTree;
  TEfficiency * e_energy;
  bool         m_printDebug;

};


test::PandoraAnalyzer::PandoraAnalyzer(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset)  // ,
 // More initializers here.
{

  //create output tree
  art::ServiceHandle<art::TFileService> tfs;
  myTFile = new TFile("PandoraAnalyzerOutput.root", "RECREATE");
  myTTree = tfs->make<TTree>("pandoratree","PandoraAnalysis Tree");

  e_energy = tfs->make<TEfficiency>("e_energy",";#nu_{e} energy [GeV];",30,0,3);


  //add branches


  this->reconfigure(pset);

}

test::PandoraAnalyzer::~PandoraAnalyzer()
{

  //store output tree
  myTFile->cd();
  myTTree->Write("pandoratree");
  myTFile->Close();

  std::cout << "End!" << std::endl;
}



void test::PandoraAnalyzer::analyze(art::Event const & evt)
{

  //do the analysis
  art::InputTag pandoraNu_tag { "pandoraNu" };
  art::InputTag generator_tag { "generator" };


  auto const& generator_handle = evt.getValidHandle< std::vector< simb::MCTruth > >( generator_tag );
  auto const& generator(*generator_handle);

  std::vector<simb::MCParticle> nu_mcparticles;
  for (int i = 0; i < generator[0].NParticles(); i++) {
    if (generator[0].Origin() == 1) {
      nu_mcparticles.push_back(generator[0].GetParticle(i));
    }
  }

  double proton_energy = 0;
  double electron_energy = 0;
  int protons = 0;
  bool is_pion = false;
  bool is_electron = false;
  double proton_energy_threshold = 0.06;
  double electron_energy_threshold = 0.035;

  for (auto& mcparticle: nu_mcparticles) {
    if (mcparticle.Process() == "primary" and mcparticle.T() != 0 and mcparticle.StatusCode() == 1) {

      double position[3] = {mcparticle.Position().X(),mcparticle.Position().Y(),mcparticle.Position().Z()};
      double end_position[3] = {mcparticle.EndPosition().X(),mcparticle.EndPosition().Y(),mcparticle.EndPosition().Z()};

      if (mcparticle.PdgCode() == 2212 && (mcparticle.E()-mcparticle.Mass()) > proton_energy_threshold && is_fiducial(position,fidvol) && is_fiducial(end_position,fidvol)) {
        protons++;
        proton_energy = std::max(mcparticle.E()-mcparticle.Mass(), proton_energy);
      }

      if (mcparticle.PdgCode() == 11 && (mcparticle.E()-mcparticle.Mass()) > electron_energy_threshold && is_fiducial(position,fidvol)) {
        is_electron = true;
        electron_energy = std::max(mcparticle.E()-mcparticle.Mass(), electron_energy);
      }

      if (mcparticle.PdgCode() == 211 || mcparticle.PdgCode() == 111) {
        is_pion = true;
      }

    }
  }

  double nu_energy = generator[0].GetNeutrino().Nu().E();

  if (is_electron && !is_pion && protons == 1 && nu_energy > 0.2) {
    std::cout << "CCQE 1e1p event" << std::endl;

    bool reco_shower = false;
    bool reco_track = false;
    int showers = 0;
    int tracks = 0;

    try {
      auto const& pfparticle_handle = evt.getValidHandle< std::vector< recob::PFParticle > >( pandoraNu_tag );
      auto const& pfparticles(*pfparticle_handle);

      art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, evt, pandoraNu_tag);


      for (size_t ipf = 0; ipf < pfparticles.size(); ipf++) {

        bool is_neutrino = (abs(pfparticles[ipf].PdgCode()) == 12 || abs(pfparticles[ipf].PdgCode()) == 14) && pfparticles[ipf].IsPrimary();

        // Is a nu_e or nu_mu PFParticle?
        if (!is_neutrino) continue;

        double closest_distance = std::numeric_limits<double>::max();

        double neutrino_vertex[3];

        auto const& neutrino_vertex_obj = vertex_per_pfpart.at(ipf);
        neutrino_vertex_obj->XYZ(neutrino_vertex); // PFParticle neutrino vertex coordinates

        // Is the vertex within fiducial volume?
        if (!is_fiducial(neutrino_vertex, fidvol)) continue;

        closest_distance = std::min(distance(neutrino_vertex,neutrino_vertex),closest_distance);
        //cout << pfparticles[ipf].PdgCode() << " " << distance(neutrino_vertex,correct_neutrino_vertex) << endl;

        // Loop over the neutrino daughters and check if there is a shower and a track
        for (auto const& pfdaughter: pfparticles[ipf].Daughters()) {

          // auto const& daughter_vertex_obj = vertex_per_pfpart.at(pfdaughter);
          // double daughter_vertex[3];
          // daughter_vertex_obj->XYZ(daughter_vertex);
          // double distance_daugther = distance(neutrino_vertex,daughter_vertex);

          if (pfparticles[pfdaughter].PdgCode() == 11) {
            reco_shower = true;
            showers++;
          }

          if (pfparticles[pfdaughter].PdgCode() == 13) {
            reco_track = true;
            tracks++;
          }

        } // end for pfparticle daughters


      } // end for pfparticles


    } catch (...) {
      std::cout << "NO RECO DATA PRODUCTS" << std::endl;
    }

    e_energy->Fill(showers == 1 && tracks == 1, nu_energy);

  } // end CCQE if


} // end analyze function

//------------------------------------------------------------------------------------------------------------------------------------


void test::PandoraAnalyzer::reconfigure(fhicl::ParameterSet const & pset)
{

  //TODO: add an external fcl file to change configuration
  //add what you want to read, and default values of your labels etc. example:
  //  m_particleLabel = pset.get<std::string>("PFParticleModule","pandoraNu");

  m_printDebug = pset.get<bool>("PrintDebug",false);

}

//---------------------------------------------------------------------------------------------------------------------------
//add other functions here

DEFINE_ART_MODULE(test::PandoraAnalyzer)
