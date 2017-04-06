////////////////////////////////////////////////////////////////////////
// Class:       PandoraAnalyzer
// Module Type: analyzer
// File:        PandoraAnalyzer_module.cc
//
// Generated at Thu Jun 23 00:24:52 2016 by Lorena Escudero Sanchez using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////
// Holaaaa
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

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TTree.h"
#include "TFile.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "PandoraAnalysis/PandoraAnalysis.hh"

using namespace lar_pandora;

namespace geo { class Geometry; }

namespace test {
  class PandoraAnalyzer;
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

  auto const& generator_handle = evt.getValidHandle< std::vector< simb::MCTruth > >( "generator" );
  auto const& generator(*generator_handle);
  
  std::vector<simb::MCParticle> nu_mcparticles;
  for (int i = 0; i < generator[0].NParticles(); i++) {
    if (generator[0].Origin() == 1) {
      nu_mcparticles.push_back(generator[0].GetParticle(i));
    }
  }

}

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
