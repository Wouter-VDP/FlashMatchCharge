////////////////////////////////////////////////////////////////////////
// Class:       FitSimPhotons
// File:        FitSimPhotons_module.cc
// Date:        May 4, 2017
// Author:      Wouter Van De Pontseele
////////////////////////////////////////////////////////////////////////

#ifndef FIT_SIM_PHOTONS_H
#define FIT_SIM_PHOTONS_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "fhiclcpp/ParameterSet.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/Simulation/SimPhotons.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightCharge.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"

#include "TTree.h"
#include "TVector3.h"
#include <numeric>


class FitSimPhotons;

class FitSimPhotons : public art::EDAnalyzer
{
public:

    explicit FitSimPhotons(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    FitSimPhotons(FitSimPhotons const &) = delete;
    FitSimPhotons(FitSimPhotons &&) = delete;
    FitSimPhotons & operator = (FitSimPhotons const &) = delete;
    FitSimPhotons & operator = (FitSimPhotons &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

private:

    //Functions
    void clearTreeVar();
    void fillTree             (art::Event const & e);
    void fillPandoraTree      (art::Event const & e);
    flashana::Flash_t fillOticalTree       (art::Event const & e);
    void calculateChargeCenter(art::Event const & e);

    flashana::QCluster_t collect3DHitsZ( 	size_t pfindex,
                                            const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                                            art::Event const & e);

    flashana::Flash_t makeMatch(art::Event const & e, flashana::Flash_t & flashReco);
    void makeMatch(flashana::Flash_t const & flash, flashana::QCluster_t const & );




    //Variables
    std::vector<flashana::FlashMatch_t> m_result;
    ::flashana::FlashMatchManager       m_mgr;
    art::ServiceHandle<geo::Geometry>   m_geo;
    TTree*                              m_tree;

    /* TREE VARIABLES*/
    //Run Subrun Event
    unsigned short            run;
    unsigned short            subrun;
    unsigned int              event;

    //PandoraNu information
    float                     q_z_total;              ///< The total charge on the collectionplane in beamwindow
    float                     q_z_sps;                ///< The charge on the collectionplane in beamwindow from PFPtree spacepoints
    float                     q_z_hit;                ///< The charge on the collectionplane in beamwindow from PFPtree hits
    std::vector<double>       flashhypo_channel;      ///< Array with PMT channel of the hypothetical flash made by the flashmatcher (DOUBLE)
    float                     flashhypo_time;         ///< time of the hypothetical flash made by the flashmatcher
    float                     center_of_charge_x;     ///< x Center of deposited charge
    float                     center_of_charge_y;     ///< y Center of deposited charge
    float                     center_of_charge_z;     ///< z Center of deposited charge

    //MC information
   	unsigned short            true_pdg;				  ///< The true particle pdgcode, for single particle generation
    float                     true_energy;            ///< The true particle energy, for single particle generation
    float                     true_time;              ///< The true particle interaction time, for single particle generation
    float                     true_x;                 ///< The true particle interaction positon
    float                     true_z;
    float                     true_y;
    std::vector<float>        simphot_time;           ///< Array with times of the simphotons
    std::vector<int>          simphot_channel;        ///< Array with PMT channel of the simphotons

    //Reconstructed photon information
    std::vector<int>          recphot_channel;        ///< Array with PMT channel of the flash
    float                     recphot_time;           ///< time of the simplebeamflash
    float                     center_of_flash_x;      ///< x Center of opFlash
    float                     center_of_flash_y;      ///< y Center of opFlash
    float                     center_of_flash_z;      ///< z Center of opFlash
    float                     matchscore;             ///< Matchscore of the single opflash with the qcluster containing all the pfparticles

    /* FCL VARIABLES */
    bool                      m_debug;
    float                     m_startbeamtime;
    float                     m_endbeamtime;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FitSimPhotons::FitSimPhotons(fhicl::ParameterSet const & p):EDAnalyzer(p)
{
    simphot_time.resize(m_geo->NOpDets(),0.00);
    simphot_channel.resize(m_geo->NOpDets(),0);
    recphot_channel.resize(m_geo->NOpDets(),0);

    //initialize fcl parameters
    m_debug         = p.get<bool> ("DebugMode"          ,false );
    m_startbeamtime = p.get<float>("FlashVetoTimeStart" ,3.2   );
    m_endbeamtime   = p.get<float>("FlashVetoTimeEnd"   ,4.8   );

    m_mgr.Configure(  p.get<flashana::Config_t>("FlashMatchConfig"));

    //initialize output tree
    art::ServiceHandle<art::TFileService> tfs;
    m_tree  = tfs->make<TTree>("FitSimPhotons","FitSimPhotons Tree");

    //Set branches for (Run Subrun Event)
    m_tree->Branch("run",          &run,           "run/s"       );
    m_tree->Branch("subrun",       &subrun,        "subrun/s"    );
    m_tree->Branch("event",        &event,         "event/i"     );

    //Set branches for PandoraNU information
    m_tree->Branch("q_z_total",     &q_z_total,    "q_z_total/F" );
    m_tree->Branch("q_z_sps",       &q_z_sps,      "q_z_sps/F"   );
    m_tree->Branch("q_z_hit",       &q_z_hit,      "q_z_hit/F"   );
    m_tree->Branch("flashhypo_channel","std::vector<double>",  &flashhypo_channel );
    m_tree->Branch("flashhypo_time",   &flashhypo_time,        "flashhypo_time/F" );
    m_tree->Branch("center_of_charge_x",   &center_of_charge_x,        "center_of_charge_x/F" );
    m_tree->Branch("center_of_charge_y",   &center_of_charge_y,        "center_of_charge_y/F" );
    m_tree->Branch("center_of_charge_z",   &center_of_charge_z,        "center_of_charge_z/F" );


    //Set branches for MC information
    m_tree->Branch("true_pdg",        &true_pdg,              "true_pdg/s"        );
    m_tree->Branch("true_energy",     &true_energy,           "true_energy/F"     );
    m_tree->Branch("true_time",       &true_time,             "true_time/F"       );
    m_tree->Branch("true_x",          &true_x,                "true_x/F"          );
    m_tree->Branch("true_y",          &true_y,                "true_y/F"          );
    m_tree->Branch("true_z",          &true_z,                "true_z/F"          );
    m_tree->Branch("simphot_time",    "std::vector<float>",   &simphot_time       );
    m_tree->Branch("simphot_channel", "std::vector<int>",     &simphot_channel    );

    //Reconstructed photon information
    m_tree->Branch("recphot_channel",     "std::vector<int>",        &recphot_channel      );
    m_tree->Branch("recphot_time",        &recphot_time,             "recphot_time/F"      );
    m_tree->Branch("center_of_flash_x",   &center_of_flash_x,        "center_of_flash_x/F" );
    m_tree->Branch("center_of_flash_y",   &center_of_flash_y,        "center_of_flash_y/F" );
    m_tree->Branch("center_of_flash_z",   &center_of_flash_z,        "center_of_flash_z/F" );
    m_tree->Branch("matchscore",          &matchscore,               "matchscore/F"        );

}

void FitSimPhotons::clearTreeVar()
{
    run               = 0;
    subrun            = 0;
    event             = 0;

    //PandoraNu information
    q_z_total         = 0;
    q_z_sps           = 0;
    q_z_hit           = 0;
    flashhypo_channel.clear();
    flashhypo_time    = 0;
    center_of_charge_x= 0;
    center_of_charge_y= 0;
    center_of_charge_z= 0;

    //MC information
    true_energy       = 0;
    true_time         = 0;
    true_x            = 0;
    true_y            = 0;
    true_z            = 0;
    simphot_time.clear();
    std::fill(simphot_channel.begin(), simphot_channel.end(), 0);

    //Reconstructed photon information
    std::fill(recphot_channel.begin(), recphot_channel.end(), 0);
    recphot_time      = 0;
    matchscore        =-1;
    center_of_flash_x = 0;
    center_of_flash_y = 0;
    center_of_flash_z = 0;
}

#endif // FIT_SIM_PHOTONS_H
