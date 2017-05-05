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
#include "lardataobj/Simulation/SimPhotons.h"

#include "TTree.h"
#include "TVector3.h"


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

        void fillTree       (art::Event const & e);
        void fillPandoraTree(art::Event const & e);
        void fillOticalTree (art::Event const & e);


        //Variables

        TTree*                    m_tree;

        /* TREE VARIABLES*/
        //Run Subrun Event
        unsigned short            run;
        unsigned short            subrun;
        unsigned int              event;

        //PandoraNu information
        float                     q_z_total;              ///< The total charge on the collectionplane in beamwindow 
        float                     q_z_sps;                ///< The charge on the collectionplane in beamwindow from PFPtree spacepoints
        float                     q_z_hit;                ///< The charge on the collectionplane in beamwindow from PFPtree hits

        //MC photon information
        std::vector<float>        simphot_time;           ///< Array with times of the simphotons 
        std::vector<short>        simphot_channel;        ///< Array with PMT channel of the simphotons 


        /* FCL VARIABLES */
        bool                      m_debug;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FitSimPhotons::FitSimPhotons(fhicl::ParameterSet const & p):EDAnalyzer(p) 
{
    //initialize fcl parameters
    m_debug         = p.get<bool>("DebugMode",      false);

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

    //Set branches for MC photon information
    m_tree->Branch("simphot_time",    "std::vector<float>",   &simphot_time       );
    m_tree->Branch("simphot_channel", "std::vector<short>",   &simphot_channel    );
}

#endif // FIT_SIM_PHOTONS_H