#include "FitSimPhotons_module.h"


void FitSimPhotons::analyze(art::Event const & e)
{
	//try {
	fillTree(e);
	//} catch(...) {std::cerr<<"Something went wrong filling root tree"<<std::endl;}
	return;
}

void FitSimPhotons::fillTree(art::Event const & e)
{
	// Fill run information
	run    = e.run();
	subrun = e.subRun();
	event  = e.event();

	std::cout<<"\n Begin filling variables of (run,subrun,event) \t ("<< run <<"," <<subrun <<"," <<event<< ")"<< std::endl;

	fillPandoraTree(e);       // Fill PandoraNu information
	fillOticalTree(e);        // Fill optical information

	std::cout<<"variables filled, fill tree"<<std::endl;
    m_tree->Fill();
}

void FitSimPhotons::fillPandoraTree(art::Event const & e)
{
	std::cout << "Filling PandoraNu information " << std::endl;

	q_z_total = 0;
	q_z_sps   = 0;
	q_z_hit   = 0;

	auto const& pfparticle_handle = e.getValidHandle< std::vector< recob::PFParticle > >( "pandoraNu" );
	art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, e, "pandoraNu");
    art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle,   e, "pandoraNu");
}

void FitSimPhotons::fillOticalTree(art::Event const & e)
{
	simphot_time.clear();
	simphot_channel.clear();
}
