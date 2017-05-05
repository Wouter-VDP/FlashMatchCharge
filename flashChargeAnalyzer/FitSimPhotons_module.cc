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
	auto const& spacepoint_handle = e.getValidHandle< std::vector< recob::SpacePoint > >( "pandoraNu" );

	art::FindOneP< recob::Shower >      shower_per_pfpart  (pfparticle_handle, e, "pandoraNu");
    art::FindOneP< recob::Track >       track_per_pfpart   (pfparticle_handle, e, "pandoraNu");
    art::FindManyP< recob::SpacePoint > spcpnts_per_pfpart (pfparticle_handle, e, "pandoraNu");
    art::FindManyP< recob::Hit >         hits_per_spcpnts ( spacepoint_handle, e, "pandoraNu");
    art::FindManyP< recob::Hit >        hits_per_pfpart    (pfparticle_handle, e, "pandoraNu");


    for (size_t pfpindex =0; pfpindex < pfparticle_handle->size() ; ++pfpindex )
  	{
  		// Total charge through spacepoints
    	std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfpindex);
    	for (auto & _sps : spcpnts) 
    	{
      		std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
      		for (auto & hit : hits) {
        		if (hit->View() == geo::kZ) 
        		{
        			q_z_sps += hit->Integral();
        		}
        	}
        }
        // Total charge through hits directly, no 3D info this way!
        std::vector<art::Ptr < recob::Hit > > hits = hits_per_pfpart.at(pfpindex);
    	for (auto & _hit : hits) 
    	{
    		if (_hit->View() == geo::kZ) 
    		{
    			q_z_hit += _hit->Integral();
    		}
    	}
    }
}

void FitSimPhotons::fillOticalTree(art::Event const & e)
{
	simphot_time.clear();
	simphot_channel.clear();

    art::ServiceHandle<geo::Geometry> geo; 

    auto const& simphot_handle = e.getValidHandle< std::vector< sim::SimPhotons > >( "largeant" );

    for(size_t opdet=0; opdet<geo->NOpDets(); ++opdet) 
    {
        sim::SimPhotons const& simph = simphot_handle->at(opdet);
        for(auto const& oneph : simph) 
        {
            simphot_channel.emplace_back(opdet);
            simphot_time.emplace_back(oneph.Time);
        }
    }
}

DEFINE_ART_MODULE(FitSimPhotons)

