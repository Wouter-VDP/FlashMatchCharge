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
    auto const& cluster_handle    = e.getValidHandle< std::vector< recob::Cluster    > >( "pandoraNu" );
	auto const& spacepoint_handle = e.getValidHandle< std::vector< recob::SpacePoint > >( "pandoraNu" );
    auto const& hit_handle        = e.getValidHandle< std::vector< recob::Hit        > >( "gaushit"   );  

	art::FindOneP< recob::Shower >      shower_per_pfpart  (pfparticle_handle, e, "pandoraNu");
    art::FindOneP< recob::Track >       track_per_pfpart   (pfparticle_handle, e, "pandoraNu");
    art::FindManyP< recob::SpacePoint > spcpnts_per_pfpart (pfparticle_handle, e, "pandoraNu");
    art::FindManyP< recob::Hit >        hits_per_spcpnts   (spacepoint_handle, e, "pandoraNu");
    art::FindManyP< recob::Cluster >    clusters_per_pfpart(pfparticle_handle, e, "pandoraNu");
    art::FindManyP< recob::Hit >        hits_per_cluster   (cluster_handle,    e, "pandoraNu");

    for (size_t pfpindex =0; pfpindex < pfparticle_handle->size() ; ++pfpindex )
  	{
        int PDGcode = pfparticle_handle->at(pfpindex).PdgCode();
  		// Total charge through spacepoints
        if(PDGcode!=11 && PDGcode!=13) continue;
    	std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfpindex);

        if(m_debug)
        {
            std::cout << "PFParticle" << pfpindex <<", with PDGcode " <<  PDGcode <<", has " << spcpnts.size() << " spacepoints " << std::endl;
        }

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
        // Total charge through clusters
        std::vector<art::Ptr < recob::Cluster > > clusters = clusters_per_pfpart.at(pfpindex);
        if(m_debug)
        {
            std::cout << "PFParticle" << pfpindex <<", with PDGcode " <<  PDGcode <<", has " << clusters.size() << " clusters " << std::endl;
        }
    	for (auto & _cluster : clusters) 
    	{
    		if (_cluster->View() == geo::kZ) 
    		{
    			q_z_hit += _cluster->Integral();
    		}
    	}
    }
    // Total charge through hits
    for (size_t hitindex =0; hitindex <hit_handle->size() ; ++hitindex )
    {
        if (hit_handle->at(hitindex).View() == geo::kZ){
            q_z_total+=   hit_handle->at(hitindex).Integral();
        }
    }

    if(m_debug)
    {
        std::cout << "Number PFParticles in event: \t" <<   pfparticle_handle->size() << std::endl;
        std::cout << "q_z_total: \t" <<   q_z_total << std::endl;
        std::cout << "q_z_sps : \t" <<   q_z_sps << std::endl;
        std::cout << "q_z_hit: \t" <<   q_z_hit << std::endl;
    }
}

void FitSimPhotons::fillOticalTree(art::Event const & e)
{
	simphot_time.clear();
	std::fill(simphot_channel.begin(), simphot_channel.end(), 0);
    std::fill(recphot_channel.begin(), recphot_channel.end(), 0); 

    auto const& simphot_handle = e.getValidHandle< std::vector< sim::SimPhotons > >( "largeant" );
    auto const& optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>("simpleFlashBeam");

    for(auto const& pmtsimvec :  *simphot_handle) 
    {

        simphot_channel[pmtsimvec.OpChannel()]+=pmtsimvec.size();
        //simphot_time.emplace_back(oneph.Time);
    }

    for(auto const& flash : *optical_handle)
    {
        for(unsigned int ipmt=0; ipmt<m_geo->NOpDets() ; ++ipmt)
        {
            recphot_channel[ipmt]+=flash.PE(ipmt);
        }
    }
    if(m_debug)
    {
        std::cout << "Number of simpleFlashBeam: \t" <<  optical_handle->size() << std::endl;
        std::cout << "Number of simphotons: \t" <<  simphot_handle->size() << std::endl;
        std::cout << "Total simphotons: \t" <<    std::accumulate(simphot_channel.begin(), simphot_channel.end(), 0) << std::endl;
        std::cout << "Total PE: \t" <<    std::accumulate(recphot_channel.begin(), recphot_channel.end(), 0) << std::endl;
    }
}


flashana::QCluster_t FitSimPhotons::collect3DHitsZ(    std::vector<flashana::Hit3D_t> & hitlist, 
                                        size_t pfindex, 
                                        const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                                        art::Event const & e)
{
    auto const& spacepoint_handle = e.getValidHandle<std::vector<recob::SpacePoint>>("pandoraNu");
    art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart   ( pfparticles,       e, "pandoraNu" );
    art::FindManyP<recob::Hit > hits_per_spcpnts            ( spacepoint_handle, e, "pandoraNu" );

    std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfindex);

    // Loop over the spacepoints and get the associated hits:
    for (auto & _sps : spcpnts) {
        std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
        // Add the hits to the weighted average, if they are collection hits:
        for (auto & hit : hits) {
            if (hit->View() == geo::kZ) {
            // Collection hits only
                auto xyz = _sps->XYZ();
                flashana::Hit3D_t hit3D; //Collection plane hits
                hit3D.x = xyz[0];
                hit3D.y = xyz[1];
                hit3D.z = xyz[2];
                hit3D.plane =  2;
                hit3D.q = hit->Integral();

                hitlist.emplace_back(hit3D);
            }
        }
    }
    //the conversion number is a pure guess! will this depend on x?
    flashana::QCluster_t result = ((flashana::LightCharge*)(m_mgr.GetCustomAlgo("LightCharge")))->FlashHypothesisCharge(hitlist, 1);
    return result;
}

flashana::Flash_t FitSimPhotons::Matching(  size_t pfindex, 
                                            const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                                            art::Event const & e) 
{
    std::vector<flashana::Hit3D_t> hitlist;
    flashana::QCluster_t qcluster = collect3DHitsZ(hitlist, pfindex, pfparticles, e);

    flashana::Flash_t flashHypo;
    flashHypo.pe_v.resize(32);

    return flashHypo;
}



DEFINE_ART_MODULE(FitSimPhotons)

