// Some notes:
// Check if nfls = 0 skip event because there was no in time flash
// Check min_x_sps_temp!=9999 otherwise, only do Z matching


#include "VertexFlashMatch_module.h"


void VertexFlashMatch::analyze(art::Event const & e)
{
	resetTreeVar();
    //try {
    fillTree(e);
    //} catch(...) {std::cerr<<"Something went wrong filling root tree"<<std::endl;}
    return;
}




/* LEVEL 0 FUCNTIONS */

void VertexFlashMatch::fillTree(art::Event const & e)
{
    std::cout << "---------------------------------------------------------------------" << std::endl;

    // Fill run information
    run    = e.run();
    subrun = e.subRun();
    event  = e.event();

    std::cout << "Run: " << run << ", SubRun: " << subrun << ", Event: " << event;

    if(e.isRealData() )
    {
        std::cout << ", This is data"<< std::endl;
    }
    else
    {
        std::cout << ", This is MC, fill the truth tree"<< std::endl;
        fillTrueTree(e);
    }
    std::cout << std::endl;
	std::vector<flashana::QCluster_t> qcvec      = fillPandoraTree(e);
	flashana::Flash_t flash                      = fillOticalTree (e);
	std::vector<flashana::FlashMatch_t> matchvec = makeMatch(qcvec,flash);
	fillMatchTree(matchvec);

    std::cout<<"variables filled, fill tree"<<std::endl;
    m_tree->Fill();
    std::cout << "---------------------------------------------------------------------" << std::endl;
}




/* LEVEL 1 FUCNTIONS */

void VertexFlashMatch::fillTrueTree(art::Event const & e)
{
    auto const& truth_handle = e.getValidHandle< std::vector< simb::MCTruth > >( "generator" );
    if (truth_handle->size() > 0)
    {
    	simb::MCParticle mcpart;
        if (truth_handle->at(0).NeutrinoSet())
        {
            std::cout << "Neutrino interaction" << std::endl;
            mcpart = truth_handle->at(0).GetNeutrino().Nu();

            true_ccnc = truth_handle->at(0).GetNeutrino().CCNC();
            true_mode = truth_handle->at(0).GetNeutrino().Mode();

        } //neutrino
        else
        {
            std::cout << "No neutrino interaction." << std::endl;
            auto const& mcparticle_handle = e.getValidHandle< std::vector< simb::MCParticle > >( "largeant" );
            mcpart = mcparticle_handle->at(0);
         }

        true_pdg    = mcpart.PdgCode();
        true_energy = mcpart.E();
        true_time   = mcpart.T();
        true_x      = mcpart.Vx();
        true_y      = mcpart.Vy();
        true_z      = mcpart.Vz();
        true_px     = mcpart.Px();
        true_py     = mcpart.Py();
        true_pz     = mcpart.Pz();
        true_sce_x  = true_x - m_sce->GetPosOffsets(true_x, true_y, true_z)[0]+0.7;
        true_sce_y  = true_y + m_sce->GetPosOffsets(true_x, true_y, true_z)[1];
        true_sce_z  = true_z + m_sce->GetPosOffsets(true_x, true_y, true_z)[2];
    }
    auto const& simphot_handle = e.getValidHandle< std::vector< sim::SimPhotons > >("largeant");

    for(auto const& pmtsimvec :  *simphot_handle)
    {
        simphot_spectrum[pmtsimvec.OpChannel()]+=pmtsimvec.size();
    }
}


std::vector<flashana::QCluster_t> VertexFlashMatch::fillPandoraTree(art::Event const & e)
{
	std::vector<flashana::QCluster_t> qcvec;

    auto const& pfparticle_handle =   e.getValidHandle< std::vector< recob::PFParticle > >("pandoraNu");
    art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, e, "pandoraNu");

    nr_pfp = pfparticle_handle->size();
    std::cout << "Fill the Pandora tree, PFP in event: " << nr_pfp << std::endl;

    //The first qvec element should correspond to the sum of all pfparticles in the event to compare 
    std::vector<size_t> allpfplist(nr_pfp);
    std::iota(std::begin(allpfplist), std::end(allpfplist), 0); //0 is the starting number
    calculateChargeCenter(e,allpfplist);
    qcvec.emplace_back(collect3DHits(e,allpfplist));

    for(size_t pfpindex =0; pfpindex< nr_pfp; pfpindex++)
    {
        
        // We are lokking for neutrinos 
        Short_t pdgcode = pfparticle_handle->at(pfpindex).PdgCode();
        //std::cout<< "PFP" << pfpindex << " has code " << pdgcode << std::endl;

        if(pdgcode!=12 and pdgcode!=14) continue;
        nr_nupfp++;

        auto const& vertex_obj = vertex_per_pfpart.at(pfpindex);
        std::vector<double> vertex(3);
        vertex_obj->XYZ(&vertex[0]);

        nuvtxx.emplace_back(vertex[0]);
        nuvtxy.emplace_back(vertex[1]);
        nuvtxz.emplace_back(vertex[2]);
        nupfp_pdg.emplace_back(pdgcode);

        std::vector<size_t> unordered_daugthers;
        traversePFParticleTree(pfpindex,unordered_daugthers,pfparticle_handle); 

        //Fill nr_trck and nr_shwr
        Short_t nr_trck_temp=0;
        Short_t nr_shwr_temp=0;
        for(auto childi : unordered_daugthers)
        {
            Short_t pdgcode = pfparticle_handle->at(childi).PdgCode();
            if(pdgcode==13)
            {
                nr_trck_temp++;
            }
            if(pdgcode==11)
            {
                nr_shwr_temp++;
            }
        }
        nr_trck.emplace_back(nr_trck_temp);
        nr_shwr.emplace_back(nr_shwr_temp);

        calculateChargeCenter(e,unordered_daugthers);
        qcvec.emplace_back(collect3DHits(e,unordered_daugthers));
        std::cout << "PFP neutrino with PDG" << pdgcode <<", " << nr_trck_temp << " tracks and " << nr_shwr_temp <<  " showers."<<std::endl;

    }

	return qcvec;
}


flashana::Flash_t VertexFlashMatch::fillOticalTree(art::Event const & e)
{
	flashana::Flash_t matchflash;
	
	auto const& optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>("simpleFlashBeam");

    // Loop over the beamflashes, convert the one that is in beamtime and has the highest PE to the flash_t
    int maxPE = 0;
    for(auto const& flash : *optical_handle)
    {
        if(flash.Time()>m_startbeamtime && flash.Time()<m_endbeamtime)
        {
        	nr_flash+=1;
            int thisPE =0;
            ::flashana::Flash_t f;
            f.x = f.x_err = 0;
            f.y = flash.YCenter();
            f.z = flash.ZCenter();
            f.y_err = flash.YWidth();
            f.z_err = flash.ZWidth();
            f.pe_v.resize(m_geo->NOpDets());
            f.pe_err_v.resize(m_geo->NOpDets());
            f.time = flash.Time();
            for(unsigned int ipmt=0; ipmt<m_geo->NOpDets() ; ++ipmt)
            {
                unsigned int opdet = m_geo->OpDetFromOpChannel(ipmt);
                reco_spectrum[opdet]+=flash.PE(ipmt);
                f.pe_v[opdet] = flash.PE(ipmt);
                f.pe_err_v[opdet] = sqrt(flash.PE(ipmt));
                thisPE+=flash.PE(ipmt);
            }

            if(thisPE>maxPE)
            {
                recphot_time = flash.Time();
                matchflash=f;
                center_of_flash_y = flash.YCenter();
                center_of_flash_z = flash.ZCenter();
                width_of_flash_y = flash.YWidth();
                width_of_flash_z = flash.ZWidth();
            }
        }

    }
    std::cout << nr_flash << " Flashes in event." << std::endl;
	return matchflash;
}


std::vector<flashana::FlashMatch_t> VertexFlashMatch::makeMatch(std::vector<flashana::QCluster_t> const  & clustervec,  flashana::Flash_t const & flashReco)
{
	std::vector<flashana::FlashMatch_t> matchvec;
	
    ::flashana::Flash_t flashRecoCopy = flashReco; 

    m_mgr.Reset();
    m_mgr.Emplace(std::move(flashRecoCopy));

    for (auto cluster : clustervec)
    {
        flashana::QCluster_t clusterCopy = cluster;
        m_mgr.Emplace(std::move(clusterCopy));
    }
    
    matchvec = m_mgr.Match();

	return matchvec;
}


void VertexFlashMatch::fillMatchTree(std::vector<flashana::FlashMatch_t> const & matchvec)
{
	if(matchvec.size()==0)
    {
        std::cout << "Something must be wrong, no matching results." <<std::endl;
    }
    else
    {
        for(auto match : matchvec)
        {
            tpc_id.emplace_back(match.tpc_id);
            center_of_flash_x.emplace_back(match.tpc_point.x);
            width_of_flash_x.emplace_back(match.tpc_point_err.x);
            matchscore.emplace_back(match.score);
            hypo_spectrum.emplace_back(match.hypothesis);
        }
    }
}




/* LEVEL 2 FUCNTIONS */

std::vector<double> VertexFlashMatch::calculateChargeCenter (art::Event const & e, std::vector<size_t> const & pfplist)
{
	std::vector<double> center(3,0);
	double min_x_sps_temp = 9999; //random big value that will be overwritten

    //std::cerr << "Calculating the center of charge" << std::endl;

    // Get the associations from pfparticle to spacepoint
    auto const& spacepoint_handle  =   e.getValidHandle< std::vector<recob::SpacePoint > > ( "pandoraNu" );
    auto const& pfparticle_handle =   e.getValidHandle< std::vector< recob::PFParticle > >("pandoraNu"  );

    art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart ( pfparticle_handle, e , "pandoraNu" );
    art::FindManyP<recob::Hit >        hits_per_spcpnts   ( spacepoint_handle , e , "pandoraNu" );

    // Variables for the total weight and center of charge
    double totalweight = 0;

    // Loop over the pfparticles, get their space points, and compute the weighted average:
    for (size_t pfpindex : pfplist)
    {

        // Get the associated spacepoints:
        std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfpindex);

        // Loop over the spacepoints and get the associated hits:
        for (auto & _sps : spcpnts)
        {
            auto xyz = _sps->XYZ();
            if(xyz[0]>0 && xyz[0]<min_x_sps_temp)
            {
                min_x_sps_temp=xyz[0];
            }

            std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
            // Add the hits to the weighted average, if they are collection hits:
            for (auto & hit : hits)
            {
                if (hit->View() == geo::kZ)
                {
                    // Collection hits only
                    float weight = hit->Integral();
                    center[0] += (xyz[0]) * weight;
                    center[1] += (xyz[1]) * weight;
                    center[2] += (xyz[2]) * weight;
                    totalweight += weight;
                } // if collection

            } // hits

        } // spacepoints

    } // pfparticles


    // Normalize;
    center[0] /= totalweight;
    center[1] /= totalweight;
    center[2] /= totalweight;

    // Store the data:
    min_x_sps.emplace_back(min_x_sps_temp);
    q_Y_sps.emplace_back(totalweight);
    center_of_charge_x.emplace_back(center[0]);
    center_of_charge_y.emplace_back(center[1]);
    center_of_charge_z.emplace_back(center[2]);

	return center;
}


flashana::QCluster_t VertexFlashMatch::collect3DHits(art::Event const & e, std::vector<size_t> const & pfplist)
{
	flashana::QCluster_t cluster;
	
	auto const& pfparticle_handle = e.getValidHandle< std::vector< recob::PFParticle > >("pandoraNu");
	auto const& spacepoint_handle = e.getValidHandle<std::vector<recob::SpacePoint>>("pandoraNu");

	art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart   ( pfparticle_handle, e, "pandoraNu" );
	art::FindManyP<recob::Hit > hits_per_spcpnts            ( spacepoint_handle , e, "pandoraNu" );

	for(auto & pfpindex : pfplist)
	{
		unsigned short pdgcode   = pfparticle_handle->at(pfpindex).PdgCode();
	    double lycoef = m_ly_map[pdgcode];

	    std::vector<flashana::Hit3D_t> hitlist;
	    std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfpindex);

	    // Loop over the spacepoints and get the associated hits:
	    for (auto & _sps : spcpnts)
	    {
	        std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
	        // Add the hits to the weighted average, if they are collection hits:
	        for (auto & hit : hits)
	        {
	            if (hit->View() == geo::kZ)
	            {
	                // Collection hits only
	                auto xyz = _sps->XYZ();
	                flashana::Hit3D_t hit3D; //Collection plane hits
	                hit3D.x = xyz[0];
	                hit3D.y = xyz[1];
	                hit3D.z = xyz[2];
	                hit3D.plane =  2;
	                double q = hit->Integral();
                    if (m_normalized)
                    {
                        q*=lycoef;
                    }
                    hit3D.q = q;

	                hitlist.emplace_back(hit3D);
	            }
	        }
	    }
    	cluster+= ((flashana::LightCharge*)(m_mgr.GetCustomAlgo("LightCharge")))->FlashHypothesisCharge(hitlist, lycoef);   
	}
	return cluster;
}




/* LEVEL 3 FUCNTIONS */

void VertexFlashMatch::traversePFParticleTree (size_t top_index, std::vector<size_t> & unordered_daugthers,const art::ValidHandle<std::vector<recob::PFParticle> > pfparticle_handle)
{
	unordered_daugthers.push_back(top_index);
	if (pfparticle_handle->at(top_index).Daughters().size() == 0) 
	{
		return;
	}

	// Else, recursive:
	for (size_t i = 0; i < pfparticle_handle->at(top_index).Daughters().size(); i ++) 
	{
		traversePFParticleTree(pfparticle_handle->at(top_index).Daughters().at(i), unordered_daugthers ,pfparticle_handle);
	}
}


DEFINE_ART_MODULE(VertexFlashMatch)