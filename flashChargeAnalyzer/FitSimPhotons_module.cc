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
    clearTreeVar();

    // Fill run information
    run    = e.run();
    subrun = e.subRun();
    event  = e.event();

    auto const& mcparticle_handle = e.getValidHandle< std::vector< simb::MCParticle > >( "largeant" );
    true_pdg    = mcparticle_handle->at(0).PdgCode();
    true_energy = mcparticle_handle->at(0).E();
    true_time   = mcparticle_handle->at(0).T();
    true_x      = mcparticle_handle->at(0).Vx();
    true_y      = mcparticle_handle->at(0).Vy();
    true_z      = mcparticle_handle->at(0).Vz();

    std::cout<<"\n Begin filling variables of (run,subrun,event) \t ("<< run <<"," <<subrun <<"," <<event<< ") Particle energy: "<< true_energy << "GeV, Time:" << true_time << "?" << std::endl;

    fillPandoraTree(e);       // Fill PandoraNu information
    ::flashana::Flash_t flashReco = fillOticalTree(e);        // Fill optical information
    calculateChargeCenter(e);
    if(flashReco.pe_v.size()==m_geo->NOpDets()) 
    {
        flashana::Flash_t flashHypo = makeMatch(e,flashReco);
        flashhypo_channel = flashHypo.pe_v;    // Fill the resonstructed hypothetical flash corresponding to all pfps
        flashhypo_time    = flashHypo.time;
    }
    else
    {
        std::cout<<"vsomething went wrong, no recoflash object found in beamtime!"<<std::endl;
    }

    std::cout<<"variables filled, fill tree"<<std::endl;
    m_tree->Fill();
}

void FitSimPhotons::fillPandoraTree(art::Event const & e)
{
    std::cout << "Filling PandoraNu information " << std::endl;

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
            for (auto & hit : hits)
            {
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
        if (hit_handle->at(hitindex).View() == geo::kZ)
        {
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

flashana::Flash_t FitSimPhotons::fillOticalTree(art::Event const & e)
{
    std::cout << "Filling OpticalTree information " << std::endl;
    std::fill(recphot_channel.begin(), recphot_channel.end(), 0);

    auto const& simphot_handle = e.getValidHandle< std::vector< sim::SimPhotons > >( "largeant" );
    auto const& optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>("simpleFlashBeam");

    for(auto const& pmtsimvec :  *simphot_handle)
    {

        simphot_channel[pmtsimvec.OpChannel()]+=pmtsimvec.size();
        //simphot_time.emplace_back(oneph.Time);
    }

    ::flashana::Flash_t f;
    for(auto const& flash : *optical_handle)
    {
        if(flash.Time()>m_startbeamtime && flash.Time()<m_endbeamtime)    //Expect that this happens only one time, cetainly of with single particles for now
        {
            center_of_flash_y = flash.YCenter();
            center_of_flash_z = flash.ZCenter();
            f.x = f.x_err = 0;
            f.y = flash.YCenter();
            f.z = flash.ZCenter();
            f.y_err = flash.YWidth();
            f.z_err = flash.ZWidth();
            f.pe_v.resize(m_geo->NOpDets());
            f.pe_err_v.resize(m_geo->NOpDets());
            for(unsigned int ipmt=0; ipmt<m_geo->NOpDets() ; ++ipmt)
            {
                unsigned int opdet = m_geo->OpDetFromOpChannel(ipmt);
                recphot_channel[opdet]+=flash.PE(ipmt);
                f.pe_v[opdet] = flash.PE(ipmt);
                f.pe_err_v[opdet] = sqrt(flash.PE(ipmt));
            }
            std::cout << std::endl;
            recphot_time = flash.Time();
        }

    }
    if(m_debug)
    {
        std::cout << "Number of simpleFlashBeam: \t" <<  optical_handle->size() << std::endl;
        std::cout << "Number of simphotons: \t" <<  simphot_handle->size() << std::endl;
        std::cout << "Total simphotons: \t" <<    std::accumulate(simphot_channel.begin(), simphot_channel.end(), 0) << std::endl;
        std::cout << "Total PE: \t" <<    std::accumulate(recphot_channel.begin(), recphot_channel.end(), 0) << std::endl;
    }
    return f;
}


flashana::QCluster_t FitSimPhotons::collect3DHitsZ( size_t pfindex,
        const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
        art::Event const & e)
{
    std::cout << "Collecting hits information " << std::endl;
    std::vector<flashana::Hit3D_t> hitlist;

    auto const& spacepoint_handle = e.getValidHandle<std::vector<recob::SpacePoint>>("pandoraNu");
    art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart   ( pfparticles,       e, "pandoraNu" );
    art::FindManyP<recob::Hit > hits_per_spcpnts            ( spacepoint_handle, e, "pandoraNu" );

    std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(pfindex);

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
                hit3D.q = hit->Integral();

                hitlist.emplace_back(hit3D);
            }
        }
    }
    //the conversion number is a pure guess! will this depend on x?
    flashana::QCluster_t result = ((flashana::LightCharge*)(m_mgr.GetCustomAlgo("LightCharge")))->FlashHypothesisCharge(hitlist, 1);
    return result;
}

flashana::Flash_t FitSimPhotons::makeMatch(art::Event const & e, flashana::Flash_t & flashReco)
{
    std::cout << "Making the flash hypothesis and matching the qcluster to the opflash" << std::endl;
    auto const& pfparticle_handle = e.getValidHandle< std::vector< recob::PFParticle > >( "pandoraNu" );
    flashana::QCluster_t summed_qcluster;
    summed_qcluster.clear();
    

    for (size_t pfpindex =0; pfpindex < pfparticle_handle->size() ; ++pfpindex )
    {
        summed_qcluster += collect3DHitsZ(pfpindex, pfparticle_handle, e);
    }

    flashana::Flash_t flashHypo;
    flashHypo.pe_v.resize(32);
    ((flashana::PhotonLibHypothesis*)(m_mgr.GetAlgo(flashana::kFlashHypothesis)))->FillEstimate(summed_qcluster,flashHypo);
    double hyposum=0;
    double recosum=0;

    for (unsigned int ipmt = 0; ipmt < m_geo->NOpDets(); ipmt++)
    {
        hyposum+=flashHypo.pe_v[ipmt];
        recosum+=flashReco.pe_v[ipmt];
        std::cout << "for pmt " <<  ipmt << ", hypoPE: " << flashHypo.pe_v[ipmt] << ", recoPE: " << flashReco.pe_v[ipmt] << std::endl;
    }

    if(hyposum>0.01 && recosum>0.01)
    {

        m_mgr.Reset();
        m_mgr.Emplace(std::move(flashReco));
        m_mgr.Emplace(std::move(summed_qcluster));
        m_result = m_mgr.Match();

        if(m_result.size()!=1)
        {
            std::cout << "Something must be wrong, flash result has size " << m_result.size() <<std::endl;
        }
        else
        {
            center_of_flash_x = m_result[0].tpc_point.x;
            matchscore        = m_result[0].score;
        }
    }
    else
    {
         std::cout << "Something must be wrong,the flash hypothesis was 0!" <<std::endl;
    }
    return flashHypo;
}


// Method to calculate the total the center for a parent particle (index of neutrino pfp)
void FitSimPhotons::calculateChargeCenter(const art::Event & e)
{
    std::cout << "Calculating the center of charge" << std::endl;
    // Get the associations from pfparticle to spacepoint
    auto const& spacepoint_handle = e.getValidHandle< std::vector<recob::SpacePoint > > ( "pandoraNu" );
    auto const& pfparticle_handle = e.getValidHandle< std::vector< recob::PFParticle > >( "pandoraNu" );

    art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart ( pfparticle_handle, e, "pandoraNu" );
    art::FindManyP<recob::Hit > hits_per_spcpnts          ( spacepoint_handle, e, "pandoraNu" );


    // Variables for the total weight and center of charge
    double totalweight = 0;
    std::vector<double> chargecenter;
    chargecenter.resize(3);

    // Loop over the pfparticles, get their space points, and compute the weighted average:

    for (size_t pfpindex =0; pfpindex < pfparticle_handle->size() ; ++pfpindex )
    {

        // Get the associated spacepoints:
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
                    double weight = hit->Integral();
                    auto xyz = _sps->XYZ();
                    chargecenter[0] += (xyz[0]) * weight;
                    chargecenter[1] += (xyz[1]) * weight;
                    chargecenter[2] += (xyz[2]) * weight;
                    totalweight += weight;
                    break; // Exit the loop over hits
                } // if collection

            } // hits

        } // spacepoints

    } // pfparticles


    // Normalize;
    chargecenter[0] /= totalweight;
    chargecenter[1] /= totalweight;
    chargecenter[2] /= totalweight;

    // Store the data:
    center_of_charge_x =chargecenter[0];
    center_of_charge_y =chargecenter[1];
    center_of_charge_z =chargecenter[2];
}

DEFINE_ART_MODULE(FitSimPhotons)

