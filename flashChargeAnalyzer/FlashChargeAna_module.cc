#include "FlashChargeAna_module.h"

namespace flashcharge {


void FlashChargeAna::analyze(art::Event const & e)
{
  try {
  fillTree(e);
  } catch(...) {std::cerr<<"Something went wrong filling root tree"<<std::endl;}
  return;
}











void collect3DHitsZ(    std::vector<flashana::Hit3D_t> &hitlist, 
                        size_t pfindex, 
                        const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
                        const art::Event & evt)
{
  art::InputTag pandoraNu_tag { "pandoraNu" };
  auto const& spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(pandoraNu_tag);

  art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart ( pfparticles, evt, pandoraNu_tag );
  art::FindManyP<recob::Hit > hits_per_spcpnts ( spacepoint_handle, evt, pandoraNu_tag );

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
}


void FlashChargeAna::fillTree(art::Event const & e)
{
  // Fill run information
  run    = e.run();
  subrun = e.subRun();
  event  = e.event();

  std::cout<<"\n Begin filling variables of (run,subrun,event) \t ("<< run <<"," <<subrun <<"," <<event<< ")"<< std::endl;
  // Fill truth information
  if(!e.isRealData() && bool_truth){fillTruthTree(e);}
  // Fill PandoraNu information
  if(bool_pandora){fillPandoraTree(e);}
  // Fill optical information
  if(bool_optical){fillOticalTree(e);}

  std::cout<<"variables filled, fill tree"<<std::endl;
  fTree->Fill();
}

bool FlashChargeAna::is_fiducial(const std::vector<double> & x) const {
  if (x.size() != 3) {
    return false;
  }

  art::ServiceHandle<geo::Geometry> geo;
  std::vector<double> bnd = {0., 2.*geo->DetHalfWidth(), -geo->DetHalfHeight(), geo->DetHalfHeight(), 0., geo->DetLength()};

  bool is_x = x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
  bool is_y = x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
  bool is_z = x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
  return is_x && is_y && is_z;
}

void FlashChargeAna::traversePFParticleTree(
  const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
  size_t top_index,
  std::vector<size_t> & unordered_daugthers )
{

  // This is a tree-traversal algorithm.  It returns the index of the top particle, plus the index
  // of all daughter particles.

  // This is a recursive algorithm, so it needs a break clause:

  if (pfparticles->at(top_index).Daughters().size() == 0) {
    unordered_daugthers.push_back(top_index);
    return;
  }

  // Else, go through the tree:
  for (size_t i = 0; i < pfparticles->at(top_index).Daughters().size(); i ++) {
    traversePFParticleTree(pfparticles, pfparticles->at(top_index).Daughters().at(i), unordered_daugthers );
  }
  return;
}

// Method to calculate the total the center for a parent particle (index of neutrino pfp)
TVector3 FlashChargeAna::calculateChargeCenter(
  size_t top_particle_index,
  const art::ValidHandle<std::vector<recob::PFParticle> > pfparticles,
  const art::Event & evt)
{

  // First, get the indexes of pfparticles that are in the hierarchy of this particle:
  std::vector<size_t> daughters;
  daughters.reserve(50);
  traversePFParticleTree( pfparticles, top_particle_index, daughters);


  // Get the associations from pfparticle to spacepoint
  art::InputTag pandoraNu_tag { "pandoraNu" };
  auto const& spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(pandoraNu_tag);

  art::FindManyP<recob::SpacePoint > spcpnts_per_pfpart ( pfparticles, evt, pandoraNu_tag );
  art::FindManyP<recob::Hit > hits_per_spcpnts ( spacepoint_handle, evt, pandoraNu_tag );


  // Variables for the total weight and center of charge
  double totalweight = 0;
  std::vector<double> chargecenter;
  chargecenter.resize(3);

  // Loop over the pfparticles, get their space points, and compute the weighted average:

  for (auto & _i_pfp : daughters) {

    // Get the associated spacepoints:
    std::vector<art::Ptr < recob::SpacePoint > > spcpnts = spcpnts_per_pfpart.at(_i_pfp);

    // Loop over the spacepoints and get the associated hits:
    for (auto & _sps : spcpnts) {
      std::vector<art::Ptr<recob::Hit> > hits = hits_per_spcpnts.at(_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto & hit : hits) {
        if (hit->View() == geo::kZ) {
          // Collection hits only
          double weight = hit->Integral();
          // std::cout << "Hit Integral: " << hit->Integral() << std::endl;
          // std::cout << "Hit PeakAmplitude: " << hit->PeakAmplitude() << std::endl;
          // std::cout << "Hit SummedADC: " << hit->SummedADC() << std::endl;
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
  TVector3 _center_of_charge;
  _center_of_charge.SetX(chargecenter[0]);
  _center_of_charge.SetY(chargecenter[1]);
  _center_of_charge.SetZ(chargecenter[2]);

  return _center_of_charge;

}


double FlashChargeAna::trackEnergy(const art::Ptr<recob::Track>& track, const art::Event & evt)
{
  art::InputTag pandoraNu_tag { "pandoraNu" };
  auto const& track_handle = evt.getValidHandle< std::vector< recob::Track > >( pandoraNu_tag );
  art::FindManyP<anab::Calorimetry> calo_track_ass(track_handle, evt, "pandoraNucalo");
  const std::vector<art::Ptr<anab::Calorimetry>> calos = calo_track_ass.at(track.key());

  for (size_t ical = 0; ical < calos.size(); ++ical)
  {
    if (!calos[ical]) continue;
    if (!calos[ical]->PlaneID().isValid) continue;
    int planenum = calos[ical]->PlaneID().Plane;
    if (planenum < 0 || planenum > 2) continue;
    if (planenum != 2) continue;                           // Use informartion from collection plane only

    double dedxsum=0;
    double dedx=0;
    double counter=0;
    for (size_t iTrkHit = 0; iTrkHit < calos[ical]->dEdx().size(); ++iTrkHit)
    {
      dedx = calos[ical]->dEdx()[iTrkHit];
      if (dedx > 0 && dedx < 10)
      {
        dedxsum+=dedx;
        counter++;
      }
    }
    return dedxsum/counter;
  }
  return 0.0;
}


void FlashChargeAna::fillTruthTree(art::Event const & e)
{
  std::cout << "Filling truth information " << std::endl;
  art::InputTag truth_tag { "generator" };
  auto const& truth_handle = e.getValidHandle< std::vector< simb::MCTruth > >( truth_tag );

  true_fid = false;
  mcevts_truth=0;
  nuPDG_truth.clear();
  ccnc_truth.clear();
  mode_truth.clear();
  enu_truth.clear();
  nuvtxx_truth.clear();
  nuvtxy_truth.clear();
  nuvtxz_truth.clear();
  nuvtxx_sc.clear();
  nuvtxy_sc.clear();
  nuvtxz_sc.clear();

  if (truth_handle->size() > 0) {
    for(unsigned int iList = 0; iList < truth_handle->size() ; ++iList){
      if (truth_handle->at(iList).NeutrinoSet())
      {
        simb::MCNeutrino const& neutrino = truth_handle->at(iList).GetNeutrino();
        mcevts_truth++;
        nuPDG_truth.emplace_back(neutrino.Nu().PdgCode());
        ccnc_truth.emplace_back(neutrino.CCNC());
        mode_truth.emplace_back(neutrino.Mode());
        enu_truth.emplace_back(neutrino.Nu().E());

        double vx = neutrino.Nu().Vx();
        double vy = neutrino.Nu().Vy();
        double vz = neutrino.Nu().Vz();

        std::vector<double> vertex = {vx,vy,vz};
        if(is_fiducial(vertex))
        {
          true_fid = true;
        }

        nuvtxx_truth.emplace_back(vx);
        nuvtxy_truth.emplace_back(vy);
        nuvtxz_truth.emplace_back(vz);

        nuvtxx_sc.emplace_back(vx-SCE.GetPosOffsets(vx, vy, vz)[0]+0.7);
        nuvtxy_sc.emplace_back(vy+SCE.GetPosOffsets(vx, vy, vz)[1]);
        nuvtxz_sc.emplace_back(vz+SCE.GetPosOffsets(vx, vy, vz)[2]);
      }
    }
  }
}


void FlashChargeAna::fillPandoraTree(art::Event const & e)
{
  std::cout << "Filling PandoraNu information " << std::endl;

  nnuvtx=0;
  nuvtxx.clear();
  nuvtxy.clear();
  nuvtxz.clear();
  center_of_charge_x.clear();
  center_of_charge_y.clear();
  center_of_charge_z.clear();
  nuvtxpdg.clear();

  shwr_dir.clear();
  shwr_en.clear();
  shwr_angle.clear();

  trck_dir.clear();
  trck_dedxavg.clear();
  trck_len.clear();

  art::InputTag pandoraNu_tag { "pandoraNu" };
  auto const& pfparticle_handle = e.getValidHandle< std::vector< recob::PFParticle > >( pandoraNu_tag );
  art::FindOneP< recob::Shower > shower_per_pfpart(pfparticle_handle, e, pandoraNu_tag);
  art::FindOneP< recob::Track > track_per_pfpart(pfparticle_handle,   e, pandoraNu_tag);
  art::FindOneP< recob::Vertex > vertex_per_pfpart(pfparticle_handle, e, pandoraNu_tag);

  for (size_t pfpindex =0; pfpindex < pfparticle_handle->size() ; ++pfpindex )
  {
    recob::PFParticle const& pfp = pfparticle_handle->at(pfpindex);
    auto const& vertex_obj = vertex_per_pfpart.at(pfpindex);
    std::vector<double> vertex(3);
    vertex_obj->XYZ(&vertex[0]);

    if(!is_fiducial(vertex)) continue;
    if(pfp.PdgCode()!=12 and pfp.PdgCode()!=14) continue;

    nnuvtx++;

    nuvtxx.emplace_back(vertex[0]);
    nuvtxy.emplace_back(vertex[1]);
    nuvtxz.emplace_back(vertex[2]);

    TVector3 center_of_charge = calculateChargeCenter(pfpindex,pfparticle_handle,e);
    center_of_charge_x.emplace_back(center_of_charge.X());
    center_of_charge_y.emplace_back(center_of_charge.Y());
    center_of_charge_z.emplace_back(center_of_charge.Z());
    
    nuvtxpdg.emplace_back(pfp.PdgCode());


    // Track Shower info

    std::vector<TVector3> shwr_dir_primary;
    std::vector<double> shwr_angle_primary;
    std::vector<double> shwr_en_primary;

    std::vector<TVector3> trck_dir_primary;
    std::vector<double> trck_dedxavg_primary;
    std::vector<double> trck_len_primary;

    std::vector<size_t> daugthers; 
    traversePFParticleTree( pfparticle_handle, pfpindex,daugthers );

    for (size_t childi : daugthers)
    {
      recob::PFParticle const& child = pfparticle_handle->at(childi);

      if (child.PdgCode()==11)
      {
        auto const& shower_obj = shower_per_pfpart.at(childi);
        shwr_dir_primary.emplace_back(shower_obj->Direction());
        shwr_en_primary.emplace_back(*std::max_element(shower_obj->Energy().begin(),shower_obj->Energy().end()));
        shwr_angle_primary.emplace_back(shower_obj->OpenAngle());
      }
      if (child.PdgCode()==13)
      {
        auto const& track_obj   = track_per_pfpart.at(childi);
        trck_dir_primary.emplace_back(TVector3 (track_obj->Direction().second.x(),track_obj->Direction().second.y(),track_obj->Direction().second.z()));
        trck_len_primary.emplace_back(track_obj->Length());
        trck_dedxavg_primary.emplace_back(trackEnergy(track_obj, e));
      }

    }

    shwr_dir.emplace_back(shwr_dir_primary);
    shwr_en.emplace_back(shwr_en_primary);
    shwr_angle.emplace_back(shwr_angle_primary);

    trck_dir.emplace_back(trck_dir_primary);
    trck_len.emplace_back(trck_len_primary);
    trck_dedxavg.emplace_back(trck_dedxavg_primary);
  }
}


void FlashChargeAna::fillOticalTree(art::Event const & e){
  std::cout << "Filling optical information " << std::endl;
  art::InputTag optical_tag{"simpleFlashBeam"};
  auto const& optical_handle = e.getValidHandle<std::vector<recob::OpFlash>>(optical_tag);
  if(optical_handle->size())
  {
    nfls=0;
    flsTime.clear();
    flsPe.clear();
    flsYcenter.clear();
    flsZcenter.clear();
    flsYwidth.clear();
    flsZwidth.clear();
    flsPePMT.clear();

    nfls = optical_handle->size();
    for(int ifl=0; ifl< nfls; ++ifl)
    {
      recob::OpFlash const& flash = optical_handle->at(ifl);
      flsTime.emplace_back(flash.Time());
      flsPe.emplace_back(flash.TotalPE());
      flsYcenter.emplace_back(flash.YCenter());
      flsZcenter.emplace_back(flash.ZCenter());
      flsYwidth.emplace_back(flash.YWidth());
      flsZwidth.emplace_back(flash.ZWidth());
      for(int ipmt=0; ipmt < m_nrPMT ; ++ipmt)
      {
      	flsPePMT.emplace_back(flash.PE(ipmt));
      }
    }
  }
}

}

DEFINE_ART_MODULE(flashcharge::FlashChargeAna)
