//R__LOAD_LIBRARY(libfmt.so)
//#include "fmt/core.h"
R__LOAD_LIBRARY(libDDG4IO.so)
//
//#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
//#include "DDRec/CellIDPositionConverter.h"
//#include "DDRec/SurfaceManager.h"
//#include "DDRec/Surface.h"
#include "ROOT/RDataFrame.hxx"
//
//#include "lcio2/MCParticleData.h"
//#include "lcio2/ReconstructedParticleData.h"

//#include "Math/Vector3D.h"
//#include "Math/Vector4D.h"
//#include "Math/VectorUtil.h"
#include "TCanvas.h"
//#include "TLegend.h"
//#include "TMath.h"
//#include "TRandom3.h"
//#include "TFile.h"
//#include "TH1F.h"
//#include "TH1D.h"
//#include "TTree.h"
#include "TChain.h"
//#include "TF1.h"
#include <random>
//#include "lcio2/TrackerRawDataData.h"
//#include "lcio2/TrackerRawData.h"

void simple_checking(const char* fname = "./sim_output/output_zdc_photons.root"){
std::cout << "testing 1\n";
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel
  //using namespace lcio2;
  double degree = TMath::Pi()/180.0;

  TChain* t = new TChain("EVENT");
  t->Add(fname);

  ROOT::RDataFrame d0(*t);//, {"GEMTrackerHintits","MCParticles"});

  //std::cout << t->GetBranch("GEMTrackerHits")->GetClassName() << std::endl;
  //std::vector<dd4hep::sim::Geant4Tracker::Hit*>
  
  // -------------------------
  // Get the DD4hep instance
  // Load the compact XML file
  // Initialize the position converter tool
  //dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  //detector.fromCompact("gem_tracker_disc.xml");
  //dd4hep::rec::CellIDPositionConverter cellid_converter(detector);

  //// -------------------------
  //// Get the surfaces map
  //dd4hep::rec::SurfaceManager& surfMan = *detector.extension<dd4hep::rec::SurfaceManager>() ;
  //auto surfMap = surfMan.map( "world" ) ;
  
  auto nhits = [] (std::vector<dd4hep::sim::Geant4Calorimeter::Hit*>& hits){ return (int) hits.size(); };
  //auto hit_position = [&](const std::vector<dd4hep::sim::Geant4Tracker::Hit*>& hits){
  //for(const auto& h: hits){
  //  //std::cout << (h->position/10.0) << std::endl;
  //  //std::cout << cellid_converter.position(h->cellID) << std::endl;
  //  //dd4hep::rec::SurfaceMap::const_iterator
  //  const auto si = _surfMap.find( cellid_converter.findContext(h->cellID)->identifier ); //identifier=volumeID
  //  dd4hep::rec::ISurface* surf = (si != _surfMap.end() ?  si->second  : 0);
  //  dd4hep::rec::Vector3D pos =  surf->origin();//fit_global(pivot[0],pivot[1],pivot[2]);
  //  //std::cout << pos.x() << ", " << pos.y() << ", " << pos.z()<< std::endl;
  //  // transform lcio units to dd4hep units, see documentation for other functions              
  //  //DDSurfaces::Vector2D fit_local = surf->globalToLocal( dd4hep::mm * fit_global );
  //}
  //  return hits.size(); };

  //auto digitize_gem_hits = 
  //  [&](const std::vector<dd4hep::sim::Geant4Tracker::Hit*>& hits) {
  //    std::vector<lcio2::TrackerRawDataData> digi_hits;
  //    std::normal_distribution<> time_dist(0,2.0);
  //    std::normal_distribution<> adc_dist(5.0,3.0);

  //    std::map<int64_t,lcio2::TrackerRawDataData> hits_by_id;
  //    for(const auto& h: hits) {
  //      //lcio2::TrackerRawDataData ahit;
  //      auto& ahit = hits_by_id[(int64_t)h->cellID];
  //      auto pos = h->position/10.0; //cm

  //      ahit.cellID0   = h->cellID;
  //      ahit.cellID1   = h->cellID;
  //      ahit.channelID = h->cellID;
  //      //fmt::print("{} vs {} vs {}\n", id1, id2,id3);
  //      fmt::print("{} vs {}\n", h->cellID, ahit.cellID0);
  //      // time is not kept from dd4hep hit, instead using z position as crude substitute
  //      ahit.time = pos.z() + time_dist(gen);
  //      ahit.adc  = adc_dist(gen);
  //      //digi_hits.push_back(ahit);
  //    }
  //    for (auto& [cell, hit] : hits_by_id) {
  //      //fmt::print("{} vs {}\n", cell, hit.cellID0);
  //      digi_hits.push_back(hit);
  //    }
  //    return digi_hits;
  //  };

  auto d1 = d0.Define("nhits", nhits, {"ZDCHits"})
              //.Filter([](int n){ return (n>4); },{"nhits"})
              //.Define("delta",hit_position, {"GEMTrackerHits"})
              //.Define("RawTrackerHits", digitize_gem_hits, {"GEMTrackerHits"})
              ;

  auto h0 = d1.Histo1D(TH1D("h0", "nhits; ", 20, 0,20), "nhits");

  auto n0 = d1.Filter([](int n){ return (n>0); },{"nhits"}).Count();

  TCanvas* c = new TCanvas();

  //d1.Snapshot("digitized_EVENT","test_gem_tracker_digi.root");
  h0->DrawClone();
  std::cout << *n0 << " events with nonzero hits\n";

  if(*n0<5) {
    std::quick_exit(1);
  }

}

