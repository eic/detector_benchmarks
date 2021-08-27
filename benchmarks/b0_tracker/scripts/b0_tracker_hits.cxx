R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"

#include "ROOT/RDataFrame.hxx"
#include "Math/Vector3D.h"
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

#include "dd4pod/TrackerHitCollection.h"

#include "common_bench/particles.h"
#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"


#include "dd4pod/TrackerHitCollection.h"

void b0_tracker_hits(const char* fname = "./sim_output/sim_forward_protons.root"){

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel
  double degree = TMath::Pi()/180.0;

  TChain* t = new TChain("events");
  t->Add(fname);

  ROOT::RDataFrame d0(*t);

  auto hits_theta = [&](const std::vector<dd4pod::TrackerHitData>& hits) {
    std::vector<double> result;
    for (const auto& h : hits) {
      ROOT::Math::XYZVector vec(h.position.x,h.position.y,h.position.z);
      result.push_back(1000*vec.theta());
      std::cout << 1000*vec.theta() << "\n";
    }
    return result;
  };

  auto local_position = [&](const std::vector<dd4pod::TrackerHitData>& hits) {
    std::vector<std::array<double, 2>> result;
    for (const auto& h : hits) {

      auto pos0 = (h.position);
	  result.push_back({pos0.x , pos0.y});
    }
    return result;
  };


  auto x_pos = [&](const std::vector<std::array<double, 2>>& xypos) {
    std::vector<double> result;
    for (const auto& h : xypos) {

		result.push_back(h.at(0));
    }
    return result;
  };

  auto y_pos = [&](const std::vector<std::array<double, 2>>& xypos) {
    std::vector<double> result;
    for (const auto& h : xypos) {

        result.push_back(h.at(1));
    }
    return result;
  };


  auto d1 = d0.Define("nhits", hits_theta, {"B0TrackerHits"})
			  .Define("xy_hit_pos", local_position, {"B0TrackerHits"})
              .Define("x_pos", x_pos, {"xy_hit_pos"})
              .Define("y_pos", y_pos, {"xy_hit_pos"});

  auto h_local_pos = d1.Histo2D({"h_local_pos", ";x [mm]; y [mm] ", 100,  -100.0, -200.0, 100, -100.0, 100.0}, "x_pos", "y_pos");

  auto d2 = d0.Define("hits_theta", hits_theta, {"B0TrackerHits"});

  auto h1 = d2.Histo1D({"h1", "hits_theta", 100, 0,20}, "hits_theta");
  TCanvas* c = new TCanvas();
  h1->DrawCopy();
  c->SaveAs("results/b0_tracker_hits_theta.png");
  c->SaveAs("results/b0_tracker_hits_theta.pdf");
  
  h_local_pos->DrawCopy("colz");
  c->SaveAs("results/b0_tracker_hits_occupancy_disk_1.png");
  c->SaveAs("results/b0_tracker_hits_occupancy_disk_1.pdf");

  auto n1 = h1->GetMean();
  std::cout << "Polar angle of hits: " << n1 << std::endl;

  //if (n1 < 5) {
  //        std::quick_exit(1);
  //}

}

