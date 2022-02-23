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

void b0_tracker_hits(const char* fname = "./sim_output/sim_forward_protons.edm4hep.root"){

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel
  double degree = TMath::Pi()/180.0;

  TChain* t = new TChain("events");
  t->Add(fname);

  ROOT::RDataFrame d0(*t);

  auto hits_eta = [&](const std::vector<dd4pod::TrackerHitData>& hits) {
    std::vector<double> result;
    for (const auto& h : hits) {
      ROOT::Math::XYZVector vec(h.position.x,h.position.y,h.position.z);
      result.push_back(vec.eta());
      std::cout << vec.eta() << "\n";
    }
    return result;
  };

  auto d1 = d0.Define("hits_eta", hits_eta, {"B0TrackerHits"});

  auto h1 = d1.Histo1D({"h1", "hits_eta", 100, 0,20}, "hits_eta");
  TCanvas* c = new TCanvas();
  h1->DrawCopy();
  c->SaveAs("results/b0_tracker_hits_eta.png");
  c->SaveAs("results/b0_tracker_hits_eta.pdf");
  auto n1 = h1->GetMean();
  std::cout << "Pseudorapidity of hits: " << n1 << std::endl;

  //if (n1 < 5) {
  //        std::quick_exit(1);
  //}

}

