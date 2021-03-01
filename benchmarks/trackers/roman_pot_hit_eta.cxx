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

void roman_pot_hit_eta(const char* fname = "./sim_output/roman_pot_out.root"){

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel
  double degree = TMath::Pi()/180.0;

  TChain* t = new TChain("EVENT");
  t->Add(fname);

  ROOT::RDataFrame d0(*t);
  auto hits_eta = [&](const std::vector<dd4hep::sim::Geant4Tracker::Hit*>& hits){
	  std::vector<double> result;
	  for (const auto& h: hits){
		  result.push_back(h->momentum.eta());
	  }
	  return result;
  };

  auto d1 = d0.Define("hits_eta", hits_eta, {"ForwardRomanPotHits"});

  auto h1 = d1.Histo1D(TH1D("h1", "hits_eta", 300, 0,20), "hits_eta");
  auto n1 = h1->GetMean();
  std::cout << "Pseudorapidity of hits: " << n1 << std::endl;
  TCanvas* c = new TCanvas();
  h1->DrawClone();

  if (n1 < 5) {
	  std::quick_exit(1);
  }

}

