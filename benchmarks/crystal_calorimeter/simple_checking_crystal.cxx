R__LOAD_LIBRARY(libDDG4IO.so)
#include "DDG4/Geant4Data.h"
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include <random>

void simple_checking_crystal(const char* fname = "sim_output/output_emcal_electrons.edm4hep.root"){
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel
  double degree = TMath::Pi()/180.0;

  TChain* t = new TChain("events");
  t->Add(fname);

  ROOT::RDataFrame d0(*t);//, {"EcalHits","MCParticles"});
  
  auto nhits = [] (const std::vector<dd4pod::CalorimeterHit>& hits){ return (int) hits.size(); };

  auto d1 = d0.Define("nhits", nhits, {"CrystalEcalHits"});
  auto h0 = d1.Histo1D(TH1D("h0", "nhits; ", 20, 0,20), "nhits");

  auto n0 = d1.Filter([](int n){ return (n>0); },{"nhits"}).Count();

  TCanvas* c = new TCanvas();

  std::cout << *n0 << " events with nonzero hits\n";
  if(*n0<5) {
    std::quick_exit(1);
  }

}

