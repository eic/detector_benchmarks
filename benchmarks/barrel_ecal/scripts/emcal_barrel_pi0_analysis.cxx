////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"

R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFitResult.h"
#include <nlohmann/json.hpp>

using ROOT::RDataFrame;
using namespace ROOT::VecOps;
using json = nlohmann::json;

void emcal_barrel_pi0_analysis(
                                const char* input_fname = "sim_output/sim_emcal_barrel_pi0.edm4hep.root"
                                //const char* input_fname = "../sim_output/sim_emcal_barrel_uniform_pi0.edm4hep.root"
                                )
{
  // Setting for graphs
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  gStyle->SetLineWidth(2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.14);

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d0("events", input_fname);

  // Sampling Fraction grabbed from json file
  // Note that this value is derived from electrons
  json j;
  std::ifstream prev_steps_ifstream("results/emcal_barrel_calibration.json");
  prev_steps_ifstream >> j;
  double samp_frac = j["electron"]["sampling_fraction"];

  // Thrown Energy [GeV]
  auto Ethr = [](std::vector<edm4hep::MCParticleData> const& input) {
    return TMath::Sqrt(input[2].momentum.x*input[2].momentum.x + input[2].momentum.y*input[2].momentum.y + input[2].momentum.z*input[2].momentum.z + input[2].mass*input[2].mass);
  };

  // Number of hits
  auto nhits = [] (const std::vector<edm4hep::SimCalorimeterHitData>& evt) {return (int) evt.size(); };

  // Energy deposition [GeV]
  auto Esim = [](const std::vector<edm4hep::SimCalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt){
      total_edep += i.energy;
    }
    return total_edep;
  };

  // Sampling fraction = Esampling / Ethrown
  auto fsam = [](const double& sampled, const double& thrown) {
    return sampled / thrown;
  };

  // Energy Resolution = Esampling/Sampling_fraction - Ethrown
  auto eResol = [&](const double& sampled, const double& thrown){
    return sampled / samp_frac - thrown;
  };

  // Relative Energy Resolution = (Esampling/Sampling fraction - Ethrown)/Ethrown
  auto eResol_rel = [&](const double& sampled, const double& thrown){
    return (sampled / samp_frac - thrown) / thrown;
  };

  // Returns the pdgID of the particle
  auto getpid = [](std::vector<edm4hep::MCParticleData> const& input) {
    return input[2].PDG;
  };

  // Returns number of particle daughters
  auto getdau = [](std::vector<edm4hep::MCParticleData> const& input) {
    return input[2].daughters_begin;
  };

  // Define variables
  auto d1 = ROOT::RDF::RNode(
    d0.Define("Ethr", Ethr, {"MCParticles"})
      .Define("pid", getpid, {"MCParticles"})
      .Define("dau", getdau, {"MCParticles"})
      .Define("nhits", nhits, {"EcalBarrelHits"})
  );

  auto Ethr_max = 7.5;
  auto fsam_est = 1.0;
  if (d1.HasColumn("EcalBarrelScFiHits")) {
    d1 = d1.Define("EsimImg", Esim, {"EcalBarrelHits"})
           .Define("EsimScFi", Esim, {"EcalBarrelScFiHits"})
           .Define("Esim", "EsimImg+EsimScFi")
           .Define("fsamImg", fsam, {"EsimImg", "Ethr"})
           .Define("fsamScFi", fsam, {"EsimScFi", "Ethr"})
           .Define("fsam", fsam, {"Esim", "Ethr"});
    fsam_est = 0.1;
  } else {
    d1 = d1.Define("Esim", Esim, {"EcalBarrelHits"})
           .Define("fsam", fsam, {"Esim", "Ethr"});
    fsam_est = 1.0;
  }
  d1 = d1.Define("dE",         eResol,        {"Esim","Ethr"})
         .Define("dE_rel",     eResol_rel,    {"Esim","Ethr"});

  // Define Histograms
  std::vector <std::string> titleStr = {
                                        "Thrown Energy; Thrown Energy [GeV]; Events",
                                        "Number of hits per events; Number of hits; Events",
                                        "Energy Deposit; Energy Deposit [GeV]; Events",
                                        "dE Relative; dE Relative; Events"
  };

  std::vector<std::vector<double>> range = {{0, Ethr_max}, {0, 2000}, {0, fsam_est * Ethr_max}, {-3, 3}};
  std::vector<std::string> col           = {"Ethr",   "nhits",   "Esim", "dE_rel"}; 

  double meanE  = 5; 
  int nCol = range.size();
  for (int i = 0; i < nCol; i++){
    int binNum = 100;
    auto h = d1.Histo1D({"hist", titleStr[i].c_str(), binNum, range[i][0], range[i][1]}, col[i].c_str());
    if (col[i] == "Ethr"){
      meanE = h->GetMean();
    }
    auto *c = new TCanvas("c", "c", 700, 500);
    c->SetLogy(1);
    auto h1 = h->DrawCopy();
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->SetLineWidth(2);
    h1->SetLineColor(kBlue);
    c->SaveAs((fmt::format("results/emcal_barrel_pi0_{}.png", col[i])).c_str());
    c->SaveAs((fmt::format("results/emcal_barrel_pi0_{}.pdf", col[i])).c_str());
    //std::printf("Generic %d\n", i);
  }

  // Resolution Plots
  titleStr = {
              "Sampling Fraction; Sampling Fraction; Events",
              "dE; dE[GeV]; Events"
  };
  range = {{0,fsam_est}, {-3, 3}};
  col   = {"fsam",   "dE"}; 
  nCol  = range.size();
  std::printf("Here %d\n", 10);
  std::vector<std::vector<double>> fitRange = {{0.005, fsam_est}, {-3, 3}};
  double sigmaOverE = 0;

  auto hr = d1.Histo1D({"histr", titleStr[0].c_str(), 150, range[0][0], range[0][1]}, col[0].c_str());
  auto *c = new TCanvas("c", "c", 700, 500);
  c->SetLogy(1);
  auto h2 = hr->DrawCopy();
  h2->GetYaxis()->SetTitleOffset(1.4);
  h2->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->Fit("gaus","","", fitRange[0][0], fitRange[0][1]);
  h2->GetFunction("gaus")->SetLineWidth(2);
  h2->GetFunction("gaus")->SetLineColor(kRed);
  
  c->SaveAs((fmt::format("results/emcal_barrel_pi0_{}.png", col[0])).c_str());
  c->SaveAs((fmt::format("results/emcal_barrel_pi0_{}.pdf", col[0])).c_str());
  std::printf("Resolution %d\n", 0);


  auto hs = d1.Histo1D({"hists", titleStr[1].c_str(), 100, range[1][0], range[1][1]}, col[1].c_str());
  auto *c1 = new TCanvas("c1", "c1", 700, 500);
  c1->SetLogy(1);
  auto h3 = hs->DrawCopy();
  h3->GetYaxis()->SetTitleOffset(1.4);
  h3->SetLineWidth(2);
  h3->SetLineColor(kBlue);
  auto fit = h3->Fit("gaus","","", fitRange[1][0], fitRange[1][1]);
  double* res = h3->GetFunction("gaus")->GetParameters();
  sigmaOverE  = res[2] / meanE;
  
  c1->SaveAs((fmt::format("results/emcal_barrel_pi0_{}.png", col[1])).c_str());
  c1->SaveAs((fmt::format("results/emcal_barrel_pi0_{}.pdf", col[1])).c_str());
  std::printf("Resolution %d\n", 1);
 
  // Energy Resolution Calculation
  std::string test_tag = "Barrel_emcal_pi0";// TODO: Change test_tag to something else
  std:string detEle    = "Barrel_emcal";

  // Energy resolution in the barrel region (-1 < eta < 1)
  // Taken from : Initial considerations for EMCal of the EIC detector by A. Bazilevsky
  // sigma_E / E = 12% / E^0.5 convoluted with 2%
  // sigma_E / E = [ (0.12/E^0.5)^2 + 0.02^2]^0.5, with E in [GeV]
  double resolutionTarget = TMath::Sqrt(0.12 * 0.12 / meanE + 0.02 * 0.02);

  common_bench::Test pi0_energy_resolution{
   {{"name", fmt::format("{}_energy_resolution", test_tag)},
   {"title", "Pi0 Energy resolution"},
   {"description",
    fmt::format("Pi0 energy resolution for {}, estimated using a Gaussian fit.", detEle)},
   {"quantity", "resolution (in %)"},
   {"target", std::to_string(resolutionTarget)}}
  };

  //// Pass/Fail
  sigmaOverE <= resolutionTarget ? pi0_energy_resolution.pass(sigmaOverE) : pi0_energy_resolution.fail(sigmaOverE);
  common_bench::write_test({pi0_energy_resolution}, fmt::format("results/{}_pi0.json", detEle));
}
