////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"

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

void emcal_barrel_pi0_analysis(const char* input_fname = "sim_output/sim_emcal_barrel_pi0.root")
{
  //input_fname = "../sim_output/sim_emcal_barrel_pi0.root";
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
  //double meanE = 5; // Calculated later
  auto Ethr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz + input[2].mass*input[2].mass);
  };

  // Number of hits
  auto nhits = [] (const std::vector<dd4pod::CalorimeterHitData>& evt) {return (int) evt.size(); };

  // Energy deposition [GeV]
  auto Esim = [](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt){
      total_edep += i.energyDeposit;
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
  auto getpid = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return input[2].pdgID;
  };

  // Returns number of particle daughters
  auto getdau = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return input[2].daughters_begin;
  };

  // Define variables
  auto d1 = d0.Define("Ethr",   Ethr,       {"mcparticles"})
              .Define("nhits",  nhits,      {"EcalBarrelHits"})
              .Define("EsimImg", Esim,      {"EcalBarrelHits"})
              .Define("EsimScFi", Esim,      {"EcalBarrelScFiHits"})
              .Define("Esim", "EsimImg+EsimScFi")
              .Define("fsam",   fsam,       {"Esim","Ethr"})
              .Define("pid",    getpid,     {"mcparticles"})
              .Define("dau",    getdau,     {"mcparticles"})
              .Define("dE",     eResol,     {"Esim","Ethr"})
              .Define("dE_rel", eResol_rel, {"Esim","Ethr"})
              ;

  // Define Histograms
  auto hEthr  = d1.Histo1D({"hEthr",  "Thrown Energy; Thrown Energy [GeV]; Events",        100,  0.0,    7.5}, "Ethr");
  auto hNhits = d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events", 100,  0.0, 2000.0}, "nhits");
  auto hEsim  = d1.Histo1D({"hEsim",  "Energy Deposit; Energy Deposit [GeV]; Events",      100,  0.0,    2.0}, "Esim");
  auto hfsam  = d1.Histo1D({"hfsam",  "Sampling Fraction; Sampling Fraction; Events",      150,  0.0,    0.15}, "fsam");
  auto hpid   = d1.Histo1D({"hpid",   "PID; PID; Count",                                   100,  -220,   220}, "pid");
  auto hdau   = d1.Histo1D({"hdau",   "Number of Daughters; Number of Daughters; Count",   10,   0,      10},  "dau");

  // put this line above GetMean below so it can be lazily evaluated with the histgorams.
  auto nevents_thrown = d1.Count();

  const double meanE     = hEthr->GetMean(); 

  // Event Counts
  std::cout << "Number of Thrown Events: " << (*nevents_thrown) << "\n";

  // Draw Histograms
  TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
  c1->SetLogy(1);
  fmt::print("0\n");
  auto h1 = hEthr->DrawCopy();
  //h1->GetYaxis()->SetTitleOffset(1.4);
  h1->SetLineWidth(2);
  h1->SetLineColor(kBlue);
  c1->SaveAs("results/emcal_barrel_pi0_Ethr.png");
  c1->SaveAs("results/emcal_barrel_pi0_Ethr.pdf");

  fmt::print("1\n");

  //TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);
  //c2->SetLogy(1);
  //h1 = hNhits->DrawCopy();
  //h1->GetYaxis()->SetTitleOffset(1.4);
  //h1->SetLineWidth(2);
  //h1->SetLineColor(kBlue);
  //c2->SaveAs("results/emcal_barrel_pi0_nhits.png");
  //c2->SaveAs("results/emcal_barrel_pi0_nhits.pdf");

  //TCanvas *c3 = new TCanvas("c3", "c3", 700, 500);
  //c3->SetLogy(1);
  //h1 = hEsim->DrawCopy();
  //h1->GetYaxis()->SetTitleOffset(1.4);
  //h1->SetLineWidth(2);
  //h1->SetLineColor(kBlue);
  //c3->SaveAs("results/emcal_barrel_pi0_Esim.png"); 
  //c3->SaveAs("results/emcal_barrel_pi0_Esim.pdf");

  //TCanvas *c4 = new TCanvas("c4", "c4", 700, 500);
  //c4->SetLogy(1);
  //h1 = hfsam->DrawCopy();
  //h1->GetYaxis()->SetTitleOffset(1.4);
  //h1->SetLineWidth(2);
  //h1->SetLineColor(kBlue);
  //h1->Fit("gaus","","",0.005,0.1);
  //h1->GetFunction("gaus")->SetLineWidth(2);
  //h1->GetFunction("gaus")->SetLineColor(kRed);
  //c4->SaveAs("results/emcal_barrel_pi0_fsam.png");
  //c4->SaveAs("results/emcal_barrel_pi0_fsam.pdf");

  //TCanvas *c5 = new TCanvas("c5", "c5", 700, 500);
  //c5->SetLogy(1);
  //h1 = hpid->DrawCopy();
  //h1->GetYaxis()->SetTitleOffset(1.4);
  //h1->SetLineWidth(2);
  //h1->SetLineColor(kBlue);
  //c5->SaveAs("results/emcal_barrel_pi0_pid.png");
  //c5->SaveAs("results/emcal_barrel_pi0_pid.pdf");

  //TCanvas *c6 = new TCanvas("c6", "c6", 700, 500);
  //c5->SetLogy(1);
  //h1 = hdau->DrawCopy();
  //h1->GetYaxis()->SetTitleOffset(1.4);
  //h1->SetLineWidth(2);
  //h1->SetLineColor(kBlue);
  //c6->SaveAs("results/emcal_barrel_pi0_dau.png");
  //c6->SaveAs("results/emcal_barrel_pi0_dau.pdf");

  //// Energy Resolution Calculation
  //std::string test_tag = "Barrel_emcal_pi0";// TODO: Change test_tag to something else
  //std:string detEle    = "Barrel_emcal";

  //// Energy resolution in the barrel region (-1 < eta < 1)
  //// Taken from : Initial considerations for EMCal of the EIC detector by A. Bazilevsky
  //// sigma_E / E = 12% / E^0.5 convoluted with 2%
  //// sigma_E / E = [ (0.12/E^0.5)^2 + 0.02^2]^0.5, with E in [GeV]
  //
  //double resolutionTarget = TMath::Sqrt(0.12 * 0.12 / meanE + 0.02 * 0.02);

  //common_bench::Test pi0_energy_resolution{
  // {{"name", fmt::format("{}_energy_resolution", test_tag)},
  // {"title", "Pi0 Energy resolution"},
  // {"description",
  //  fmt::format("Pi0 energy resolution for {}, estimated using a Gaussian fit.", detEle)},
  // {"quantity", "resolution (in %)"},
  // {"target", std::to_string(resolutionTarget)}}
  //};

  //// Histograms and Fitting
  //auto hdE          = d1.Histo1D({"hdE",      "dE; dE[GeV]; Events",              100, -3.0, 3.0}, "dE");
  //auto hdE_rel      = d1.Histo1D({"hdE_rel",  "dE Relative; dE Relative; Events", 100, -3.0, 3.0}, "dE_rel");
  //auto hdEcopy      = hdE->DrawCopy();
  //hdEcopy->Fit("gaus", "", "", -3.0,  3.0);
  //double* res       = hdEcopy->GetFunction("gaus")->GetParameters();
  //double sigmaOverE = res[2] / meanE;

  //// Pass/Fail
  //sigmaOverE <= resolutionTarget ? pi0_energy_resolution.pass(sigmaOverE) : pi0_energy_resolution.fail(sigmaOverE);
  //std::printf("Energy Resolution is %f\n", res[2]);

  //// Energy Resolution Histogram Plotting
  //auto *cdE = new TCanvas("cdE", "cdE", 700, 500);
  //cdE->SetLogy(1);
  //h1 = hdEcopy->DrawCopy();
  //hdEcopy->GetYaxis()->SetTitleOffset(1.4);
  //hdEcopy->SetLineWidth(2);
  //hdEcopy->SetLineColor(kBlue);
  //hdEcopy->GetFunction("gaus")->SetLineWidth(2);
  //hdEcopy->GetFunction("gaus")->SetLineColor(kRed);
  //cdE->SaveAs("results/emcal_barrel_pi0_dE.png");
  //cdE->SaveAs("results/emcal_barrel_pi0_dE.pdf");

  //auto *cdE_rel = new TCanvas("cdE_rel", "cdE_rel", 700, 500);
  //h1 = hdE_rel->DrawCopy();
  //hdE_rel->GetYaxis()->SetTitleOffset(1.4);
  //hdE_rel->SetLineWidth(2);
  //hdE_rel->SetLineColor(kBlue);
  //cdE_rel->SaveAs("results/emcal_barrel_pi0_dE_rel.png");
  //cdE_rel->SaveAs("results/emcal_barrel_pi0_dE_rel.pdf");

  //common_bench::write_test({pi0_energy_resolution}, fmt::format("results/{}_pi0.json", detEle));
}
