////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"

#include "benchmark.h"
#include "mt.h"
#include "util.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFitResult.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void emcal_barrel_pions_analysis(const char* input_fname = "sim_output/sim_emcal_barrel_uniform_pions.root")
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

  //Tests
  std::string test_tag = "Barrel_emcal_pi0";
  //TODO: Change test_tag to something else
  std:string detector = "Barrel_emcal";
  // Energy resolution in the barrel region(-1 < eta < 1)
  // Taken from : Initial considerations for EMCal of the EIC detector by A. Bazilevsky
  // sigma_E / E = 12% / E^0.5 convoluted with 2%
  // sigma_E / E = [ (0.12/E^0.5)^2 + 0.02^2]^0.5, with E in [GeV]
  double thrown_energy = 5; // Current thrown energy, will need to grab from json file
  double resolutionTarget = TMath::Sqrt(0.12 * 0.12 / thrown_energy + 0.02 * 0.02);

  eic::util::Test pi0_energy_resolution{
      {{"name", fmt::format("{}_energy_resolution", test_tag)},
       {"title", "Pion0 Energy resolution"},
       {"description",
        fmt::format("Pion0 energy resolution with {}, estimated using a Gaussian fit.", detector)},
       {"quantity", "resolution (in %)"},
       {"target", std::to_string(resolutionTarget)}}};


  ROOT::EnableImplicitMT();
  ROOT::RDataFrame d0("events", input_fname);

  // Sampling Fraction
  double samp_frac = 0.0136;

  // Thrown Energy [GeV]
  auto Ethr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    std::vector<double> result;
    result.push_back(TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz + input[2].mass*input[2].mass));
  return result;
  };

  // Number of hits
  auto nhits = [] (const std::vector<dd4pod::CalorimeterHitData>& evt) {return (int) evt.size(); };

  // Energy deposition [GeV]
  auto Esim = [](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    std::vector<double> result;
    auto total_edep = 0.0;
    for (const auto& i: evt)
      total_edep += i.energyDeposit;
    result.push_back(total_edep);
  return result;
  };

  // Sampling fraction = Esampling / Ethrown
  auto fsam = [](const std::vector<double>& sampled, const std::vector<double>& thrown) {
    std::vector<double> result;
    auto it_sam = sampled.cbegin();
    auto it_thr = thrown.cbegin();
    for (; it_sam != sampled.end() && it_thr != thrown.end(); ++it_sam, ++it_thr) {
        result.push_back(*it_sam / *it_thr);
    }
    return result;
  };

  // Energy Resolution = Esampling/Sampling_fraction - Ethrown
  auto eResol = [samp_frac](const std::vector<double>& sampled, const std::vector<double>& thrown) {
    std::vector<double> result;
    auto it_sam = sampled.cbegin();
    auto it_thr = thrown.cbegin();
    for (; it_sam != sampled.end() && it_thr != thrown.end(); ++it_sam, ++it_thr) {
        result.push_back(*it_sam / samp_frac - *it_thr);
    }
    return result;
  };

  // Relative Energy Resolution = (Esampling/Sampling fraction - Ethrown)/Ethrown
  auto eResol_rel = [samp_frac](const std::vector<double>& sampled, const std::vector<double>& thrown) {
    std::vector<double> result;
    auto it_sam = sampled.cbegin();
    auto it_thr = thrown.cbegin();
    for (; it_sam != sampled.end() && it_thr != thrown.end(); ++it_sam, ++it_thr) {
        result.push_back((*it_sam / samp_frac - *it_thr) / *it_thr);
    }
    return result;
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
              .Define("Esim",   Esim,       {"EcalBarrelHits"})
              .Define("fsam",   fsam,       {"Esim","Ethr"})
              .Define("pid",    getpid,     {"mcparticles"})
              .Define("dau",    getdau,     {"mcparticles"})
              ;

  // Define Histograms
  auto hEthr  = d1.Histo1D({"hEthr",  "Thrown Energy; Thrown Energy [GeV]; Events",        100,  0.0,    7.5}, "Ethr");
  auto hNhits = d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events", 100,  0.0, 2000.0}, "nhits");
  auto hEsim  = d1.Histo1D({"hEsim",  "Energy Deposit; Energy Deposit [GeV]; Events",      100,  0.0,    1.0}, "Esim");
  auto hfsam  = d1.Histo1D({"hfsam",  "Sampling Fraction; Sampling Fraction; Events",      100,  0.0,    0.1}, "fsam");
  auto hpid   = d1.Histo1D({"hpid",   "PID; PID; Count",                                   100,  -220,   220}, "pid");
  auto hdau   = d1.Histo1D({"hdau",   "Number of Daughters; Number of Daughters; Count",   10,   0,      10},  "dau");

  // Set sampling Fraction, ideally this will be taken from a json file 
  samp_frac = hfsam -> GetMean();

  auto d2 = d1.Define("dE",     eResol,     {"Esim","Ethr"})
              .Define("dE_rel", eResol_rel, {"Esim","Ethr"})
              ;

  // Event Counts
  auto nevents_thrown      = d1.Count();
  std::cout << "Number of Thrown Events: " << (*nevents_thrown) << "\n";

  // Draw Histograms
  TCanvas *c1 = new TCanvas("c1", "c1", 700, 500);
  c1->SetLogy(1);
  hEthr->GetYaxis()->SetTitleOffset(1.4);
  hEthr->SetLineWidth(2);
  hEthr->SetLineColor(kBlue);
  hEthr->DrawClone();
  c1->SaveAs("results/emcal_barrel_pions_Ethr.png");
  c1->SaveAs("results/emcal_barrel_pions_Ethr.pdf");

  TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);
  c2->SetLogy(1);
  hNhits->GetYaxis()->SetTitleOffset(1.4);
  hNhits->SetLineWidth(2);
  hNhits->SetLineColor(kBlue);
  hNhits->DrawClone();
  c2->SaveAs("results/emcal_barrel_pions_nhits.png");
  c2->SaveAs("results/emcal_barrel_pions_nhits.pdf");

  TCanvas *c3 = new TCanvas("c3", "c3", 700, 500);
  c3->SetLogy(1);
  hEsim->GetYaxis()->SetTitleOffset(1.4);
  hEsim->SetLineWidth(2);
  hEsim->SetLineColor(kBlue);
  hEsim->DrawClone();
  c3->SaveAs("results/emcal_barrel_pions_Esim.png"); 
  c3->SaveAs("results/emcal_barrel_pions_Esim.pdf");

  TCanvas *c4 = new TCanvas("c4", "c4", 700, 500);
  c4->SetLogy(1);
  hfsam->GetYaxis()->SetTitleOffset(1.4);
  hfsam->SetLineWidth(2);
  hfsam->SetLineColor(kBlue);
  hfsam->Fit("gaus","","",0.005,0.1);
  hfsam->GetFunction("gaus")->SetLineWidth(2);
  hfsam->GetFunction("gaus")->SetLineColor(kRed);
  hfsam->DrawClone();
  c4->SaveAs("results/emcal_barrel_pions_fsam.png");
  c4->SaveAs("results/emcal_barrel_pions_fsam.pdf");

  TCanvas *c5 = new TCanvas("c5", "c5", 700, 500);
  c5->SetLogy(1);
  hpid->GetYaxis()->SetTitleOffset(1.4);
  hpid->SetLineWidth(2);
  hpid->SetLineColor(kBlue);
  hpid->DrawClone();
  c5->SaveAs("results/emcal_barrel_pions_pid.png");
  c5->SaveAs("results/emcal_barrel_pions_pid.pdf");

  TCanvas *c6 = new TCanvas("c6", "c6", 700, 500);
  c5->SetLogy(1);
  hdau->GetYaxis()->SetTitleOffset(1.4);
  hdau->SetLineWidth(2);
  hdau->SetLineColor(kBlue);
  hdau->DrawClone();
  c6->SaveAs("results/emcal_barrel_pions_dau.png");
  c6->SaveAs("results/emcal_barrel_pions_dau.pdf");

  //Energy Resolution Calculation
  auto hdE          = d2.Histo1D({"hdE",      "dE; dE[GeV]; Events",              100, -3.0, 3.0}, "dE");
  auto hdE_rel      = d2.Histo1D({"hdE_rel",  "dE Relative; dE Relative; Events", 100, -3.0, 3.0}, "dE_rel");
  hdE->Fit("gaus", "", "", -3.0,  3.0);
  double* res       = hdE->GetFunction("gaus")->GetParameters();
  double sigmaOverE = res[2] / thrown_energy;

  //Pass/Fail
  if (sigmaOverE <= resolutionTarget) {
    pi0_energy_resolution.pass(sigmaOverE);
  } else {
    pi0_energy_resolution.fail(sigmaOverE);
  }
  //std::printf("Energy Resolution is %f\n", res[2]);

  //Energy Resolution Histogram Plotting
  auto *cdE = new TCanvas("cdE", "cdE", 700, 500);
  cdE->SetLogy(1);
  hdE->GetYaxis()->SetTitleOffset(1.4);
  hdE->SetLineWidth(2);
  hdE->SetLineColor(kBlue);
  hdE->GetFunction("gaus")->SetLineWidth(2);
  hdE->GetFunction("gaus")->SetLineColor(kRed);
  hdE->DrawClone();
  cdE->SaveAs("results/emcal_barrel_pi0_dE.png");
  cdE->SaveAs("results/emcal_barrel_pi0_dE.pdf");

  auto *cdE_rel = new TCanvas("cdE_rel", "cdE_rel", 700, 500);
  hdE_rel->GetYaxis()->SetTitleOffset(1.4);
  hdE_rel->SetLineWidth(2);
  hdE_rel->SetLineColor(kBlue);
  hdE_rel->DrawClone();
  cdE_rel->SaveAs("results/emcal_barrel_pi0_dE_rel.png");
  cdE_rel->SaveAs("results/emcal_barrel_pi0_dE_rel.pdf");

  eic::util::write_test({pi0_energy_resolution}, fmt::format("{}_pions.json", detector));


}
