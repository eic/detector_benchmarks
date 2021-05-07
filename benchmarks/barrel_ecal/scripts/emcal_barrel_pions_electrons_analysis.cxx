////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFitResult.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void emcal_barrel_pions_electrons_analysis(const char* input_fname = "sim_output/sim_emcal_barrel_uniform_pions_electrons.root")
{
  //input_fname = "temp_pions_electrons.root";
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


}
