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

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void zdc_neutrons_analysis(const char* input_fname = "sim_output/sim_zdc_uniform_neutrons.root")
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

  // Thrown Energy [GeV]
  auto Ethr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    auto p = input[2];
    auto energy = TMath::Sqrt(p.psx * p.psx + p.psy * p.psy + p.psz * p.psz + p.mass * p.mass);
    return energy;
  };

  // Number of hits
  auto nhits = [] (const std::vector<dd4pod::CalorimeterHitData>& evt) {return (int) evt.size(); };

  // Energy deposition [GeV]
  auto Esim = [](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt)
      total_edep += i.energyDeposit;
    return total_edep;
  };

  // Sampling fraction = Esampling / Ethrown
  auto fsam = [](const double sampled, const double thrown) {
    return sampled / thrown;
  };

  // Define variables
  auto d1 = d0.Define("Ethr", Ethr, {"mcparticles"})
                .Define("nhits", nhits, {"EcalBarrelHits"})
                .Define("Esim", Esim, {"EcalBarrelHits"})
                .Define("fsam", fsam, {"Esim", "Ethr"});

  // Define Histograms
  auto hEthr = d1.Histo1D(
      {"hEthr", "Thrown Energy; Thrown Energy [GeV]; Events", 100, 0.0, 7.5},
      "Ethr");
  auto hNhits =
      d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events",
                  100, 0.0, 2000.0},
                 "nhits");
  auto hEsim = d1.Histo1D(
      {"hEsim", "Energy Deposit; Energy Deposit [GeV]; Events", 100, 0.0, 1.0},
      "Esim");
  auto hfsam = d1.Histo1D(
      {"hfsam", "Sampling Fraction; Sampling Fraction; Events", 100, 0.0, 0.1},
      "fsam");

  // Event Counts
  auto nevents_thrown = d1.Count();
  std::cout << "Number of Thrown Events: " << (*nevents_thrown) << "\n";

  // Draw Histograms
  {
    TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
    c1->SetLogy(1);
    auto h = hEthr->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c1->SaveAs("results/zdc_neutrons_Ethr.png");
    c1->SaveAs("results/zdc_neutrons_Ethr.pdf");
  }
  std::cout << "derp1\n";

  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c2->SaveAs("results/zdc_neutrons_nhits.png");
    c2->SaveAs("results/zdc_neutrons_nhits.pdf");
  }

  {
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500);
    c3->SetLogy(1);
    auto h = hEsim->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c3->SaveAs("results/zdc_neutrons_Esim.png");
    c3->SaveAs("results/zdc_neutrons_Esim.pdf");
  }

  std::cout << "derp4\n";
  {
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    c4->SetLogy(1);
    auto h = hfsam->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    //h->Fit("gaus", "", "", 0.01, 0.1);
    //h->GetFunction("gaus")->SetLineWidth(2);
    //h->GetFunction("gaus")->SetLineColor(kRed);
    c4->SaveAs("results/zdc_neutrons_fsam.png");
    c4->SaveAs("results/zdc_neutrons_fsam.pdf");
  }
}
