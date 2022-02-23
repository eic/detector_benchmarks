////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include <fstream>
#include <fmt/core.h>

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"

#include "hcal_barrel_common_functions.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void save_canvas(TCanvas* c, std::string label)
{
  c->SaveAs(fmt::format("results/{}.png",label).c_str());
  c->SaveAs(fmt::format("results/{}.pdf",label).c_str());
}

void save_canvas(TCanvas* c, std::string label, std::string particle_label)
{
  std::string label_with_E = fmt::format("hcal_barrel_{}_{}", particle_label, label); 
  save_canvas(c, label_with_E);
}

void hcal_barrel_particles_analysis(std::string particle_name = "electron")
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
  std::string input_fname = fmt::format("sim_output/sim_hcal_barrel_{}.edm4hep.root", particle_name);
  ROOT::RDataFrame d0("events", input_fname);

  // Thrown Energy [GeV]
  auto Ethr = [](std::vector<edm4hep::MCParticleData> const& input) {
    auto p = input[2];
    auto energy = TMath::Sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y + p.momentum.z * p.momentum.z + p.mass * p.mass);
    return energy;
  };

  // Number of hits
  auto nhits = [] (const std::vector<edm4hep::SimCalorimeterHitData>& evt) {return (int) evt.size(); };

  // Energy deposition [GeV]
  auto Esim = [](const std::vector<edm4hep::SimCalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt)
      total_edep += i.energy;
    return total_edep;
  };

  // Sampling fraction = Esampling / Ethrown
  auto fsam = [](const double sampled, const double thrown) {
    return sampled / thrown;
  };

  // Define variables
  auto d1 = d0.Define("Ethr", Ethr, {"MCParticles"})
                .Define("nhits", nhits, {"HcalBarrelHits"})
                .Define("Esim", Esim, {"HcalBarrelHits"})
                .Define("fsam", fsam, {"Esim", "Ethr"});

  // Define Histograms
  auto hEthr = d1.Histo1D(
      {"hEthr", "Thrown Energy; Thrown Energy [GeV]; Events", 100, 0.0, 25.0},
      "Ethr");
  auto hNhits =
      d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events", 100, 0.0, 2000.0},
      "nhits");
  auto hEsim = d1.Histo1D(
      {"hEsim", "Energy Deposit; Energy Deposit [GeV]; Events", 500, 0.0, 0.5},
      "Esim");
  auto hfsam = d1.Histo1D(
      {"hfsam", "Sampling Fraction; Sampling Fraction; Events", 200, 0.0, 0.05},
      "fsam");

  // Event Counts
  auto nevents_thrown = d1.Count();
  std::cout << "Number of Thrown Events: " << (*nevents_thrown) << "\n";

  // Draw Histograms
  {
    TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
    c1->SetLogy(1);
    auto h = hEthr->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    save_canvas(c1,"Ethr",particle_name);
  }

  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    save_canvas(c2,"nhits",particle_name);
  }

  {
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500);
    c3->SetLogy(1);
    auto h = hEsim->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    double up_fit = h->GetMean() + 5*h->GetStdDev();
    double down_fit = h->GetMean() - 5*h->GetStdDev();
    h->GetXaxis()->SetRangeUser(0.,up_fit);
    save_canvas(c3,"Esim",particle_name);
  }

  {
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    c4->SetLogy(1);
    auto h = hfsam->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    double up_fit = h->GetMean() + 5*h->GetStdDev();
    double down_fit = h->GetMean() - 5*h->GetStdDev();
    h->Fit("gaus", "", "", down_fit, up_fit);
    h->GetXaxis()->SetRangeUser(0.,up_fit);
    TF1 *gaus = h->GetFunction("gaus");
    gaus->SetLineWidth(2);
    gaus->SetLineColor(kRed);
    save_canvas(c4,"fsam",particle_name);
  }
}

