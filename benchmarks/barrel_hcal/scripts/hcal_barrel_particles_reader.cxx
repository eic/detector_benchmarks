//////////////////////////
// HCAL Barrel detector
// J.KIM 04/02/2021
// M.Zurek 05/05/2021
// W.Deconinck 05/28/21
//////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iostream>
#include <fstream>
#include <fmt/core.h>

#include "hcal_barrel_common_functions.h"

using namespace HepMC3;

void save_canvas(TCanvas* c, std::string label)
{
  c->SaveAs(fmt::format("results/{}.png",label).c_str());
  c->SaveAs(fmt::format("results/{}.pdf",label).c_str());
}

void save_canvas(TCanvas* c, std::string label, std::string particle_label)
{
  std::string label_with_E = fmt::format("input_hcal_barrel_{}_{}", particle_label, label); 
  save_canvas(c, label_with_E);
}

void hcal_barrel_particles_reader(std::string particle_name = "electron")
{
  // Setting for graphs
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
  gStyle->SetLineWidth(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.17);

  std::string in_fname = fmt::format("./data/hcal_barrel_{}.hepmc",particle_name);
  ReaderAscii hepmc_input(in_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Histograms
  TH1F* h_energy = new TH1F(fmt::format("h_{}_energy",particle_name).c_str(), fmt::format("{} energy;E [GeV];Events",particle_name).c_str(),         100, -0.5, 30.5);
  TH1F* h_eta    = new TH1F(fmt::format("h_{}_eta",particle_name).c_str(),    fmt::format("{} #eta;#eta;Events",particle_name).c_str(),              100, -10.0, 10.0);
  TH1F* h_theta  = new TH1F(fmt::format("h_{}_theta",particle_name).c_str(),  fmt::format("{} #theta;#theta [degree];Events",particle_name).c_str(), 100, -0.5, 180.5);
  TH1F* h_phi    = new TH1F(fmt::format("h_{}_phi",particle_name).c_str(),    fmt::format("{} #phi;#phi [degree];Events",particle_name).c_str(),     100, -180.5, 180.5);
  TH2F* h_pzpt   = new TH2F(fmt::format("h_{}_pzpt",particle_name).c_str(),  fmt::format("{} pt vs pz;pt [GeV];pz [GeV]",particle_name).c_str(),    100, -0.5, 30.5, 100, -30.5, 30.5);
  TH2F* h_pxpy   = new TH2F(fmt::format("h_{}_pxpy",particle_name).c_str(),  fmt::format("{} px vs py;px [GeV];py [GeV]",particle_name).c_str(),    100, -30.5, 30.5, 100, -30.5, 30.5);
  TH3F* h_p      = new TH3F(fmt::format("h_{}_p",particle_name).c_str(),     fmt::format("{} p;px [GeV];py [GeV];pz [GeV]",particle_name).c_str(),  100, -30.5, 30.5, 100, -30.5, 30.5, 100, -30.5, 30.5);

  while (!hepmc_input.failed()) {
    // Read event from input file
    hepmc_input.read_event(evt);
    // If reading failed - exit loop
    if (hepmc_input.failed())
      break;

    auto [id, mass] = extract_particle_parameters(particle_name);

    for (const auto& v : evt.vertices()) {
      for (const auto& p : v->particles_out()) {
        if (p->pid() == id) {
          h_energy->Fill(p->momentum().e());
          h_eta->Fill(p->momentum().eta());
          h_theta->Fill(p->momentum().theta() * TMath::RadToDeg());
          h_phi->Fill(p->momentum().phi() * TMath::RadToDeg());
          h_pzpt->Fill(TMath::Sqrt(p->momentum().px() * p->momentum().px() + p->momentum().py() * p->momentum().py()), p->momentum().pz());
          h_pxpy->Fill(p->momentum().px(), p->momentum().py());
          h_p->Fill(p->momentum().px(), p->momentum().py(), p->momentum().pz());
        }
      }
    }
    evt.clear();
    events_parsed++;
  }
  std::cout << "Events parsed and written: " << events_parsed << std::endl;

  TCanvas* c = new TCanvas("c", "c", 500, 500);
  h_energy->GetYaxis()->SetTitleOffset(1.8);
  h_energy->SetLineWidth(2);
  h_energy->SetLineColor(kBlue);
  h_energy->DrawClone();
  save_canvas(c, "energy", particle_name);

  TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
  h_eta->GetYaxis()->SetTitleOffset(1.9);
  h_eta->SetLineWidth(2);
  h_eta->SetLineColor(kBlue);
  h_eta->DrawClone();
  save_canvas(c1, "eta", particle_name);

  TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
  h_theta->GetYaxis()->SetTitleOffset(1.8);
  h_theta->SetLineWidth(2);
  h_theta->SetLineColor(kBlue);
  h_theta->DrawClone();
  save_canvas(c2, "theta", particle_name);

  TCanvas* c3 = new TCanvas("c3", "c3", 500, 500);
  h_phi->GetYaxis()->SetTitleOffset(1.8);
  h_phi->SetLineWidth(2);
  h_phi->GetYaxis()->SetRangeUser(0.0, h_phi->GetMaximum() + 100.0);
  h_phi->SetLineColor(kBlue);
  h_phi->DrawClone();
  save_canvas(c3, "phi", particle_name);

  TCanvas* c4 = new TCanvas("c4", "c4", 500, 500);
  h_pzpt->GetYaxis()->SetTitleOffset(1.4);
  h_pzpt->SetLineWidth(2);
  h_pzpt->SetLineColor(kBlue);
  h_pzpt->DrawClone("COLZ");
  save_canvas(c4, "pzpt", particle_name);

  TCanvas* c5 = new TCanvas("c5", "c5", 500, 500);
  h_pxpy->GetYaxis()->SetTitleOffset(1.4);
  h_pxpy->SetLineWidth(2);
  h_pxpy->SetLineColor(kBlue);
  h_pxpy->DrawClone("COLZ");
  save_canvas(c5, "pxpy", particle_name);

  TCanvas* c6 = new TCanvas("c6", "c6", 500, 500);
  h_p->GetYaxis()->SetTitleOffset(1.8);
  h_p->GetXaxis()->SetTitleOffset(1.6);
  h_p->GetZaxis()->SetTitleOffset(1.6);
  h_p->SetLineWidth(2);
  h_p->SetLineColor(kBlue);
  h_p->DrawClone();
  save_canvas(c6, "p", particle_name);
}

