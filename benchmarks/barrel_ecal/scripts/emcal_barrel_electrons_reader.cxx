//////////////////////////
// EMCAL Barrel detector
// Electron dataset
// J.KIM 04/02/2021
//////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include "TH1F.h"
#include "TStyle.h"
#include <iostream>

using namespace HepMC3;

void emcal_barrel_electrons_reader(double e_start = 0.0, double e_end = 30.0, const char* in_fname = "./data/emcal_barrel_electrons.hepmc") {
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

  ReaderAscii hepmc_input(in_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Histograms
  TH1F* h_electrons_energy = new TH1F("h_electron_energy", "electron energy;E [GeV];Events",         100, -0.5, 30.5);
  TH1F* h_electrons_eta    = new TH1F("h_electron_eta",    "electron #eta;#eta;Events",              100, -10.0, 10.0);
  TH1F* h_electrons_theta  = new TH1F("h_electron_theta",  "electron #theta;#theta [degree];Events", 100, -0.5, 180.5);
  TH1F* h_electrons_phi    = new TH1F("h_electron_phi",    "electron #phi;#phi [degree];Events",     100, -180.5, 180.5);
  TH2F* h_electrons_pzpt   = new TH2F("h_electrons_pzpt",  "electron pt vs pz;pt [GeV];pz [GeV]",    100, -0.5, 30.5, 100, -30.5, 30.5);
  TH2F* h_electrons_pxpy   = new TH2F("h_electrons_pxpy",  "electron px vs py;px [GeV];py [GeV]",    100, -30.5, 30.5, 100, -30.5, 30.5);
  TH3F* h_electrons_p      = new TH3F("h_electron_p",      "electron p;px [GeV];py [GeV];pz [GeV]",  100, -30.5, 30.5, 100, -30.5, 30.5, 100, -30.5, 30.5);

  while (!hepmc_input.failed()) {
    // Read event from input file
    hepmc_input.read_event(evt);
    // If reading failed - exit loop
    if (hepmc_input.failed())
      break;

    for (const auto& v : evt.vertices()) {
      for (const auto& p : v->particles_out()) {
        if (p->pid() == 11) {
          h_electrons_energy->Fill(p->momentum().e());
          h_electrons_eta->Fill(p->momentum().eta());
          h_electrons_theta->Fill(p->momentum().theta() * TMath::RadToDeg());
          h_electrons_phi->Fill(p->momentum().phi() * TMath::RadToDeg());
          h_electrons_pzpt->Fill(TMath::Sqrt(p->momentum().px() * p->momentum().px() + p->momentum().py() * p->momentum().py()), p->momentum().pz());
          h_electrons_pxpy->Fill(p->momentum().px(), p->momentum().py());
          h_electrons_p->Fill(p->momentum().px(), p->momentum().py(), p->momentum().pz());
        }
      }
    }
    evt.clear();
    events_parsed++;
  }
  std::cout << "Events parsed and written: " << events_parsed << std::endl;

  TCanvas* c = new TCanvas("c", "c", 500, 500);
  h_electrons_energy->GetYaxis()->SetTitleOffset(1.8);
  h_electrons_energy->SetLineWidth(2);
  h_electrons_energy->SetLineColor(kBlue);
  h_electrons_energy->DrawClone();
  c->SaveAs("results/input_emcal_barrel_electrons_energy.png");
  c->SaveAs("results/input_emcal_barrel_electrons_energy.pdf");

  TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
  h_electrons_eta->GetYaxis()->SetTitleOffset(1.9);
  h_electrons_eta->SetLineWidth(2);
  h_electrons_eta->SetLineColor(kBlue);
  h_electrons_eta->DrawClone();
  c1->SaveAs("results/input_emcal_barrel_electrons_eta.png");
  c1->SaveAs("results/input_emcal_barrel_electrons_eta.pdf");

  TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
  h_electrons_theta->GetYaxis()->SetTitleOffset(1.8);
  h_electrons_theta->SetLineWidth(2);
  h_electrons_theta->SetLineColor(kBlue);
  h_electrons_theta->DrawClone();
  c2->SaveAs("results/input_emcal_barrel_electrons_theta.png");
  c2->SaveAs("results/input_emcal_barrel_electrons_theta.pdf");

  TCanvas* c3 = new TCanvas("c3", "c3", 500, 500);
  h_electrons_phi->GetYaxis()->SetTitleOffset(1.8);
  h_electrons_phi->SetLineWidth(2);
  h_electrons_phi->GetYaxis()->SetRangeUser(0.0, h_electrons_phi->GetMaximum() + 100.0);
  h_electrons_phi->SetLineColor(kBlue);
  h_electrons_phi->DrawClone();
  c3->SaveAs("results/input_emcal_barrel_electrons_phi.png");
  c3->SaveAs("results/input_emcal_barrel_electrons_phi.pdf");

  TCanvas* c4 = new TCanvas("c4", "c4", 500, 500);
  h_electrons_pzpt->GetYaxis()->SetTitleOffset(1.4);
  h_electrons_pzpt->SetLineWidth(2);
  h_electrons_pzpt->SetLineColor(kBlue);
  h_electrons_pzpt->DrawClone("COLZ");
  c4->SaveAs("results/input_emcal_barrel_electrons_pzpt.png");
  c4->SaveAs("results/input_emcal_barrel_electrons_pzpt.pdf");

  TCanvas* c5 = new TCanvas("c5", "c5", 500, 500);
  h_electrons_pxpy->GetYaxis()->SetTitleOffset(1.4);
  h_electrons_pxpy->SetLineWidth(2);
  h_electrons_pxpy->SetLineColor(kBlue);
  h_electrons_pxpy->DrawClone("COLZ");
  c5->SaveAs("results/input_emcal_barrel_electrons_pxpy.png");
  c5->SaveAs("results/input_emcal_barrel_electrons_pxpy.pdf");

  TCanvas* c6 = new TCanvas("c6", "c6", 500, 500);
  h_electrons_p->GetYaxis()->SetTitleOffset(1.8);
  h_electrons_p->GetXaxis()->SetTitleOffset(1.6);
  h_electrons_p->GetZaxis()->SetTitleOffset(1.6);
  h_electrons_p->SetLineWidth(2);
  h_electrons_p->SetLineColor(kBlue);
  h_electrons_p->DrawClone();
  c6->SaveAs("results/input_emcal_barrel_electrons_p.png");
  c6->SaveAs("results/input_emcal_barrel_electrons_p.pdf");
}

