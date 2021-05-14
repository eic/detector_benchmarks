//////////////////////////
// EMCAL Barrel detector
// Pi0 dataset
// M. Scott 05/2021
//////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include "TH1F.h"
#include "TStyle.h"
#include <iostream>

using namespace HepMC3;

void emcal_barrel_pi0_reader(double e_start = 0.0, double e_end = 30.0, const char* in_fname = "./data/emcal_barrel_pi0.hepmc") {
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
  TH1F* h_pi0_energy = new TH1F("h_pi0_energy", "pi0 energy;E [GeV];Events",         100, -0.5, 30.5);
  TH1F* h_pi0_eta    = new TH1F("h_pi0_eta",    "pi0 #eta;#eta;Events",              100, -10.0, 10.0);
  TH1F* h_pi0_theta  = new TH1F("h_pi0_theta",  "pi0 #theta;#theta [degree];Events", 100, -0.5, 180.5);
  TH1F* h_pi0_phi    = new TH1F("h_pi0_phi",    "pi0 #phi;#phi [degree];Events",     100, -180.5, 180.5);
  TH2F* h_pi0_pzpt   = new TH2F("h_pi0_pzpt",   "pi0 pt vs pz;pt [GeV];pz [GeV]",    100, -0.5, 30.5, 100, -30.5, 30.5);
  TH2F* h_pi0_pxpy   = new TH2F("h_pi0_pxpy",   "pi0 px vs py;px [GeV];py [GeV]",    100, -30.5, 30.5, 100, -30.5, 30.5);
  TH3F* h_pi0_p      = new TH3F("h_pi0_p",      "pi0 p;px [GeV];py [GeV];pz [GeV]",  100, -30.5, 30.5, 100, -30.5, 30.5, 100, -30.5, 30.5);

  while (!hepmc_input.failed()) {
    // Read event from input file
    hepmc_input.read_event(evt);
    // If reading failed - exit loop
    if (hepmc_input.failed())
      break;

    for (const auto& v : evt.vertices()) {
      for (const auto& p : v->particles_out()) {
        if (p->pid() == 11) {
          h_pi0_energy->Fill(p->momentum().e());
          h_pi0_eta->Fill(p->momentum().eta());
          h_pi0_theta->Fill(p->momentum().theta() * TMath::RadToDeg());
          h_pi0_phi->Fill(p->momentum().phi() * TMath::RadToDeg());
          h_pi0_pzpt->Fill(TMath::Sqrt(p->momentum().px() * p->momentum().px() + p->momentum().py() * p->momentum().py()), p->momentum().pz());
          h_pi0_pxpy->Fill(p->momentum().px(), p->momentum().py());
          h_pi0_p->Fill(p->momentum().px(), p->momentum().py(), p->momentum().pz());
        }
      }
    }
    evt.clear();
    events_parsed++;
  }
  std::cout << "Events parsed and written: " << events_parsed << std::endl;

  TCanvas* c = new TCanvas("c", "c", 500, 500);
  h_pi0_energy->GetYaxis()->SetTitleOffset(1.8);
  h_pi0_energy->SetLineWidth(2);
  h_pi0_energy->SetLineColor(kBlue);
  h_pi0_energy->DrawClone();
  c->SaveAs("results/input_emcal_barrel_pi0_energy.png");
  c->SaveAs("results/input_emcal_barrel_pi0_energy.pdf");

  TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
  h_pi0_eta->GetYaxis()->SetTitleOffset(1.9);
  h_pi0_eta->SetLineWidth(2);
  h_pi0_eta->SetLineColor(kBlue);
  h_pi0_eta->DrawClone();
  c1->SaveAs("results/input_emcal_barrel_pi0_eta.png");
  c1->SaveAs("results/input_emcal_barrel_pi0_eta.pdf");

  TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
  h_pi0_theta->GetYaxis()->SetTitleOffset(1.8);
  h_pi0_theta->SetLineWidth(2);
  h_pi0_theta->SetLineColor(kBlue);
  h_pi0_theta->DrawClone();
  c2->SaveAs("results/input_emcal_barrel_pi0_theta.png");
  c2->SaveAs("results/input_emcal_barrel_pi0_theta.pdf");

  TCanvas* c3 = new TCanvas("c3", "c3", 500, 500);
  h_pi0_phi->GetYaxis()->SetTitleOffset(1.8);
  h_pi0_phi->SetLineWidth(2);
  h_pi0_phi->GetYaxis()->SetRangeUser(0.0, h_pi0_phi->GetMaximum() + 100.0);
  h_pi0_phi->SetLineColor(kBlue);
  h_pi0_phi->DrawClone();
  c3->SaveAs("results/input_emcal_barrel_pi0_phi.png");
  c3->SaveAs("results/input_emcal_barrel_pi0_phi.pdf");

  TCanvas* c4 = new TCanvas("c4", "c4", 500, 500);
  h_pi0_pzpt->GetYaxis()->SetTitleOffset(1.4);
  h_pi0_pzpt->SetLineWidth(2);
  h_pi0_pzpt->SetLineColor(kBlue);
  h_pi0_pzpt->DrawClone("COLZ");
  c4->SaveAs("results/input_emcal_barrel_pi0_pzpt.png");
  c4->SaveAs("results/input_emcal_barrel_pi0_pzpt.pdf");

  TCanvas* c5 = new TCanvas("c5", "c5", 500, 500);
  h_pi0_pxpy->GetYaxis()->SetTitleOffset(1.4);
  h_pi0_pxpy->SetLineWidth(2);
  h_pi0_pxpy->SetLineColor(kBlue);
  h_pi0_pxpy->DrawClone("COLZ");
  c5->SaveAs("results/input_emcal_barrel_pi0_pxpy.png");
  c5->SaveAs("results/input_emcal_barrel_pi0_pxpy.pdf");

  TCanvas* c6 = new TCanvas("c6", "c6", 500, 500);
  h_pi0_p->GetYaxis()->SetTitleOffset(1.8);
  h_pi0_p->GetXaxis()->SetTitleOffset(1.6);
  h_pi0_p->GetZaxis()->SetTitleOffset(1.6);
  h_pi0_p->SetLineWidth(2);
  h_pi0_p->SetLineColor(kBlue);
  h_pi0_p->DrawClone();
  c6->SaveAs("results/input_emcal_barrel_pi0_p.png");
  c6->SaveAs("results/input_emcal_barrel_pi0_p.pdf");
}
