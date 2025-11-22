////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include <fmt/core.h>

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/PhotoMultiplierHitCollection.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void mrich_analysis(const char* input_fname = "sim_output/sim_pid_backward_e-_5GeV.edm4hep.root", const char* input_pname = "e-")
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

  ROOT::EnableImplicitMT(4); // use 4 threads
  ROOT::RDataFrame d0("events", input_fname);

  // Number of hits
  auto nhits = [] (const std::vector<dd4pod::PhotoMultiplierHitData>& evt) {return (int) evt.size(); };

  // Define variables
  auto d1 = d0.Define("nhits", nhits, {"MRICHHits"});

  // Define Histograms
  auto hNhits =
      d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events",
                  100, 0.0, 2000.0},
                 "nhits");
  auto hXYhits =
      d1.Histo2D({"hXYhits", "Hit positions for events; Horizontal position [mm]; Vertical position [mm]",
                  1000, -1000.0, +1000.0, 1000, -1000.0, +1000.0},
                 "MRICHHits.position.x", "MRICHHits.position.y");

  // Event Counts
  auto nevents_thrown = d1.Count();
  std::cout << "Number of Thrown Events: " << (*nevents_thrown) << "\n";

  // Draw Histograms
  {
    TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
    auto h = hXYhits->DrawCopy();
    c1->SaveAs(fmt::format("results/pid/backward/mrich/mrich_{}_hits_xy.png",input_pname).c_str());
    c1->SaveAs(fmt::format("results/pid/backward/mrich/mrich_{}_hits_xy.pdf",input_pname).c_str());
  }
  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c2->SaveAs(fmt::format("results/pid/backward/mrich/mrich_{}_nhits.png",input_pname).c_str());
    c2->SaveAs(fmt::format("results/pid/backward/mrich/mrich_{}_nhits.pdf",input_pname).c_str());
  }
}
