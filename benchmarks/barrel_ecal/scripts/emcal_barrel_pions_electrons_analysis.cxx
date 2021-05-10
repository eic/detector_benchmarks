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

R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"
#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"


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

  //Using the detector layers 
  std::string detector_path = "";
  std::string detector_name = "topside";
  if(std::getenv("DETECTOR_PATH")) {
    detector_path = std::getenv("DETECTOR_PATH");
  }
  if(std::getenv("JUGGLER_DETECTOR")) {
    detector_name = std::getenv("JUGGLER_DETECTOR");
  }

  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(fmt::format("{}/{}.xml", detector_path,detector_name));
  //dd4hep::rec::CellIDPositionConverter cellid_converter(detector);

  auto decoder = detector.readout("EcalBarrelHits").idSpec().decoder();
  fmt::print("{}\n", decoder->fieldDescription());
  auto layer_index = decoder->index("layer");
  fmt::print(" layer index is {}.\n", layer_index);



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

  // Energy deposition electron
  auto Esim_ele = [](const std::vector<dd4pod::CalorimeterHitData>& evt, std::vector<dd4pod::Geant4ParticleData> const& input) {
    std::vector<double> result;
    auto total_edep = 0.0;
    int count = 0;
    if (input[2].pdgID == 11)// Electron
    { 
      for (const auto& i: evt)
        total_edep += i.energyDeposit;
        count++;
    }
    result.push_back(total_edep);
  return result;
  };

  // Energy deposition [GeV] pion
  auto Esim_pi = [](const std::vector<dd4pod::CalorimeterHitData>& evt, std::vector<dd4pod::Geant4ParticleData> const& input) {
    std::vector<double> result;
    auto total_edep = 0.0;
    int count = 0;
    if (input[2].pdgID == -211)// Negative pion
    { 
      for (const auto& i: evt)
        total_edep += i.energyDeposit;
        count++;
    }
    result.push_back(total_edep);
  return result;
  };

  auto Esim_front = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      fmt::print("cell id {}, layer {}\n",i.cellID, decoder->get(i.cellID, layer_index));
      if( decoder->get(i.cellID, layer_index) < 4 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
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
    //std::vector<int> result = {input[2].pdgID, input[3].pdgID};
    return input[2].pdgID;
  };

  // Returns number of particle daughters
  auto getdau = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return input[2].daughters_begin;
  };

  auto is_electron = [](std::vector<dd4pod::Geant4ParticleData> const& input){
    if (input[2].pdgID == 11){return true;}
    else {return false;}
  };

  // Define variables
  auto d1 = d0.Define("Ethr",   Ethr,       {"mcparticles"})
              .Define("nhits",  nhits,      {"EcalBarrelHits"})
              .Define("Esim",   Esim,       {"EcalBarrelHits"})
              .Define("fsam",   fsam,       {"Esim","Ethr"})
              .Define("pid",    getpid,     {"mcparticles"})
              .Define("dau",    getdau,     {"mcparticles"})
              .Define("Esim_ele", Esim_ele, {"EcalBarrelHits", "mcparticles"})
              .Define("Esim_pi", Esim_pi,   {"EcalBarrelHits", "mcparticles"})
              .Define("Esim_ele_red", Esim_ele, {"EcalBarrelHits", "mcparticles"})
              .Define("Esim_pi_red", Esim_pi,   {"EcalBarrelHits", "mcparticles"})
              .Define("Esim_front", Esim_front, {"EcalBarrelHits"})
              ;

  auto d2 = d1.Filter(is_electron, {"mcparticles"});

  // Define Histograms
  auto hEthr  = d1.Histo1D({"hEthr",  "Thrown Energy; Thrown Energy [GeV]; Events",        100,  0.0,    7.5}, "Ethr");
  auto hNhits = d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events", 100,  0.0, 2000.0}, "nhits");
  auto hEsim  = d1.Histo1D({"hEsim",  "Energy Deposit; Energy Deposit [GeV]; Events",      10,  0.0,    0.5}, "Esim");
  auto hfsam  = d1.Histo1D({"hfsam",  "Sampling Fraction; Sampling Fraction; Events",      100,  0.0,    0.1}, "fsam");
  auto hpid   = d1.Histo1D({"hpid",   "PID; PID; Count",                                   100,  -220,   220}, "pid");
  auto hdau   = d1.Histo1D({"hdau",   "Number of Daughters; Number of Daughters; Count",   10,   0,      10},  "dau");

  auto hEsim_ele        = d2.Histo1D({"hEsim_ele",        "Energy Deposit Electron; Energy Deposit [GeV]; Events",      10,  0.0,    0.5}, "Esim");
  auto hEsim_ele_front  = d2.Histo1D({"hEsim_ele_front",  "Energy Deposit Front Electron; Energy Deposit [GeV]; Events",      10,  0.0,    0.5}, "Esim_front");
  auto hEsim_front      = d1.Histo1D({"hEsim_front",      "Energy Deposit Front Electron; Energy Deposit [GeV]; Events",      10,  0.0,    0.5}, "Esim_front");

  TH1D* hElePurity_initial = (TH1D *)hEsim_ele -> Clone();
  hElePurity_initial -> Divide(hEsim.GetPtr());

  TH1D* hElePurity_final = (TH1D *)hEsim_ele_front -> Clone();
  hElePurity_initial -> Divide(hEsim_front.GetPtr());


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
  c1->SaveAs("results/emcal_barrel_pions_electrons_Ethr.png");
  c1->SaveAs("results/emcal_barrel_pions_electrons_Ethr.pdf");

  TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);
  c2->SetLogy(1);
  hNhits->GetYaxis()->SetTitleOffset(1.4);
  hNhits->SetLineWidth(2);
  hNhits->SetLineColor(kBlue);
  hNhits->DrawClone();
  c2->SaveAs("results/emcal_barrel_pions_electrons_nhits.png");
  c2->SaveAs("results/emcal_barrel_pions_electrons_nhits.pdf");

  TCanvas *c3 = new TCanvas("c3", "c3", 700, 500);
  c3->SetLogy(1);
  hEsim->GetYaxis()->SetTitleOffset(1.4);
  hEsim->SetLineWidth(2);
  hEsim->SetLineColor(kBlue);
  hEsim->DrawClone();
  c3->SaveAs("results/emcal_barrel_pions_electrons_Esim.png"); 
  c3->SaveAs("results/emcal_barrel_pions_electrons_Esim.pdf");

  TCanvas *c4 = new TCanvas("c4", "c4", 700, 500);
  c4->SetLogy(1);
  hEsim_ele->GetYaxis()->SetTitleOffset(1.4);
  hEsim_ele->SetLineWidth(2);
  hEsim_ele->SetLineColor(kBlue);
  hEsim_ele->DrawClone();
  c4->SaveAs("results/emcal_barrel_pions_electrons_Esim_ele.png");
  c4->SaveAs("results/emcal_barrel_pions_electrons_Esim_ele.pdf");

  TCanvas *c5 = new TCanvas("c5", "c5", 700, 500);
  c5->SetLogy(1);
  hEsim_pi->GetYaxis()->SetTitleOffset(1.4);
  hEsim_pi->SetLineWidth(2);
  hEsim_pi->SetLineColor(kBlue);
  hEsim_pi->DrawClone();
  c5->SaveAs("results/emcal_barrel_pions_electrons_Esim_pi.png");
  c5->SaveAs("results/emcal_barrel_pions_electrons_Esim_pi.pdf");

  TCanvas *c6 = new TCanvas("c6", "c6", 700, 500);
  c6->SetLogy(1);
  hpid->GetYaxis()->SetTitleOffset(1.4);
  hpid->SetLineWidth(2);
  hpid->SetLineColor(kBlue);
  hpid->DrawClone();
  c6->SaveAs("results/emcal_barrel_pions_electrons_pid.png");
  c6->SaveAs("results/emcal_barrel_pions_electrons_pid.pdf");

  TCanvas *c7 = new TCanvas("c7", "c7", 700, 500);
  //c7->SetLogy(1);
  hElePurity_initial->GetYaxis()->SetTitleOffset(1.4);
  hElePurity_initial->SetLineWidth(2);
  hElePurity_initial->SetLineColor(kBlue);
  hElePurity_initial->DrawClone();
  hElePurity_final->SetLineWidth(2);
  hElePurity_final->SetLineColor(kRed);
  hElePurity_final->DrawClone("SAME");
  c7->SaveAs("results/emcal_barrel_pions_electrons_purity.png");
  c7->SaveAs("results/emcal_barrel_pions_electrons_purity.pdf");


}
