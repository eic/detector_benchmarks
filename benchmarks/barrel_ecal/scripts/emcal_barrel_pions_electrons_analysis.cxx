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
#include "TLegend.h"
#include "TString.h"

R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"
#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "emcal_barrel_common_functions.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void emcal_barrel_pions_electrons_analysis(const char* input_fname1 = "sim_output/sim_emcal_barrel_piminus.root", const char* input_fname2 = "sim_output/sim_emcal_barrel_uniform_electrons.root")
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
  ROOT::RDataFrame d0("events", {input_fname1, input_fname2});

  // Environment Variables
  std::string detector_path = "";
  std::string detector_name = "topside";
  if(std::getenv("DETECTOR_PATH")) {
    detector_path = std::getenv("DETECTOR_PATH");
  }
  if(std::getenv("JUGGLER_DETECTOR")) {
    detector_name = std::getenv("JUGGLER_DETECTOR");
  }

  // Using the detector layers 
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(fmt::format("{}/{}.xml", detector_path,detector_name));
  //dd4hep::rec::CellIDPositionConverter cellid_converter(detector);

  auto decoder = detector.readout("EcalBarrelHits").idSpec().decoder();
  fmt::print("{}\n", decoder->fieldDescription());
  auto layer_index = decoder->index("layer");
  fmt::print(" layer index is {}.\n", layer_index);

  // Rejection Value [GeV] based upon Energy deposit in the first 4 layers.
  // Currently constructed from looking at tthe pion and electron Energy deposit plots
  // should eventually grab a value from a json file
  // Suggested cut at 0.005
  double rejectCut = 0.0;

  // Thrown Energy [GeV]
  auto Ethr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz + input[2].mass*input[2].mass);
  };

  // Thrown Momentum [GeV]
  auto Pthr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz);
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

  // Energy deposititon [GeV] in the first 4 layers
  auto Esim_front = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      //fmt::print("cell id {}, layer {}\n",i.cellID, decoder->get(i.cellID, layer_index));
      if( decoder->get(i.cellID, layer_index) < 4 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Sampling fraction = Esampling / Ethrown
  auto fsam = [](const double& sampled, const double& thrown) {
    return sampled / thrown;
  };

  // E_front / p
  auto fEp = [](const double& E_front, const double& mom) {
    return E_front / mom;
  };

  // Returns the pdgID of the particle
  auto getpid = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return input[2].pdgID;
  };

  // Returns number of particle daughters
  auto getdau = [](std::vector<dd4pod::Geant4ParticleData> const& input){
    return input[2].daughters_begin;
  };

  // Filter function to get electrons
  auto is_electron = [](std::vector<dd4pod::Geant4ParticleData> const& input){
    if (input[2].pdgID == 11){return true;}
    else {return false;}
  };

  // Filter function to get just negative pions
  auto is_piMinus = [](std::vector<dd4pod::Geant4ParticleData> const& input){
    if (input[2].pdgID == -211){return true;}
    else {return false;}
  };

  // Define variables
  auto d1 = d0.Define("Ethr",       Ethr,       {"mcparticles"})
              .Define("Pthr",       Pthr,       {"mcparticles"})
              .Define("nhits",      nhits,      {"EcalBarrelHits"})
              .Define("EsimImg",    Esim,       {"EcalBarrelHits"})
              .Define("EsimScFi",   Esim,       {"EcalBarrelScFiHits"})
              .Define("Esim",       "EsimImg+EsimScFi")
              .Define("fsam",       fsam,       {"Esim","Ethr"})
              .Define("pid",        getpid,     {"mcparticles"})
              .Define("dau",        getdau,     {"mcparticles"})
              .Define("Esim_front", Esim_front, {"EcalBarrelHits"})
              .Define("EOverP",     fEp,        {"Esim_front", "Pthr"})
              ;

  // Particle Filters
  auto d_ele = d1.Filter(is_electron, {"mcparticles"});
  auto d_pim = d1.Filter(is_piMinus,  {"mcparticles"});
  
  // Energy Deposit Filters
  auto d_ele_rej = d_ele.Filter("Esim_front > " + to_string(rejectCut));
  auto d_pim_rej = d_pim.Filter("Esim_front > " + to_string(rejectCut));

  // Define Histograms
  auto hEthr       = d1.Histo1D({"hEthr",        "Thrown Energy; Thrown Energy [GeV]; Events",                      100,  0.0,    7.5}, "Ethr");
  auto hNhits      = d1.Histo1D({"hNhits",       "Number of hits per events; Number of hits; Events",               100,  0.0, 2000.0}, "nhits");
  auto hEsim       = d1.Histo1D({"hEsim",        "Energy Deposit; Energy Deposit [GeV]; Events",                     10,  0.0,   0.05}, "Esim");
  auto hfsam       = d1.Histo1D({"hfsam",        "Sampling Fraction; Sampling Fraction; Events",                    100,  0.0,    0.1}, "fsam");
  auto hEsim_front = d1.Histo1D({"hEsim_front",  "Energy Deposit Front; Energy Deposit [GeV]; Events",               10,  0.0,   0.05}, "Esim_front");
  addDetectorName(detector_name, hEthr.GetPtr());
  addDetectorName(detector_name, hNhits.GetPtr());
  addDetectorName(detector_name, hfsam.GetPtr());
  addDetectorName(detector_name, hEsim_front.GetPtr());

  auto hEsim_ele        = d_ele.Histo1D({"hEsim_ele",        "Energy Deposit Electron; Energy Deposit [GeV]; Events",            10,  0.0,    0.05}, "Esim");
  auto hEsim_ele_front  = d_ele.Histo1D({"hEsim_ele_front",  "Energy Deposit Front Electron; Energy Deposit [GeV]; Events",      10,  0.0,    0.05}, "Esim_front");
  auto hEsim_pim        = d_pim.Histo1D({"hEsim_pim",        "Energy Deposit Electron; Energy Deposit [GeV]; Events",            10,  0.0,    0.05}, "Esim");
  auto hEsim_pim_front  = d_pim.Histo1D({"hEsim_pim_front",  "Energy Deposit Front Pion-; Energy Deposit [GeV]; Events",         10,  0.0,    0.05}, "Esim_front");
  addDetectorName(detector_name, hEsim_ele.GetPtr());
  addDetectorName(detector_name, hEsim_pim_front.GetPtr());
  addDetectorName(detector_name, hEsim_pim.GetPtr());
  addDetectorName(detector_name, hEsim_pim_front.GetPtr());

  auto hEsim_ele_front_rej  = d_ele_rej.Histo1D({"hEsim_ele_front_rej",  "Energy Deposit Front Electron; Energy Deposit [GeV]; Events",      10,  0.0,    0.05}, "Esim_front");
  auto hEsim_pim_front_rej  = d_pim_rej.Histo1D({"hEsim_pim_front_rej",  "Energy Deposit Front Pion-; Energy Deposit [GeV]; Events",         10,  0.0,    0.05}, "Esim_front");
  addDetectorName(detector_name, hEsim_ele_front_rej.GetPtr());
  addDetectorName(detector_name, hEsim_pim_front_rej.GetPtr()); 

  auto hEpvp_ele = d_ele.Histo2D({"hEpvp_ele", "Energy Deposit 1st 4 Layers/P vs P; E/P; P [GeV]", 20, 0.0, 0.01, 20, 1.0, 7.0}, "EOverP", "Pthr");
  auto hEpvp_pim = d_pim.Histo2D({"hEpvp_pim", "Energy Deposit 1st 4 Layers/P vs P; E/P; P [GeV]", 20, 0.0, 0.01, 20, 1.0, 7.0}, "EOverP", "Pthr");
  addDetectorName(detector_name, hEpvp_ele.GetPtr());
  addDetectorName(detector_name, hEpvp_pim.GetPtr());


  TH1D* hElePurity_initial = (TH1D *)hEsim_ele -> Clone();
  hElePurity_initial -> Divide(hEsim.GetPtr());
  hElePurity_initial -> SetTitle("Electron/Pion Rejection");
  addDetectorName(detector_name, hElePurity_initial);

  TH1D* hElePurity_ele = (TH1D *)hEsim_ele_front -> Clone();
  hElePurity_ele -> Divide(hEsim_front.GetPtr());
  hElePurity_ele -> SetTitle("Electron/Pion Rejection : Electron");
  addDetectorName(detector_name, hElePurity_ele);

  TH1D* hElePurity_pim = (TH1D *)hEsim_pim_front -> Clone();
  hElePurity_pim -> Divide(hEsim_front.GetPtr());
  hElePurity_pim -> SetTitle("Electron/Pion Rejection : Pion Minus");
  addDetectorName(detector_name, hElePurity_pim);

  // Rejection
  TH1D* hElePurity_rej = (TH1D *)hEsim_ele_front_rej -> Clone();
  hElePurity_rej -> Add(hEsim_pim_front_rej.GetPtr());
  hElePurity_rej -> SetTitle("Electron/Pion Rejection");
  addDetectorName(detector_name, hElePurity_rej);

  TH1D* hElePurity_ele_rej = (TH1D *)hEsim_ele_front_rej -> Clone();
  hElePurity_ele_rej -> Divide(hElePurity_rej);
  hElePurity_ele_rej -> SetTitle("Electron/Pion Rejection : Electron");
  addDetectorName(detector_name, hElePurity_ele_rej);

  TH1D* hElePurity_pim_rej = (TH1D *)hEsim_pim_front_rej -> Clone();
  hElePurity_pim_rej -> Divide(hElePurity_rej);
  hElePurity_pim_rej -> SetTitle("Electron/Pion Rejection : Pi Minus");
  addDetectorName(detector_name, hElePurity_pim_rej);

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
  //c5->SetLogy(1);
  hElePurity_initial->GetYaxis()->SetTitleOffset(1.4);
  hElePurity_initial->SetLineWidth(2);
  hElePurity_initial->SetLineColor(kBlue);
  hElePurity_initial->DrawClone();
  c5->SaveAs("results/emcal_barrel_pions_electrons_rejection_initial.png");
  c5->SaveAs("results/emcal_barrel_pions_electrons_rejection_initial.pdf");

  TCanvas *c6 = new TCanvas("c6", "c6", 700, 500);
  //c6->SetLogy(1);
  hElePurity_ele->GetYaxis()->SetTitleOffset(1.4);
  hElePurity_ele->SetLineWidth(2);
  hElePurity_ele->SetLineColor(kBlue);
  hElePurity_ele->DrawClone();
  hElePurity_ele_rej->SetLineWidth(2);
  hElePurity_ele_rej->SetLineColor(kRed);
  hElePurity_ele_rej->SetLineStyle(10);
  hElePurity_ele_rej->DrawClone("Same");
  c6->SaveAs("results/emcal_barrel_pions_electrons_rejection_ele.png");
  c6->SaveAs("results/emcal_barrel_pions_electrons_rejection_ele.pdf");

  auto leg = new TLegend(0.7, 0.8, 0.8, 0.9);
  leg->AddEntry(hElePurity_initial, "Initial", "l");
  leg->AddEntry(hElePurity_ele, "Final", "l");

  TCanvas *c7 = new TCanvas("c7", "c7", 700, 500);
  //c6->SetLogy(1);
  hElePurity_pim->GetYaxis()->SetTitleOffset(1.4);
  hElePurity_pim->SetLineWidth(2);
  hElePurity_pim->SetLineColor(kBlue);
  hElePurity_pim->DrawClone();
  hElePurity_pim_rej->SetLineWidth(2);
  hElePurity_pim_rej->SetLineColor(kRed);
  hElePurity_pim_rej->SetLineStyle(10);
  hElePurity_pim_rej->DrawClone("Same");
  c7->SaveAs("results/emcal_barrel_pions_electrons_rejection_pim.png");
  c7->SaveAs("results/emcal_barrel_pions_electrons_rejection_pim.pdf");

  TCanvas *c8 = new TCanvas("c8", "c8", 700, 500);
  //c6->SetLogy(1);
  hEpvp_ele->GetYaxis()->SetTitleOffset(1.4);
  hEpvp_ele->DrawClone("COLZ");
  c8->SaveAs("results/emcal_barrel_pions_electrons_Epvp_ele.png");
  c8->SaveAs("results/emcal_barrel_pions_electrons_Epvp_ele.pdf");

  TCanvas *c9 = new TCanvas("c9", "c9", 700, 500);
  //c6->SetLogy(1);
  hEpvp_pim->GetYaxis()->SetTitleOffset(1.4);
  hEpvp_pim->DrawClone("COLZ");
  c9->SaveAs("results/emcal_barrel_pions_electrons_Epvp_pim.png");
  c9->SaveAs("results/emcal_barrel_pions_electrons_Epvp_pim.pdf");

}
