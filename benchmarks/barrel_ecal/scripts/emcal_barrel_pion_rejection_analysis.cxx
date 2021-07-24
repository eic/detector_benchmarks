////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include <algorithm>
#include <string>

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"

#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"

#include <boost/range/combine.hpp>

#include "DD4hep/IDDescriptor.h"
#include "DD4hep/Readout.h"
#include "DD4hep/Segmentations.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TError.h"

R__LOAD_LIBRARY(libfmt.so)
#include "fmt/core.h"
#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "emcal_barrel_common_functions.h"
#include <nlohmann/json.hpp>

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void emcal_barrel_pion_rejection_analysis(
                                          const char* input_fname1 = "sim_output/sim_emcal_barrel_piRej_electron.root",
                                          const char* input_fname2 = "sim_output/sim_emcal_barrel_piRej_piminus.root"
                                          )
{
  // Error Ignore Level Set
  gErrorIgnoreLevel = kFatal;

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
  std::string detector_name = "athena";//athena
  if(std::getenv("DETECTOR_PATH")) {
    detector_path = std::getenv("DETECTOR_PATH");
  }
  if(std::getenv("JUGGLER_DETECTOR")) {
    detector_name = std::getenv("JUGGLER_DETECTOR");
  }

  /*
  // Sampling Fraction grabbed from json file
  // Note that this value is derived from electrons
  json j;
  std::ifstream prev_steps_ifstream("results/emcal_barrel_calibration.json");
  prev_steps_ifstream >> j;

  // Sampling Fraction
  double samp_frac = j["electron"]["sampling_fraction"];
  */

  // Detector Layer Variables
  int layerNum; 
  int dep_min = 1;
  int dep_max = 6;

  // DD4HEP interface 
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(fmt::format("{}/{}.xml", detector_path,detector_name));

  auto decoder         = detector.readout("EcalBarrelHits").idSpec().decoder();
  auto decoderScFi     = detector.readout("EcalBarrelScFiHits").idSpec().decoder();
  auto layer_index     = decoder->index("layer");
  auto layer_indexScFi = decoderScFi->index("layer");

  // Thrown Energy [GeV]
  auto Ethr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz + input[2].mass*input[2].mass);
  };

  // Thrown Momentum [GeV]
  auto Pthr = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz);
  };

  // Thrown Eta 
  auto Eta = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    double E  = TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz + input[2].mass*input[2].mass);
    return 0.5*TMath::Log((E + input[2].psz) / (E - input[2].psz));
  };

  // Thrown pT [GeV]
  auto pT = [](std::vector<dd4pod::Geant4ParticleData> const& input) {
    return TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy);
  };

  // Number of hits
  auto nhits = [] (const std::vector<dd4pod::CalorimeterHitData>& evt) {return (int) evt.size(); };

  // Energy deposition [GeV]
  auto Esim = [](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    double total_edep = 0.0;
    for (const auto& i: evt){
      total_edep += i.energyDeposit;
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 2 layers
  auto Esim_dep2 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 3 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 3 layers
  // Same as Esim_front from previous codes
  auto Esim_dep3 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 4 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 4 layers
  auto Esim_dep4 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 5 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 5 layers
  auto Esim_dep5 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 6 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 6 layers
  auto Esim_dep6 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 7 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 6 layers
  auto Esim_dep6_ScFi = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoderScFi->get(i.cellID, layer_indexScFi) < 7 ){
        total_edep += i.energyDeposit;
      } 
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 7 layers
  auto Esim_dep7 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 8 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };
    // Energy deposititon [GeV] in the first 8 layers
  auto Esim_dep8 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 9 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV] in the first 9 layers
  auto Esim_dep9 = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    auto total_edep = 0.0;
    for (const auto& i: evt) {
      if( decoder->get(i.cellID, layer_index) < 10 ){
        total_edep += i.energyDeposit;
      }
    }
    return total_edep;
  };

  // Energy deposititon [GeV], returns array
  // Note Layer_index = 0 does not exist at the wrting of this code
  auto Esim_dep = [=](const std::vector<dd4pod::CalorimeterHitData>& evt) {
    std::vector<double>res(20);
    for (const auto& i: evt) {
      res[decoder->get(i.cellID, layer_index)] += i.energyDeposit;
    }
    return res;
  };

  // Sum of Energy deposititon [GeV]
  auto Esim_dep_sum = [&dep_min, &dep_max](const std::vector<double>& dep) {
    double res = 0;
    for (int i = dep_min; i < dep_max + 1; i++) {
      if (i >= dep.size()) continue;
      res += dep[i];
    }
    return res;
  };

  // Energy deposititon in a layer [GeV]
  auto Esim_depN = [&layerNum](const std::vector<double>& dep) {
    return (layerNum < dep.size() ? dep[layerNum] : 0.0);
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
    return (input[2].pdgID == 11 ? true : false);
  };

  // Filter function to get just negative pions
  auto is_piMinus = [](std::vector<dd4pod::Geant4ParticleData> const& input){
    return (input[2].pdgID == -211 ? true : false);
  };

  // Filter function Edeposit
  auto EDep_bool = [](std::vector<double> const& input){
    return (input[0] > input[1]);
  };

  auto Diff = [](const double &esim, const double &edep){
    return esim-edep;
  };


  // Define variables
  auto d1 = d0.Define("Ethr",            Ethr,                  {"mcparticles"})
              .Define("Pthr",            Pthr,                  {"mcparticles"})
              .Define("nhits",           nhits,                 {"EcalBarrelHits"})
              .Define("Esim",            Esim,                  {"EcalBarrelHits"})
              .Define("EsimScFi",        Esim,                  {"EcalBarrelScFiHits"})
              .Define("EsimOverP",       fEp,                   {"Esim", "Pthr"})
              .Define("EsimScFiOverP",   fEp,                   {"EsimScFi", "Pthr"})
              .Define("EsimTot",                                "EsimScFi+Esim")
              .Define("EsimTotOverP",    fEp,                   {"EsimTot", "Pthr"})
              .Define("fsam",            fsam,                  {"Esim","Ethr"})
              .Define("pid",             getpid,                {"mcparticles"})
              .Define("EDep",            Esim_dep,              {"EcalBarrelHits"})
              .Define("EDepSum",         Esim_dep_sum,          {"EDep"})
              .Define("EDepN",           Esim_depN,             {"EDep"})
              .Define("EDep2",           Esim_dep2,             {"EcalBarrelHits"})
              .Define("EDep3",           Esim_dep3,             {"EcalBarrelHits"})
              .Define("EDep4",           Esim_dep4,             {"EcalBarrelHits"})
              .Define("EDep5",           Esim_dep5,             {"EcalBarrelHits"})
              .Define("EDep6",           Esim_dep6,             {"EcalBarrelHits"})
              .Define("EDep6OverP",      fEp,                   {"EDep6", "Pthr"})
              .Define("EOverP",          fEp,                   {"EDep3", "Pthr"})
              .Define("Eta",             Eta,                   {"mcparticles"})
              .Define("pT",              pT,                    {"mcparticles"})
              .Define("EDepOverP",       fEp,                   {"EDepN", "Pthr"})
              .Define("EDepOverPT",      fEp,                   {"EDepN", "pT"})
              .Define("EDepSumOverP",    fEp,                   {"EDepSum", "Pthr"})
              .Define("EDepSumOverPT",   fEp,                   {"EDepSum", "pT"})
              .Define("EDepFrac",        fEp,                   {"EDepSum", "Esim"})
              ;
  
  // Particle Filters
  dep_min = 1;
  dep_max = 6;
  auto d_ele = d1.Filter(is_electron, {"mcparticles"});
  auto d_pim = d1.Filter(is_piMinus,  {"mcparticles"});

  // Cut Filter
  std::string currentCut = "(EDep6OverP>2.5e-3)&&(EDep6>5e-3)";// Good athena cut, that is changed by cutEE later

  // Generic 1D Histogram Plots Comparing Electons and Pions w/o cuts
  // Edep first 6 layers(EDep6), EDep/p, pT, eta
  std::vector<std::string> var              = {"Esim [GeV];", "EsimTot [GeV];", "EDep6 [GeV];", "EDep6/p;",     "pT [GeV];", "#eta;",  "EsimScFi [GeV]",  "EsimScFi/p"};
  std::vector<std::string> var_save         = {"Esim",        "EsimTot",        "EDep6",        "EDep6OverP",   "pT",        "eta",    "EsimScFi",        "EsimScFiOverP"};  
  std::vector<std::string> col              = {"Esim",        "EsimTot",        "EDep6",        "EDep6OverP",   "pT",        "Eta",    "EsimScFi",        "EsimScFiOverP"};
  std::vector<std::vector<double>> h1Ranges = {{0,0.2},       {0, 0.2},         {0,0.25},       {0, 0.02},      {0, 18},     {-1, 1},  {0,0.2},            {0,0.2}};
  for (int i = 0; i < var.size(); i++){
    std::string title = "#pi^{-}, e^{-};" + var[i] + " Events";
    auto he = d_ele.Histo1D({"he", title.c_str(), 100, h1Ranges[i][0], h1Ranges[i][1]}, col[i]);
    auto hp = d_pim.Histo1D({"hp", title.c_str(), 100, h1Ranges[i][0], h1Ranges[i][1]}, col[i]);

    hp->GetYaxis()->SetTitleOffset(1.4);
    he->SetLineWidth(2);
    he->SetLineColor(kRed);
    hp->SetLineWidth(2);
    hp->SetLineColor(kBlue);
    auto c = new TCanvas("c", "c", 700, 500);
    auto leng = new TLegend(0.7, 0.7, 0.9, 0.9);
    if (var[i] != "EsimScFi/p"){ 
      hp->DrawClone();
      he->DrawClone("same");
    }
    else {
      he->DrawClone();
      hp->DrawClone("same");
    }
    c->Update();

    leng->AddEntry(he.GetPtr(),"e^{-}","l");
    leng->AddEntry(hp.GetPtr(),"#pi^{-}","l");
    leng->Draw();
    c->SaveAs(("results/emcal_barrel_pion_rej_uncut_comb_" + var_save[i] + ".png").c_str());
  }

  // Cut Generation
  // The cut breaks the Energy range in to three energy bins (EBins) and the two barrel eta bins (-1 to 0, and 0 to 1)
  // Then fits a gaussian within the range
  // The cut is then based upon Mean -2*StdDev < Mean < 3*StdDev
  std::string cutEEta;
  std::vector<std::vector<double>> EBins = {{0,2}, {2, 4}, {4, 6}, {6, 9}, {9, 12}, {12, 19}};
  for (int i = 0; i < EBins.size(); i++){
    std::string minCut = "Pthr>="+std::to_string(EBins[i][0]);
    std::string maxCut = "Pthr<"+std::to_string(EBins[i][1]);
    cutEEta += "(" + minCut + "&&" + maxCut + "&&";
    
    for (int j = -1; j < 1; j++){
      std::string title = "#pi^{-}, e^{-}";
      title += fmt::format(" : {} < E < {}", EBins[i][0], EBins[i][1]);
      title += fmt::format(" & {} < #eta < {};EDep6/p; Events", j, j+1);
      std::string etaCutMin = fmt::format("Eta>={}", j);
      std::string etaCutMax = fmt::format("Eta<{}",j+1);
      cutEEta += "(" + etaCutMin + "&&" + etaCutMax + "&&";
      auto he = d_ele.Filter(minCut).Filter(maxCut).Filter(etaCutMin).Filter(etaCutMax).Histo1D({"he", title.c_str(), 50, 0, 0.02}, "EDep6OverP");
      auto hp = d_pim.Filter(minCut).Filter(maxCut).Filter(etaCutMin).Filter(etaCutMax).Histo1D({"hp", title.c_str(), 50, 0, 0.02}, "EDep6OverP");
      auto hecopy = he->DrawCopy();
      hecopy->Fit("gaus", "", "", 0, hecopy->GetMaximum());
      double* res = hecopy->GetFunction("gaus")->GetParameters();
      cutEEta += fmt::format("EDep6OverP>={}", res[1] - 2.0*res[2]);
      cutEEta += fmt::format("&&EDep6OverP<{})||",res[1] + 3.0*res[2]);

      hp->GetYaxis()->SetTitleOffset(1.4);
      he->SetLineWidth(2);
      he->SetLineColor(kRed);
      hp->SetLineWidth(2);
      hp->SetLineColor(kBlue);
      auto c = new TCanvas("c", "c", 700, 500);
      auto leng = new TLegend(0.7, 0.7, 0.9, 0.9);
      hp->DrawClone();
      he->DrawClone("same");
      c->Update();

      leng->AddEntry(he.GetPtr(),"e^{-}","l");
      leng->AddEntry(hp.GetPtr(),"#pi^{-}","l");
      leng->Draw();
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_uncut_comb_E{}Eta{}.png", i, j+1)).c_str());
    }
    cutEEta.pop_back();
    cutEEta.pop_back();
    cutEEta += ")||";
  }
  cutEEta.pop_back();
  cutEEta.pop_back();
  currentCut = cutEEta;

  // Filtered dataframes
  auto d_pim_cut = d_pim.Filter(currentCut.c_str());
  auto d_ele_cut = d_ele.Filter(currentCut.c_str());

  // Gathering benchmarks and plotting distrubutions
  std::vector<double> E                 = {5, 10, 18};
  std::vector<double> ledges5           = {2.8, 0.4, 0.3, 0.5};
  std::vector<double> ledges10          = {1.4, 0.5, 0.6, 1.0};
  std::vector<double> ledges18          = {0.9, 0.9, 1.0, 1.8};
  std::vector<double> maxRate5          = {0.1, 100, 500, 1000};
  std::vector<double> maxRate10         = {10, 400, 800, 1000};
  std::vector<double> maxRate18         = {200, 800, 1000, 100};
  std::vector<vector<double>> maxRate   = {maxRate5, maxRate10, maxRate18};
  std::vector<vector<double>> lowEdges  = {ledges5, ledges10, ledges18};
  double suppression                    = 1e-4;
  std::vector<vector<double>> rejRatios = lowEdges;
  std::vector<vector<double>> effEle    = lowEdges;
  std::vector<vector<double>> effPim    = lowEdges;
  std::vector<std::string> etaBin       = {"Eta >= -3.5 && Eta < -2.0", "Eta >= -2.0 && Eta < -1.0", "Eta >= -1.0 && Eta < 0", "Eta >= 0 && Eta < 1.0"};
  std::vector<std::string> etaTitle     = {"-3.5 < #eta < -2.0", "-2.0 < #eta < -1.0", "-1.0 < #eta < 0", "0 < #eta < 1.0"};

  // Pion Rejection Plot that mimics that of the one in the image
  std::vector<double> pBins  = {0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,10,12,14,16,18};
  EBins = {{0,2}, {2, 4}, {4, 6}, {6, 9}, {9, 12}, {12, 19}};
  auto tg = new TGraphErrors();

  for (int i = 0; i < EBins.size(); i++){
    std::string filter = currentCut;
    filter += "&&(Pthr>=" + std::to_string(EBins[i][0]) + "&&Pthr<" + std::to_string(EBins[i][1]) + ")";
    double numer = (double)*d_ele_cut.Filter(filter.c_str()).Count();
    double denom = (double)*d_pim_cut.Filter(filter.c_str()).Count();
    double error = std::sqrt(std::pow(numer / denom, 2.0)*(1.0/numer + 1.0/denom));
    double ratio = numer / denom;
    if (denom == 0){ratio = 1; error = 1;}
    tg->SetPoint(i, 0.5*(EBins[i][0] + EBins[i][1]), ratio);
    tg->SetPointError(i, 0.5*(EBins[i][1] - EBins[i][0]), error);
  }
  double e_eff = (double)*d_ele_cut.Count() / (double)*d_ele.Count();
  tg->SetTitle(("#pi Rejection with #varepsilon_{e} = "+ std::to_string(e_eff)).c_str());
  tg->GetXaxis()->SetTitle("p [GeV]");
  tg->GetYaxis()->SetTitle("R_{e/#pi}");
  tg->SetMarkerColor(kBlue);
  tg->SetMarkerStyle(20);

  auto cp = new TCanvas("cp", "cp");
  cp->SetLogy();
  cp->SetLogx();
  tg->DrawClone("ap");
  cp->SaveAs("results/emcal_barrel_pion_rej_RatioRej.png");

  // Barrel eta cuts
  // The eta range for the barrel is -1 < eta < 1
  // Threfore the first bins are empty and will not be iterated over
  dep_min = 1;
  dep_max = 6;
  for (int i = 0; i < 3; i++){   // E loop
    for (int j = 2; j < 4; j++){ // Eta Looop
      
      // Apply eta cuts/binning and Momentum Cut
      std::string pCut = "Pthr>=" + std::to_string(lowEdges[i][j]) + "&&Pthr<" + std::to_string(E[i]);
      auto e_eta = d_ele.Filter(etaBin[j]).Filter(pCut);
      auto p_eta = d_pim.Filter(etaBin[j]).Filter(pCut);
      
      // Print out the momentum distributions for the electron and pi-
      std::string title = "e^{-} (E = " + std::to_string((int)E[i]) + " GeV) : " + etaTitle[j] + "; p [GeV]; Events";
      auto he = e_eta.Histo1D({"he", title.c_str(), 100, lowEdges[i][j], E[i]}, "Pthr");
      he->SetLineColor(kBlue);
      auto he_cut = e_eta.Filter(currentCut).Histo1D({"he_cut", title.c_str(), 100, lowEdges[i][j], E[i]}, "Pthr");
      he_cut->GetYaxis()->SetTitleOffset(1.4);
      he_cut->SetLineWidth(2);
      he_cut->SetLineColor(kRed);

      title = "#pi^{-} (E = " + std::to_string((int)E[i]) + " GeV) : " + etaTitle[j] + "; p [GeV]; Events";
      auto hp = p_eta.Histo1D({"hp", title.c_str(), 100, lowEdges[i][j], E[i]}, "Pthr");
      hp->SetLineColor(kBlue);
      auto hp_cut = p_eta.Filter(currentCut).Histo1D({"hp", title.c_str(), 100, lowEdges[i][j], E[i]}, "Pthr");
      hp_cut->GetYaxis()->SetTitleOffset(1.4);
      hp_cut->SetLineWidth(2);
      hp_cut->SetLineColor(kRed);

      auto c = new TCanvas("c", "c", 700, 500);
      c->SetLogy(1);
      he->DrawClone();
      he_cut->DrawClone("same");
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_mom_ele_E{}_eta{}.png", (int)E[i], j)).c_str());
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_mom_ele_E{}_eta{}.pdf", (int)E[i], j)).c_str());
      c->Clear();

      hp->DrawClone();
      hp_cut->DrawClone("same");
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_mom_pim_E{}_eta{}.png", (int)E[i], j)).c_str());
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_mom_pim_E{}_eta{}.pdf", (int)E[i], j)).c_str());
      c->Clear();

      // Gather ratio of pi/e and efficiencies for each Energy and eta bin
      // Then plot the distributions
      rejRatios[i][j] = (double)hp_cut->Integral() / (double)he_cut->Integral();
      effPim[i][j]    = (double)hp_cut->Integral() / (double)hp->Integral();
      effEle[i][j]    = (double)he_cut->Integral() / (double)he->Integral();

      hp_cut->Divide(he.GetPtr());
      title = "#pi^{-}/e^{-} (E = " + std::to_string((int)E[i]) + " GeV) : " + etaTitle[j];
      hp_cut->SetTitle(title.c_str());
      hp_cut->SetLineColor(kBlack);
      hp_cut->DrawClone();
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_ratio_pim_E{}_eta{}.png", (int)E[i], j)).c_str());
      c->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_ratio_pim_E{}_eta{}.pdf", (int)E[i], j)).c_str());
      c->Clear();
      
      // Print out the 1D distributions for the electron and pi- within the current Eta bin
      std::vector<std::string> endStr           = {";pT [GeV]; Events", ";EDep6/p; Events"};
      std::vector<std::string> var_save_loc     = {"pT",                "EDep6OverP"};  
      std::vector<std::string> col_loc          = {"pT",                "EDep6OverP"};
      std::vector<std::vector<double>> h1Ranges = {{0, E[i]},           {0, 0.02}};
      for (int k = 0; k < 2; k++){
        auto cl = new TCanvas("cl", "cl", 700, 500);
        title = "e^{-} (E = " + std::to_string((int)E[i]) + " GeV) : " + etaTitle[j] + endStr[k];
        auto he1 = e_eta.Histo1D({"he", title.c_str(), 100, h1Ranges[k][0], h1Ranges[k][1]}, col_loc[k]);
        he1->SetLineColor(kBlue);
        auto he1_cut = e_eta.Filter(currentCut).Histo1D({"he_cut", title.c_str(), 100, h1Ranges[k][0], h1Ranges[k][1]}, col_loc[k]);
        he1_cut->GetYaxis()->SetTitleOffset(1.4);
        he1_cut->SetLineWidth(2);
        he1_cut->SetLineColor(kRed);

        title = "#pi^{-} (E = " + std::to_string((int)E[i]) + " GeV) : " + etaTitle[j] + endStr[k];
        auto hp1 = p_eta.Histo1D({"hp", title.c_str(), 100, h1Ranges[k][0], h1Ranges[k][1]}, col_loc[k]);
        hp1->SetLineColor(kBlue);
        auto hp1_cut = p_eta.Filter(currentCut).Histo1D({"hp_cut", title.c_str(), 100, h1Ranges[k][0], h1Ranges[k][1]}, col_loc[k]);
        hp1_cut->GetYaxis()->SetTitleOffset(1.4);
        hp1_cut->SetLineWidth(2);
        hp1_cut->SetLineColor(kRed);

        cl->SetLogy(1);
        he1->DrawClone();
        he1_cut->DrawClone("same");
        cl->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_{}_ele_E{}_eta{}.png", col_loc[k], (int)E[i], j)).c_str());
        cl->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_{}_ele_E{}_eta{}.pdf", col_loc[k], (int)E[i], j)).c_str());
        cl->Clear();

        hp1->DrawClone();
        hp1_cut->DrawClone("same");
        cl->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_{}_pim_E{}_eta{}.png", col_loc[k], (int)E[i], j)).c_str());
        cl->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_{}_pim_E{}_eta{}.pdf", col_loc[k], (int)E[i], j)).c_str());
        cl->Clear();

        // Combined plots
        title = "#pi^{-}, e^{-} (E = " + std::to_string((int)E[i]) + " GeV) : " + etaTitle[j] + endStr[k];
        hp1_cut->SetLineColor(kBlue);
        he1_cut->SetLineColor(kRed);
        he1_cut->SetTitle(title.c_str());

        auto leng = new TLegend(0.7, 0.7, 0.9, 0.9);
        he1_cut->DrawClone();
        hp1_cut->DrawClone("same");
        cl->Update();

        leng->AddEntry(he1_cut.GetPtr(),"e^{-}","l");
        leng->AddEntry(hp1_cut.GetPtr(),"#pi^{-}","l");
        leng->Draw();
        cl->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_{}_comb_E{}_eta{}.png", col_loc[k], (int)E[i], j)).c_str());
        cl->SaveAs((fmt::format("results/emcal_barrel_pion_rej_cut_{}_comb_E{}_eta{}.pdf", col_loc[k], (int)E[i], j)).c_str());

      }// Generic 1d loop
    }// Eta Loop
  }// E loop

  // Writing out benchmarks
  //Tests
  std::string test_tag = "Barrel_emcal_pion_rejection";
  //TODO: Change test_tag to something else
  std:string detectorEle = "Barrel_emcal";
  
  for (int i = 0; i < etaTitle.size(); i++){
    etaTitle[i].erase(std::remove(etaTitle[i].begin(), etaTitle[i].end(), '#'), etaTitle[i].end());
    std::replace(etaTitle[i].begin(), etaTitle[i].end(), 'e', 'E');    
  }
  
  // E, Eta = 18, 2
  common_bench::Test pion_rejection_E18_Eta2{
    {{"name", fmt::format("{}_E{}_EtaBin{}", test_tag, (int)E[0], 2)},
     {"title", "Pion Rejection1"},
     {"description", fmt::format("Pion rejection with E = {}, and {}", (int)E[0], etaTitle[2])},
     {"quantity", "pi-/e-"},
     {"cut", currentCut},
     {"e- efficiency", std::to_string(effEle[0][2])},
     {"pi- efficiency", std::to_string(effPim[0][2])},
     {"target", std::to_string(suppression * maxRate[0][2])}
    }
  };
  suppression * maxRate[0][2] >= rejRatios[0][2] ? pion_rejection_E18_Eta2.pass(rejRatios[0][2]) : pion_rejection_E18_Eta2.fail(rejRatios[0][2]);   

  // E, Eta = 18, 3
  common_bench::Test pion_rejection_E18_Eta3{
    {{"name", fmt::format("{}_E{}_EtaBin{}", test_tag, (int)E[0], 3)},
     {"title", "Pion Rejection"},
     {"description", fmt::format("Pion rejection with E = {}, and {}", (int)E[0], etaTitle[3])},
     {"quantity", "pi-/e-"},
     {"cut", currentCut},
     {"e- efficiency", std::to_string(effEle[0][3])},
     {"pi- efficiency", std::to_string(effPim[0][3])},
     {"target", std::to_string(suppression * maxRate[0][3])}
   }
  };
  suppression * maxRate[0][3] >= rejRatios[0][3] ? pion_rejection_E18_Eta3.pass(rejRatios[0][3]) : pion_rejection_E18_Eta3.fail(rejRatios[0][3]);

  // E, Eta = 10, 2
  common_bench::Test pion_rejection_E10_Eta2{
    {{"name", fmt::format("{}_E{}_EtaBin{}", test_tag, (int)E[1], 2)},
     {"title", "Pion Rejection"},
     {"description", fmt::format("Pion rejection with E = {}, and {}", (int)E[1], etaTitle[2])},
     {"quantity", "pi-/e-"},
     {"cut", currentCut},
     {"e- efficiency", std::to_string(effEle[1][2])},
     {"pi- efficiency", std::to_string(effPim[1][2])},
     {"target", std::to_string(suppression * maxRate[1][2])}
    }
  };
  suppression * maxRate[1][2] >= rejRatios[1][2] ? pion_rejection_E10_Eta2.pass(rejRatios[1][2]) : pion_rejection_E10_Eta2.fail(rejRatios[1][2]);

  // E, Eta = 10, 3
  common_bench::Test pion_rejection_E10_Eta3{
    {{"name", fmt::format("{}_E{}_EtaBin{}", test_tag, (int)E[1], 3)},
     {"title", "Pion Rejection"},
     {"description", fmt::format("Pion rejection with E = {}, and {}", (int)E[1], etaTitle[3])},
     {"quantity", "pi-/e-"},
     {"e- efficiency", std::to_string(effEle[1][3])},
     {"pi- efficiency", std::to_string(effPim[1][3])},
     {"target", std::to_string(suppression * maxRate[1][3])}
    }
  };
  suppression * maxRate[1][3] >= rejRatios[1][3] ? pion_rejection_E10_Eta3.pass(rejRatios[1][3]) : pion_rejection_E10_Eta3.fail(rejRatios[1][3]);

  // E, Eta = 5, 2
  common_bench::Test pion_rejection_E5_Eta2{
    {{"name", fmt::format("{}_E{}_EtaBin{}", test_tag, (int)E[2], 2)},
     {"title", "Pion Rejection"},
     {"description", fmt::format("Pion rejection with E = {}, and {}", (int)E[2], etaTitle[2])},
     {"quantity", "pi-/e-"},
     {"cut", currentCut},
     {"e- efficiency", std::to_string(effEle[2][2])},
     {"pi- efficiency", std::to_string(effPim[2][2])},
     {"target", std::to_string(suppression * maxRate[2][2])}
    }
  };
  suppression * maxRate[2][2] >= rejRatios[2][2] ? pion_rejection_E5_Eta2.pass(rejRatios[2][2]) : pion_rejection_E5_Eta2.fail(rejRatios[2][2]);

  // E, Eta = 5, 3
  common_bench::Test pion_rejection_E5_Eta3{
    {{"name", fmt::format("{}_E{}_EtaBin{}", test_tag, (int)E[2], 3)},
     {"title", "Pion Rejection"},
     {"description", fmt::format("Pion rejection with E = {}, and {}", (int)E[2], etaTitle[3])},
     {"quantity", "pi-/e-"},
     {"cut", currentCut},
     {"e- efficiency", std::to_string(effEle[2][3])},
     {"pi- efficiency", std::to_string(effPim[2][3])},
     {"target", std::to_string(suppression * maxRate[2][3])}
    }
  };
  suppression * maxRate[2][3] >= rejRatios[2][3] ? pion_rejection_E5_Eta3.pass(rejRatios[2][3]) : pion_rejection_E5_Eta3.fail(rejRatios[2][3]);

  // Writing out all tests
  common_bench::write_test({pion_rejection_E18_Eta2, 
                         pion_rejection_E18_Eta3, 
                         pion_rejection_E10_Eta2, 
                         pion_rejection_E10_Eta3, 
                         pion_rejection_E5_Eta2, 
                         pion_rejection_E5_Eta3}, 
                         fmt::format("results/{}_pion_rej.json", detectorEle));
}
