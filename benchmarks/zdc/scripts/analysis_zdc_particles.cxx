////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHitData.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void analysis_zdc_particles(
  const std::string& input_fname = "sim_output/sim_zdc_uniform_neutrons.edm4hep.root",
  const std::string& results_path = "results/far_forward/zdc/"
) {
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

  // Hit position X
  auto hit_x_position = [&] (const std::vector<edm4hep::SimCalorimeterHitData>& hits) {
	std::vector<double> result;
	for(const auto& h: hits)
		result.push_back(h.position.x); //mm
  return result; 
  };

  // Hit position Y
  auto hit_y_position = [&] (const std::vector<edm4hep::SimCalorimeterHitData>& hits) {
        std::vector<double> result;
        for(const auto& h: hits)
		result.push_back(h.position.y); //mm
  return result;
  };

  // Hit position Z
  auto hit_z_position = [&] (const std::vector<edm4hep::SimCalorimeterHitData>& hits) {
        std::vector<double> result;
        for(const auto& h: hits)
		result.push_back(h.position.z); //mm
  return result;
  };

  // Define variables
  auto d1 = d0
    .Define("Ethr", Ethr, {"MCParticles"})
    .Define("nhits_Ecal", nhits, {"EcalFarForwardZDCHits"})
    .Define("Esim_Ecal", Esim, {"EcalFarForwardZDCHits"})
    .Define("fsam_Ecal", fsam, {"Esim_Ecal", "Ethr"})
    .Define("nhits_Hcal", nhits, {"HcalFarForwardZDCHits"})
    .Define("Esim_Hcal", Esim, {"HcalFarForwardZDCHits"})
    .Define("fsam_Hcal", fsam, {"Esim_Hcal", "Ethr"})
    .Define("hit_x_position", hit_x_position, {"HcalFarForwardZDCHits"})
	  .Define("hit_y_position", hit_y_position, {"HcalFarForwardZDCHits"})
  ;

  // Define Histograms
  auto hEthr = d1.Histo1D(
      {"hEthr", "Thrown Energy; Thrown Energy [GeV]; Events", 100, 0.0, 7.5},
      "Ethr");
  auto hNhits_Ecal =
      d1.Histo1D({"hNhits_Ecal", "Number of Ecal hits per events; Number of Ecal hits; Events",
                  100, 0.0, 2000.0},
                 "nhits_Ecal");
  auto hEsim_Ecal = d1.Histo1D(
      {"hEsim_Ecal", "Ecal Energy Deposit; Ecal Energy Deposit [GeV]; Events", 100, 0.0, 1.0},
      "Esim_Ecal");
  auto hfsam_Ecal = d1.Histo1D(
      {"hfsam_Ecal", "Ecal Sampling Fraction; Ecal Sampling Fraction; Events", 100, 0.0, 0.1},
      "fsam_Ecal");
  auto hNhits_Hcal =
      d1.Histo1D({"hNhits_Hcal", "Number of Hcal hits per events; Number of Hcal hits; Events",
                  100, 0.0, 2000.0},
                 "nhits_Hcal");
  auto hEsim_Hcal = d1.Histo1D(
      {"hEsim_Hcal", "Hcal Energy Deposit; Hcal Energy Deposit [GeV]; Events", 100, 0.0, 1.0},
      "Esim_Hcal");
  auto hfsam_Hcal = d1.Histo1D(
      {"hfsam_Hcal", "Hcal Sampling Fraction; Hcal Sampling Fraction; Events", 100, 0.0, 0.1},
      "fsam_Hcal");

  auto hXY_HitMap_Hcal = d1.Histo2D({"hXY_HitMap_Hcal", "hit position Y vs. X histogram; hit position X [mm]; hit position Y [mm]", 8,-1300,-500, 8, -400, 400}, "hit_x_position", "hit_y_position");


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
    c1->SaveAs(TString(results_path + "/zdc_Ethr.png"));
    c1->SaveAs(TString(results_path + "/zdc_Ethr.pdf"));
  }

  // Ecal
  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits_Ecal->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c2->SaveAs(TString(results_path + "/zdc_nhits_Ecal.png"));
    c2->SaveAs(TString(results_path + "/zdc_nhits_Ecal.pdf"));
  }

  {
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500);
    c3->SetLogy(1);
    auto h = hEsim_Ecal->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c3->SaveAs(TString(results_path + "/zdc_Esim_Ecal.png"));
    c3->SaveAs(TString(results_path + "/zdc_Esim_Ecal.pdf"));
  }

  std::cout << "derp4\n";
  {
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    c4->SetLogy(1);
    auto h = hfsam_Ecal->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    //h->Fit("gaus", "", "", 0.01, 0.1);
    //h->GetFunction("gaus")->SetLineWidth(2);
    //h->GetFunction("gaus")->SetLineColor(kRed);
    c4->SaveAs(TString(results_path + "/zdc_fsam_Ecal.png"));
    c4->SaveAs(TString(results_path + "/zdc_fsam_Ecal.pdf"));
  }

  // Hcal
  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits_Hcal->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c2->SaveAs(TString(results_path + "/zdc_nhits_Hcal.png"));
    c2->SaveAs(TString(results_path + "/zdc_nhits_Hcal.pdf"));
  }

  {
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500);
    c3->SetLogy(1);
    auto h = hEsim_Hcal->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    c3->SaveAs(TString(results_path + "/zdc_Esim_Hcal.png"));
    c3->SaveAs(TString(results_path + "/zdc_Esim_Hcal.pdf"));
  }

  std::cout << "derp4\n";
  {
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    c4->SetLogy(1);
    auto h = hfsam_Hcal->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    //h->Fit("gaus", "", "", 0.01, 0.1);
    //h->GetFunction("gaus")->SetLineWidth(2);
    //h->GetFunction("gaus")->SetLineColor(kRed);
    c4->SaveAs(TString(results_path + "/zdc_fsam_Hcal.png"));
    c4->SaveAs(TString(results_path + "/zdc_fsam_Hcal.pdf"));
  }
  {
    TCanvas* c5 = new TCanvas("c5", "c5", 600, 500);
    hXY_HitMap_Hcal->Draw("COLZ");
    c5->SaveAs(TString(results_path + "/zdc_xy_map_Hcal.png"));
    c5->SaveAs(TString(results_path + "/zdc_xy_map_Hcal.pdf"));
  }
}
