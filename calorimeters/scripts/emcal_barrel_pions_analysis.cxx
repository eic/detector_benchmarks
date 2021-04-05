////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>

#include "dd4pod/Geant4ParticleCollection.h"
#include "dd4pod/CalorimeterHitCollection.h"
#include "dd4pod/TrackerHitCollection.h"
#include "eicd/RawCalorimeterHitCollection.h"
#include "eicd/RawCalorimeterHitData.h"
#include "eicd/CalorimeterHitCollection.h"
#include "eicd/CalorimeterHitData.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ClusterData.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1D.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

void emcal_barrel_pions_analysis(const char* input_fname = "sim_output/rec_emcal_barrel_uniform_pions.root")
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
    std::vector<double> result;
    result.push_back(TMath::Sqrt(input[2].psx*input[2].psx + input[2].psy*input[2].psy + input[2].psz*input[2].psz + input[2].mass*input[2].mass));
  return result;
  };

  // Reconstructed Energy [GeV] in XY merger
  auto ErecXY = [] (const std::vector<eic::CalorimeterHitData> & evt) {
    std::vector<double> result;
    auto total_eng = 0.0;
    for (const auto& i: evt)
      total_eng += i.energy;
    result.push_back(total_eng / 1.e+3);
    return result;
  };

  // Reconstructed Energy [GeV] in Z merger
  auto ErecZ = [] (const std::vector<eic::CalorimeterHitData> & evt) {
    std::vector<double> result;
    auto total_eng = 0.0;
    for (const auto& i: evt)
      total_eng += i.energy;
    result.push_back(total_eng / 1.e+3);
    return result;
  };

  // Number of Clusters
  auto ncluster = [] (const std::vector<eic::ClusterData>& evt) {return (int) evt.size(); };

  // Cluster Energy [GeV]
  auto Ecluster = [] (const std::vector<eic::ClusterData>& evt) {
    std::vector<double> result;
    for (const auto& i: evt)
      result.push_back(i.energy / 1.e+3);
    return result;
  };

  // Sampling fraction = Esampling / Ethrown
  auto fsam = [](const std::vector<double>& sampled, const std::vector<double>& thrown) {
    std::vector<double> result;
    for (const auto& E1 : thrown) {
      for (const auto& E2 : sampled)
        result.push_back(E2/E1);
    }
    return result;
  };

  // Define variables
  auto d1 = d0.Define("Ethr",      Ethr,       {"mcparticles2"})
	      .Define("ErecXY",    ErecXY,     {"RecoEcalBarrelHitsXY"})
	      .Define("ErecZ",     ErecZ,      {"RecoEcalBarrelHitsZ"})
	      .Define("ncluster",  ncluster,   {"EcalBarrelClusters"})
	      .Define("Ecluster",  Ecluster,   {"EcalBarrelClusters"})
	      .Define("fsam",      fsam,       {"Ecluster","Ethr"})
	      ;

  // Define Histograms
  auto hEthr     = d1.Histo1D({"hEthr",     "Thrown Energy; Thrown Energy [GeV]; Events",                            100, -0.5, 10.5}, "Ethr");
  auto hErecXY   = d1.Histo1D({"hErecXY",   "Reconstructed Energy in XY merger; Reconstructed Energy [GeV]; Events", 100, -0.5, 10.5}, "ErecXY");
  auto hErecZ    = d1.Histo1D({"hErecZ",    "Reconstructed Energy in Z merger; Reconstructed Energy [GeV]; Events",  100, -0.5, 10.5}, "ErecZ");
  auto hNCluster = d1.Histo1D({"hNCluster", "Number of Clusters; # of Clusters; Events",                              20, -0.5, 20.5}, "ncluster");
  auto hEcluster = d1.Histo1D({"hEcluster", "Cluster Energy; Cluster Energy [GeV]; Events",                          100, -0.5, 10.5}, "Ecluster");
  auto hfsam     = d1.Histo1D({"hfsam",     "Sampling Fraction; Sampling Fraction; Events",                          100,  0.0,  1.0}, "fsam");

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
  c1->SaveAs("results/emcal_pions_Ethr.png");
  c1->SaveAs("results/emcal_pions_Ethr.pdf");

  TCanvas *c2 = new TCanvas("c2", "c2", 700, 500);
  c2->SetLogy(1);
  hErecXY->GetYaxis()->SetTitleOffset(1.4);
  hErecXY->SetLineWidth(2);
  hErecXY->SetLineColor(kBlue);
  hErecXY->DrawClone();
  c2->SaveAs("results/emcal_pions_ErecXY.png");
  c2->SaveAs("results/emcal_pions_ErecXY.pdf");

  TCanvas *c3 = new TCanvas("c3", "c3", 700, 500);
  c3->SetLogy(1);
  hErecZ->GetYaxis()->SetTitleOffset(1.4);
  hErecZ->SetLineWidth(2);
  hErecZ->SetLineColor(kBlue);
  hErecZ->DrawClone();
  c3->SaveAs("results/emal_pions_ErecZ.png"); 
  c3->SaveAs("results/emal_pions_ErecZ.pdf");

  TCanvas *c4 = new TCanvas("c4", "c4", 700, 500);
  c4->SetLogy(1); 
  hNCluster->GetYaxis()->SetTitleOffset(1.6);
  hNCluster->SetLineWidth(2);
  hNCluster->SetLineColor(kBlue);
  hNCluster->DrawClone();
  c4->SaveAs("results/emcal_pions_ncluster.png");
  c4->SaveAs("results/emcal_pions_ncluster.pdf");

  TCanvas *c5 = new TCanvas("c5", "c5", 700, 500); 
  c5->SetLogy(1);
  hEcluster->GetYaxis()->SetTitleOffset(1.4);
  hEcluster->SetLineWidth(2);
  hEcluster->SetLineColor(kBlue);
  hEcluster->DrawClone();
  c5->SaveAs("results/emcal_pions_Ecluster.png");
  c5->SaveAs("results/emcal_pions_Ecluster.pdf");

  TCanvas *c6 = new TCanvas("c6", "c6", 700, 500);
  c6->SetLogy(1);
  hfsam->GetYaxis()->SetTitleOffset(1.4);
  hfsam->SetLineWidth(2);
  hfsam->SetLineColor(kBlue);
  hfsam->Fit("gaus","","",0.1,1.0);
  hfsam->GetFunction("gaus")->SetLineWidth(2);
  hfsam->GetFunction("gaus")->SetLineColor(kRed);
  hfsam->DrawClone();
  c6->SaveAs("results/emcal_pions_fsam.png");
  c6->SaveAs("results/emcal_pions_fsam.pdf");
}
