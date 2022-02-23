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
#include "TGraphErrors.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;


void save_canvas(TCanvas* c, std::string label)
{
  c->SaveAs(fmt::format("results/energy_scan/{}.png",label).c_str());
  c->SaveAs(fmt::format("results/energy_scan/{}.pdf",label).c_str());
}

void save_canvas(TCanvas* c, std::string label, double E)
{
  std::string label_with_E = fmt::format("{}/{}", E, label); 
  save_canvas(c, label_with_E);
}
void save_canvas(TCanvas* c, std::string label, std::string E_label)
{
  std::string label_with_E = fmt::format("{}/{}", E_label, label); 
  save_canvas(c, label_with_E);
}
void save_canvas(TCanvas* c, std::string var_label, std::string E_label, std::string particle_label)
{
  std::string label_with_E = fmt::format("{}/hcal_barrel_{}_{}", E_label, particle_label, var_label); 
  save_canvas(c, label_with_E);
}

std::tuple <double, double, double, double> extract_sampling_fraction_parameters(std::string particle_label, std::string E_label)
{
  std::string input_fname = fmt::format("sim_output/energy_scan/{}/sim_hcal_barrel_{}.edm4hep.root", E_label, particle_label);
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

  // Define variables
  auto d1 = d0.Define("Ethr", Ethr, {"MCParticles"})
                .Define("nhits", nhits, {"EcalBarrelHits"})
                .Define("Esim", Esim, {"EcalBarrelHits"})
                .Define("fsam", fsam, {"Esim", "Ethr"});

  // Define Histograms
  auto hEthr = d1.Histo1D(
      {"hEthr", "Thrown Energy; Thrown Energy [GeV]; Events", 100, 0.0, 25.0},
      "Ethr");
  auto hNhits =
      d1.Histo1D({"hNhits", "Number of hits per events; Number of hits; Events",
                  100, 0.0, 2000.0},
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
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    save_canvas(c1, "Ethr", E_label, particle_label);  
  }

  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    save_canvas(c2, "nhits", E_label, particle_label);
  }

  {
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500);
    c3->SetLogy(1);
    auto h = hEsim->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    double up_fit = h->GetMean() + 5*h->GetStdDev();
    double down_fit = h->GetMean() - 5*h->GetStdDev();
    h->GetXaxis()->SetRangeUser(0.,up_fit);
    save_canvas(c3, "Esim", E_label, particle_label);
  }

  {
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    //c4->SetLogy(1);
    auto h = hfsam->DrawCopy();
    //h->GetYaxis()->SetTitleOffset(1.4);
    h->SetLineWidth(2);
    h->SetLineColor(kBlue); 
    double up_fit = h->GetMean() + 5*h->GetStdDev();
    double down_fit = h->GetMean() - 5*h->GetStdDev();
    h->Fit("gaus", "", "", down_fit, up_fit);
    h->GetXaxis()->SetRangeUser(0.,up_fit);
    TF1 *gaus = h->GetFunction("gaus");
    gaus->SetLineWidth(2);
    gaus->SetLineColor(kRed);    
    double mean = gaus->GetParameter(1);
    double sigma = gaus->GetParameter(2);
    double mean_err = gaus->GetParError(1);
    double sigma_err = gaus->GetParError(2);
    save_canvas(c4, "fsam", E_label, particle_label);
    return std::make_tuple(mean, sigma, mean_err, sigma_err);
  }
}

std::vector<std::string> read_scanned_energies(std::string input_energies_fname)
{
    std::vector<std::string> scanned_energies;
    std::string E_label;
    ifstream E_file (input_energies_fname);
    if (E_file.is_open())
    {
        while (E_file >> E_label)
        {
            scanned_energies.push_back(E_label);
        }
        E_file.close();
        return scanned_energies;
    }
    else
    {
        std::cout << "Unable to open file " << input_energies_fname << std::endl;
        abort();
    }
}

void hcal_barrel_energy_scan_analysis(std::string particle_label = "electron")
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

    auto scanned_energies = read_scanned_energies(fmt::format("sim_output/hcal_barrel_energy_scan_points_{}.txt", particle_label));

    TGraphErrors gr_fsam(scanned_energies.size()-1);
    TGraphErrors gr_fsam_res(scanned_energies.size()-1);

    for (const auto& E_label : scanned_energies) {
        auto [fsam, fsam_res, fsam_err, fsam_res_err] = extract_sampling_fraction_parameters(particle_label, E_label);
        auto E = std::stod(E_label);

        gr_fsam.SetPoint(gr_fsam.GetN(),E,100*fsam);
        gr_fsam.SetPointError(gr_fsam.GetN()-1,0., 100*fsam_err);
        gr_fsam_res.SetPoint(gr_fsam_res.GetN(),E,100.0*(fsam_res/fsam));
        auto fsam_res_rel_err = 100.0*(sqrt(pow((fsam_res_err/fsam),2)+pow((fsam_err*fsam_res)/(fsam*fsam),2)));
        gr_fsam_res.SetPointError(gr_fsam_res.GetN()-1,0.,fsam_res_rel_err);
    }

    TCanvas* c5 = new TCanvas("c5", "c5", 700, 500);
    c5->cd();
    gr_fsam.SetTitle("Sampling Fraction Scan;True Energy [GeV];Sampling Fraction [%]");
    gr_fsam.SetMarkerStyle(20);
    gr_fsam.Fit("pol0", "", "", 2., 20.);
    gr_fsam.Draw("APE");
    save_canvas(c5, fmt::format("hcal_barrel_{}_fsam_scan", particle_label));

    TCanvas* c6 = new TCanvas("c6", "c6", 700, 500);
    c6->cd();
    TF1* func_res = new TF1("func_res", "[0]/sqrt(x) + [1]", 0.25, 20.);
    func_res->SetLineWidth(2);
    func_res->SetLineColor(kRed);
    gr_fsam_res.SetTitle("Energy Resolution;True Energy [GeV];#Delta E/E [%]");
    gr_fsam_res.SetMarkerStyle(20);
    gr_fsam_res.Fit(func_res,"R");
    gr_fsam_res.Draw("APE");
    save_canvas(c6,fmt::format("hcal_barrel_{}_fsam_scan_res", particle_label));
}
