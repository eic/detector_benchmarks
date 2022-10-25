////////////////////////////////////////
// Read reconstruction ROOT output file
// Plot variables
// M.Żurek
////////////////////////////////////////

#include "ROOT/RDataFrame.hxx"
#include <iostream>
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

#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;


void save_canvas(TCanvas* c, std::string label)
{
  c->SaveAs(fmt::format("results/energy_scan/{}.png",label).c_str());
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
  std::string label_with_E = fmt::format("{}/emcal_barrel_{}_{}", E_label, particle_label, var_label); 
  save_canvas(c, label_with_E);
}

std::tuple <double, double, double, double> extract_sampling_fraction_parameters(std::string particle_label, std::string E_label, dd4hep::Detector& detector)
{
  std::string input_fname = fmt::format("sim_output/energy_scan/{}/sim_emcal_barrel_{}.edm4hep.root", E_label, particle_label);
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

  // Energy deposited in layers 
  auto decoder = detector.readout("EcalBarrelImagingHits").idSpec().decoder();
  fmt::print("{}\n", decoder->fieldDescription());
  auto layer_index = decoder->index("layer");
  fmt::print(" layer index is {}.\n", layer_index);

  // Define variables
  auto d1 = d0.Define("Ethr", Ethr, {"MCParticles"})
                .Define("nhits", nhits, {"EcalBarrelImagingHits"})
                .Define("EsimImg", Esim, {"EcalBarrelImagingHits"})
                .Define("EsimScFi", Esim, {"EcalBarrelScFiHits"})
                .Define("Esim", "EsimImg+EsimScFi")
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
      {"hfsam", "Sampling Fraction; Sampling Fraction; Events", 800, 0.0, 0.2},
      "fsam");

  // Number of layers to read the edep from
  int nlayers = 7;

  TGraphErrors gr_no_edep(nlayers);
  TGraphErrors gr_edep_mean(nlayers);

  // Draw energy per layer
  TCanvas* clayer = new TCanvas("clayer", "clayer", 1250, 1000);
  clayer->Divide(4,5);

  for(int layer=1; layer<nlayers+1; layer++) {
    auto Esim_layer = [=](const std::vector<edm4hep::SimCalorimeterHitData>& evt) {
      auto layer_edep = 0.0;
      for (const auto& i: evt) {
        if (decoder->get(i.cellID, layer_index) == layer) {
          layer_edep += i.energy;
        }
      }
      return layer_edep;
    };

    auto d2 = d0.Define(fmt::format("EsimImg_layer_{}",layer).c_str(), Esim_layer, {"EcalBarrelImagingHits"})
                  .Define(fmt::format("EsimScFi_layer_{}",layer).c_str(), Esim_layer, {"EcalBarrelScFiHits"})
                  .Define(fmt::format("Esim_layer_{}",layer).c_str(), fmt::format("EsimImg_layer_{}+EsimScFi_layer_{}",layer,layer).c_str());
    std::cout << "Layer to process: " << layer << std::endl;

    auto hEsim_layer = d2.Histo1D({fmt::format("hEsim_layer_{}",layer).c_str(),fmt::format("Energy Deposit in layer {}; Energy Deposit [Gev]; Events",layer).c_str(), 200, 0.0, 0.04}, fmt::format("Esim_layer_{}",layer));

    clayer->cd(layer);
    gPad->SetLogy();
    auto h = hEsim_layer->DrawCopy();

    auto up_range = h->GetMean() + 3*h->GetStdDev();
    auto down_range = h->GetMean() - 3*h->GetStdDev();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    
    h->GetXaxis()->SetRange(h->GetXaxis()->GetBinUpEdge(1), up_range); // skip 0th bin 
    auto mean_layer = h->GetMean();
    auto rms_layer = h->GetStdDev();
    h->GetXaxis()->SetRange(); // reset the range
    h->GetXaxis()->SetRangeUser(0.,up_range);

    auto no_edep = (h->GetBinContent(1)/h->GetEntries())*100;  
    gr_no_edep.SetPoint(gr_no_edep.GetN(),layer,no_edep);
    gr_edep_mean.SetPoint(gr_edep_mean.GetN(),layer, mean_layer);
    gr_edep_mean.SetPointError(gr_edep_mean.GetN()-1,0, rms_layer);
  }
  save_canvas(clayer, "Esim_layer", E_label, particle_label);

  // Event Counts
  auto nevents_thrown = d1.Count();
  std::cout << "Number of Thrown Events: " << (*nevents_thrown) << "\n";

  // Draw Histograms and Graphs for every energy
  TCanvas* cnoedep = new TCanvas("c_noedep", "c_noedep", 700, 500);
  cnoedep->cd();
  gr_no_edep.SetTitle("% of events with no dE;Layer;Events with no dE [%]");
  gr_no_edep.SetMarkerStyle(20);
  gr_no_edep.Draw("AP");
  save_canvas(cnoedep, "Layer_nodep", E_label, particle_label);

  TCanvas* cmean = new TCanvas("c_mean", "c_mean", 700, 500);
  cmean->cd();
  gr_edep_mean.SetTitle("Mean and RMS of energy deposit;Layer;Mean dE [GeV]");
  gr_edep_mean.GetYaxis()->SetTitleOffset(1.4);
  gr_edep_mean.SetFillColor(4);
  gr_edep_mean.SetFillStyle(3010);
  gr_edep_mean.Draw("a3P");
  save_canvas(cmean, "Layer_Esim_mean", E_label, particle_label);

  {
    TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
    c1->SetLogy(1);
    auto h = hEthr->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    save_canvas(c1, "Ethr", E_label, particle_label);  
  }

  {
    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    c2->SetLogy(1);
    auto h = hNhits->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    save_canvas(c2, "nhits", E_label, particle_label);
  }

  {
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500);
    c3->SetLogy(1);
    auto h = hEsim->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue);
    double up_fit = h->GetMean() + 5*h->GetStdDev();
    double down_fit = h->GetMean() - 5*h->GetStdDev();
    h->GetXaxis()->SetRangeUser(0.,up_fit);
    save_canvas(c3, "Esim", E_label, particle_label);
  }

  {
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    auto h = hfsam->DrawCopy();
    h->SetLineWidth(2);
    h->SetLineColor(kBlue); 
    double up_fit = h->GetMean() + 5*h->GetStdDev();
    double down_fit = h->GetMean() - 5*h->GetStdDev();
    if(down_fit <=0 ) down_fit = h->GetXaxis()->GetBinUpEdge(1);
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

void emcal_barrel_energy_scan_analysis(std::string particle_label = "electron")
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

  auto scanned_energies = read_scanned_energies(fmt::format("sim_output/emcal_barrel_energy_scan_points_{}.txt", particle_label));
  
  //Take detector layers 
  std::string detector_path = "";
  std::string detector_name = "athena";
  if(std::getenv("DETECTOR_PATH")) {
    detector_path = std::getenv("DETECTOR_PATH");
  }
  if(std::getenv("DETECTOR_CONFIG")) {
    detector_name = std::getenv("DETECTOR_CONFIG");
  }

  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(fmt::format("{}/{}.xml", detector_path, detector_name));
    
  TGraphErrors gr_fsam(scanned_energies.size()-1);
  TGraphErrors gr_fsam_res(scanned_energies.size()-1);

  for (const auto& E_label : scanned_energies) {
      auto [fsam, fsam_res, fsam_err, fsam_res_err] = extract_sampling_fraction_parameters(particle_label, E_label, detector);
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
  save_canvas(c5, fmt::format("emcal_barrel_{}_fsam_scan", particle_label));

  TCanvas* c6 = new TCanvas("c6", "c6", 700, 500);
  c6->cd();
  TF1* func_res = new TF1("func_res", "[0]/sqrt(x) + [1]", 0.25, 20.);
  func_res->SetLineWidth(2);
  func_res->SetLineColor(kRed);
  gr_fsam_res.SetTitle("Energy Resolution;True Energy [GeV];#Delta E/E [%]");
  gr_fsam_res.SetMarkerStyle(20);
  gr_fsam_res.Fit(func_res,"R");
  gr_fsam_res.Draw("APE");
  save_canvas(c6,fmt::format("emcal_barrel_{}_fsam_scan_res", particle_label));

}
