#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TProfile.h>
#include <TLegend.h>

#include "TROOT.h"

using namespace std;
using namespace TMath;

string changeExtension(const string& path, const string& new_ext)
{
    size_t pos = path.find_last_of('.');
    if (pos != string::npos)
        return path.substr(0, pos) + new_ext;
    return path + new_ext;
}

int basic_distribution_energy_resolution(const string &filename, string outname_png){
        
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.25);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);

    string outname_pdf = changeExtension(outname_png, ".pdf");

    TFile *file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) 
    {
        cerr << "Cannot open file: " << filename << endl;
        return 1;
    }

    TH2D *h_energyRes = (TH2D*)file->Get("h_energyRes");
    TProfile *p_energyRes = (TProfile*)file->Get("p_energyRes");
    TH1D *h_nHCal_hit_contrib_time   = (TH1D*)file->Get("h_nHCal_hit_contrib_time");
    TH1D *h_nHCal_hit_contrib_energy = (TH1D*)file->Get("h_nHCal_hit_contrib_energy");
    TH2D *h_nHCal_hit_contrib_energy_vs_time       = (TH2D*)file->Get("h_nHCal_hit_contrib_energy_vs_time");
    TH2D *h_nHCal_hit_contrib_energy_vs_telap      = (TH2D*)file->Get("h_nHCal_hit_contrib_energy_vs_telap");
    TH2D *h_nHCal_hit_contrib_energy_vs_time_total = (TH2D*)file->Get("h_nHCal_hit_contrib_energy_vs_time_total");

    if (!h_energyRes || !p_energyRes || !h_nHCal_hit_contrib_time || !h_nHCal_hit_contrib_energy ||
    !h_nHCal_hit_contrib_energy_vs_time || !h_nHCal_hit_contrib_energy_vs_telap ||
    !h_nHCal_hit_contrib_energy_vs_time_total) 
    {
        cerr << "Cannot retrieve histograms from file" << endl;
        file->Close();
        return 1;
    }

    TGraphErrors *g_resolution = new TGraphErrors();

    for (int i = 1; i <= p_energyRes->GetNbinsX(); i++) {
        double Ekin = p_energyRes->GetBinCenter(i);
        double mean = p_energyRes->GetBinContent(i);
        double rms = p_energyRes->GetBinError(i); 
        int entries = p_energyRes->GetBinEntries(i);
        
        if (entries < 10 || mean <= 0) continue;

        int n = g_resolution->GetN();
        g_resolution->SetPoint(n, Ekin, rms/mean);
        g_resolution->SetPointError(n, 0, rms/mean * sqrt(1.0/entries));
    }

    TF1 *fit = new TF1("fit", "[0]*pow(x, -0.5) + [1] + [2]*pow(x, 1)", 0.1, 6);
    fit->SetParameters(0.5, 0.05, 0.001);
    fit->SetParNames("a (stochastic)", "b (constant)", "c (noise)");
    g_resolution->Fit(fit, "R");

    TCanvas *c = new TCanvas("c", "Energy Resolution", 1600, 800);
    c->Divide(2, 1);
    c->cd(1); 
    
    h_energyRes->Draw("COLZ");
    p_energyRes->SetLineWidth(3); 
    p_energyRes->SetLineColor(kRed); 
    p_energyRes->SetMarkerColor(kRed);
    p_energyRes->Draw("SAME");

    c->cd(2);
    g_resolution->SetTitle("Energy Resolution;E_{kin} [GeV];#sigma/E");
    g_resolution->SetMarkerStyle(20);
    g_resolution->SetLineWidth(2); 
    g_resolution->SetLineColor(kBlue); 
    g_resolution->Draw("APE");
    fit->Draw("SAME");
    
    double par_a = fit->GetParameter(0);
    double par_b = fit->GetParameter(1);
    double par_c = fit->GetParameter(2);
    
    TLegend *leg = new TLegend(0.35, 0.65, 0.80, 0.80);
    leg->SetTextSize(0.025);
    leg->AddEntry(g_resolution, "Data", "p");
    leg->AddEntry(fit, Form("Fit: %.3f/#sqrt{E} + %.3f + %.3f#timesE", par_a, par_b, par_c), "l");
    leg->Draw();

    c->SaveAs(outname_png.c_str());
    c->SaveAs(outname_pdf.c_str());
    
    TCanvas *c_neutronThresholds = new TCanvas("c_neutronThresholds", "c_neutronThresholds", 1600, 800);
    c_neutronThresholds->Divide(3, 2);
    c_neutronThresholds->cd(1); h_nHCal_hit_contrib_time->Draw();
    c_neutronThresholds->cd(2); h_nHCal_hit_contrib_energy->Draw();
    c_neutronThresholds->cd(3); h_nHCal_hit_contrib_energy_vs_time->Draw("COLZ");
    c_neutronThresholds->cd(4); h_nHCal_hit_contrib_energy_vs_telap->Draw("COLZ");
    c_neutronThresholds->cd(5); h_nHCal_hit_contrib_energy_vs_time_total->Draw("COLZ");

    auto replaceStr = [](string& s, const string& from, const string& to) {
        size_t pos = s.find(from);
        if (pos != string::npos) s.replace(pos, from.size(), to);
    };

    string neutron_png = outname_png;
    string neutron_pdf = changeExtension(outname_png, ".pdf");

    replaceStr(neutron_png, "energy_resolution_analysis", "neutronThresholds");
    replaceStr(neutron_pdf, "energy_resolution_analysis", "neutronThresholds");

    c_neutronThresholds->SaveAs(neutron_png.c_str());
    c_neutronThresholds->SaveAs(neutron_pdf.c_str());
    delete c_neutronThresholds;

    file->Close();

    return 0;
}