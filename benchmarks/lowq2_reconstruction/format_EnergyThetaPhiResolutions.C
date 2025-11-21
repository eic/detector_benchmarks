// Macro to build a nicely formatted canvas for energy, scattering angle and filtered phi resolutions
// Reads histograms previously produced by reconstructionAnalysis and saved in a ROOT file.
// Histograms expected:
//  2D:  E_vs_E, theta_vs_theta, phi_filtered_vs_phi
//  1D:  E_res, theta_diff, phi_filtered_diff
// Usage (ROOT): .x EnergyThetaPhiResolutions.C("reconstruction_results.root", "energy_theta_phi_filtered_nice.png")
// Or compile: root -l -b -q 'EnergyThetaPhiResolutions.C("reconstruction_results.root")'

#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TSystem.h"

// Helper to annotate mean and RMS for 1D histograms
void AnnotateStats(TH1* h, const char* units="") {
    if(!h) return;
    double mean = h->GetMean();
    double rms  = h->GetRMS();
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextSize(0.04);
    lat->SetTextColor(kGray+2);
    lat->DrawLatex(0.15, 0.85, Form("Mean = %.3g %s", mean, units));
    lat->DrawLatex(0.15, 0.78, Form("RMS = %.3g %s", rms, units));
}

// Apply common style to 2D histograms
void Style2D(TH2* h, const char* title) {
    if(!h) return;
    h->SetTitle(title);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
}

// Apply common style to 1D histograms
void Style1D(TH1* h, const char* title, Color_t col) {
    if(!h) return;
    h->SetTitle(title);
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(20);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
}

void EnergyThetaPhiResolutions(const char* inFile="reconstruction_results.root",
                               const char* outPng="energy_theta_phi_filtered_nice.png") {
    // Global style tweaks
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(60);
    gStyle->SetPalette(kViridis);

    std::unique_ptr<TFile> f(TFile::Open(inFile, "READ"));
    if(!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file: " << inFile << std::endl;
        return;
    }

    // Fetch histograms
    TH2* hE_vs_E          = dynamic_cast<TH2*>(f->Get("E_vs_E"));
    TH2* hTheta_vs_theta  = dynamic_cast<TH2*>(f->Get("theta_vs_theta"));
    TH2* hPhi_filtered_vs_phi = dynamic_cast<TH2*>(f->Get("phi_filtered_vs_phi"));
    TH1* hE_res           = dynamic_cast<TH1*>(f->Get("E_res"));
    TH1* hTheta_diff      = dynamic_cast<TH1*>(f->Get("theta_diff"));
    TH1* hPhi_filtered_diff = dynamic_cast<TH1*>(f->Get("phi_filtered_diff"));

    // Basic existence check
    if(!hE_vs_E || !hTheta_vs_theta || !hPhi_filtered_vs_phi || !hE_res || !hTheta_diff || !hPhi_filtered_diff) {
        std::cerr << "ERROR: Missing one or more required histograms in file: " << inFile << std::endl;
        return;
    }

    // Scale energy resolution to percentage (fraction -> %)
    if(hE_res) hE_res->Scale(100.0);

    // Style histograms
    Style2D(hE_vs_E, "Reco vs MC Energy; E_{e} reco [GeV]; E_{e} MC [GeV]");
    Style2D(hTheta_vs_theta, "Reco vs MC Scattering Angle; #theta_{reco} [mrad]; #theta_{MC} [mrad]");
    Style2D(hPhi_filtered_vs_phi, "Reco vs MC #phi ( #theta > 1 mrad ); #phi_{reco} [deg]; #phi_{MC} [deg]");

    Style1D(hE_res, "Energy Resolution; #DeltaE/E_{MC} [%]", kAzure+2);
    Style1D(hTheta_diff, "Scattering Angle Difference; #Delta#theta [mrad]", kOrange+7);
    Style1D(hPhi_filtered_diff, "#phi Difference ( #theta > 1 mrad ); #Delta#phi [deg]", kGreen+2);
    // Create canvas with extra top space (dedicated title pad) without altering 3x2 plot layout
    TCanvas *c = new TCanvas("energy_theta_phi_filtered_nice", "Energy/ScattAngle/Phi (Filtered) Resolutions", 2600, 1600);
    TPad *padTitle = new TPad("padTitle","Title",0.0,0.95,1.0,1.0); // ~7% height for title
    padTitle->SetFillColor(0); padTitle->SetBorderSize(0); padTitle->SetMargin(0,0,0,0); padTitle->Draw();
    TPad *padMain  = new TPad("padMain","Plots",0.0,0.0,1.0,0.95);
    padMain->SetFillColor(0); padMain->SetBorderSize(0); padMain->Draw();
    padMain->cd();
    padMain->Divide(3,2,0.005,0.005);

    // Adjust margins for each of the 6 plot pads (keep plotting area proportions)
    for(int pad=1; pad<=6; ++pad) {
        padMain->cd(pad);
        gPad->SetTopMargin(0.04);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.05);
    }

    // 2D plots with log z (suppress main titles, keep axis titles)
    padMain->cd(1); hE_vs_E->SetTitle(""); hE_vs_E->Draw("COL"); gPad->SetLogz();
    padMain->cd(2); hTheta_vs_theta->SetTitle(""); hTheta_vs_theta->Draw("COL"); gPad->SetLogz();
    padMain->cd(3); hPhi_filtered_vs_phi->SetTitle(""); hPhi_filtered_vs_phi->Draw("COL"); gPad->SetLogz();

    // 1D resolution/difference plots (suppress titles)
    padMain->cd(4); hE_res->SetTitle(""); hE_res->Draw("HIST"); hE_res->GetYaxis()->SetMaxDigits(3); AnnotateStats(hE_res, "%");
    padMain->cd(5); hTheta_diff->SetTitle(""); hTheta_diff->Draw("HIST"); hTheta_diff->GetYaxis()->SetMaxDigits(3); AnnotateStats(hTheta_diff, "mrad");
    padMain->cd(6); hPhi_filtered_diff->SetTitle(""); hPhi_filtered_diff->Draw("HIST"); hPhi_filtered_diff->GetYaxis()->SetMaxDigits(3); AnnotateStats(hPhi_filtered_diff, "deg");

    // Overall title in dedicated title pad (top-left)
    padTitle->cd();

    TLatex overall; 
    overall.SetNDC(); 
    overall.SetTextFont(72); 
    overall.SetTextSize(0.6); 
    overall.SetTextAlign(13);
    overall.DrawLatex(0.02, 0.64, "ePIC Simulation 25.10.4");

    TLatex details; 
    details.SetNDC(); 
    details.SetTextFont(72); 
    details.SetTextSize(0.6); 
    details.SetTextAlign(33); // right-top alignment
    details.DrawLatex(0.98, 0.64, "Low-Q^{2} Tagger Resolutions");

    // Improve aesthetics
    // for(int pad=1; pad<=6; ++pad) {
    //     c->cd(pad);
    //     gPad->SetGrid();
    // }

    c->Update();
    c->SaveAs(outPng);

}
