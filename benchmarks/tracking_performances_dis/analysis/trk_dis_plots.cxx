#include <fstream>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

void trk_dis_plots(const std::string& config_name)
{
    // Read our configuration
    std::ifstream  config_file{config_name};
    nlohmann::json config;
    config_file >> config;
  
    const std::string hists_file    = config["hists_file"];
    const std::string detector      = config["detector"];
    const std::string output_prefix = config["output_prefix"];
    const int         ebeam         = config["ebeam"];
    const int         pbeam         = config["pbeam"];
    const int         Q2_min        = config["Min_Q2"];
    const int         nfiles        = config["nfiles"];
    
    fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
                "Plotting DIS tracking analysis...\n");
    fmt::print(" - Detector package: {}\n", detector);
    fmt::print(" - input file for histograms: {}\n", hists_file);
    fmt::print(" - output prefix for plots: {}\n", output_prefix);
    fmt::print(" - ebeam: {}\n", ebeam);
    fmt::print(" - pbeam: {}\n", pbeam);
    fmt::print(" - Minimum Q2: {}\n", Q2_min);
    fmt::print(" - nfiles: {}\n", nfiles);

    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Read file with histograms
    TFile* file = new TFile(hists_file.c_str());

    std::cout<<"Reading histograms..."<<std::endl;

    TH1D* h1a = (TH1D*) file->Get("h1a");
    TH1D* h1a1 = (TH1D*) file->Get("h1a1");
    TH1D* h1a2 = (TH1D*) file->Get("h1a2");

    TH1D* h1b = (TH1D*) file->Get("h1b");
    TH1D* h1b1 = (TH1D*) file->Get("h1b1");
    TH1D* h1b2 = (TH1D*) file->Get("h1b2");

    TH1D* h1c = (TH1D*) file->Get("h1c");
    TH1D* h1c1 = (TH1D*) file->Get("h1c1");
    TH1D* h1c2 = (TH1D*) file->Get("h1c2");

    TH1D* h2a = (TH1D*) file->Get("h2a");
    TH1D* h2b = (TH1D*) file->Get("h2b");

    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Make ratio histograms
    TH1 *h1rb1 = new TH1D("h1rb1","",100,-4,4); //Real-seeded tracks (Pt > 200 MeV/c cut)
    TH1 *h1rc1 = new TH1D("h1rc1","",100,-4,4); //Truth-seeded tracks (Pt > 200 MeV/c cut)
    TH1 *h1rb2 = new TH1D("h1rb2","",100,-4,4); //Real-seeded tracks (Pt > 500 MeV/c cut)
    TH1 *h1rc2 = new TH1D("h1rc2","",100,-4,4); //Truth-seeded tracks (Pt > 500 MeV/c cut)

    h1rb1 = (TH1*) h1b1->Clone("h1rb1");
    h1rb1->Divide(h1a1);
    h1rb1->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 200 MeV/c");
    h1rb1->GetXaxis()->SetTitle("#eta");h1rb1->GetXaxis()->CenterTitle();
    h1rb1->GetYaxis()->SetTitle("Ratio");h1rb1->GetYaxis()->CenterTitle();

    h1rc1 = (TH1*) h1c1->Clone("h1rc1");
    h1rc1->Divide(h1a1);
    h1rc1->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 200 MeV/c");
    h1rc1->GetXaxis()->SetTitle("#eta");h1rc1->GetXaxis()->CenterTitle();
    h1rc1->GetYaxis()->SetTitle("Ratio");h1rc1->GetYaxis()->CenterTitle();

    h1rb2 = (TH1*) h1b2->Clone("h1rb2");
    h1rb2->Divide(h1a2);
    h1rb2->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 500 MeV/c");
    h1rb2->GetXaxis()->SetTitle("#eta");h1rb2->GetXaxis()->CenterTitle();
    h1rb2->GetYaxis()->SetTitle("Ratio");h1rb2->GetYaxis()->CenterTitle();

    h1rc2 = (TH1*) h1c2->Clone("h1rc2");
    h1rc2->Divide(h1a2);
    h1rc2->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 500 MeV/c");
    h1rc2->GetXaxis()->SetTitle("#eta");h1rc2->GetXaxis()->CenterTitle();
    h1rc2->GetYaxis()->SetTitle("Ratio");h1rc2->GetYaxis()->CenterTitle();

    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Make plots and save to PDF file

    // Update Style
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetOptStat(0);

    std::cout<<"Making plots..."<<std::endl;

    //Generated charged particles
    TCanvas *c1a = new TCanvas("c1a");
    h1a->Draw();
    h1a1->Draw("same");
    h1a2->Draw("same");

    TLegend *leg1a = new TLegend(0.25,0.6,0.6,0.875);
    leg1a->SetBorderSize(0);leg1a->SetFillStyle(0);
    leg1a->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1a->AddEntry("h1a","All generated charged particles","l");
    leg1a->AddEntry("h1a1","+ P_{T} > 200 MeV/c","l");
    leg1a->AddEntry("h1a2","+ P_{T} > 500 MeV/c","l");
    leg1a->Draw();

    //Real-seeded tracks
    TCanvas *c1b = new TCanvas("c1b");
    h1b->Draw();
    h1b1->Draw("same");
    h1b2->Draw("same");

    TLegend *leg1b = new TLegend(0.25,0.6,0.6,0.875);
    leg1b->SetBorderSize(0);leg1b->SetFillStyle(0);
    leg1b->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1b->AddEntry("h1b","All real-seeded tracks","l");
    leg1b->AddEntry("h1b1","+ P_{T} > 200 MeV/c","l");
    leg1b->AddEntry("h1b2","+ P_{T} > 500 MeV/c","l");
    leg1b->Draw();

    //Truth-seeded tracks
    TCanvas *c1c = new TCanvas("c1c");
    h1c->Draw();
    h1c1->Draw("same");
    h1c2->Draw("same");

    TLegend *leg1c = new TLegend(0.25,0.6,0.6,0.875);
    leg1c->SetBorderSize(0);leg1c->SetFillStyle(0);
    leg1c->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1c->AddEntry("h1c","All truth-seeded tracks","l");
    leg1c->AddEntry("h1c1","+ P_{T} > 200 MeV/c","l");
    leg1c->AddEntry("h1c2","+ P_{T} > 500 MeV/c","l");
    leg1c->Draw();

    //Comparison 1
    TCanvas *c1d = new TCanvas("c1d");
    auto frame_d1 = c1d->DrawFrame(-4,0,4,1.2*h1a1->GetMaximum());
    frame_d1->GetXaxis()->SetTitle("#eta_{gen} or #eta_{rec}");frame_d1->GetXaxis()->CenterTitle();
    h1a1->Draw("same");
    h1b1->Draw("same");
    h1c1->Draw("same");

    TLegend *leg1d = new TLegend(0.15,0.675,0.6,0.875);
    leg1d->SetBorderSize(0);leg1d->SetFillStyle(0);
    leg1d->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1d->AddEntry("h1a1","Generated charged particles w/ P_{T} > 200 MeV/c","l");
    leg1d->AddEntry("h1b1","Real-seeded tracks w/ P_{T} > 200 MeV/c","l");
    leg1d->AddEntry("h1c1","Truth-seeded tracks w/ P_{T} > 200 MeV/c","l");
    leg1d->Draw();

    //Comparison 2a
    TCanvas *c1e = new TCanvas("c1e");
    auto frame_e1 = c1e->DrawFrame(-4,0,4,1.2*h1a1->GetMaximum());
    frame_e1->GetXaxis()->SetTitle("#eta_{gen} or #eta_{rec}");frame_e1->GetXaxis()->CenterTitle();
    h1a2->Draw("same");
    h1b2->Draw("P same");

    TLegend *leg1e = new TLegend(0.15,0.675,0.6,0.875);
    leg1e->SetBorderSize(0);leg1e->SetFillStyle(0);
    leg1e->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1e->AddEntry("h1a2","Generated charged particles w/ P_{T} > 500 MeV/c","fl");
    leg1e->AddEntry("h1b2","Real-seeded tracks w/ P_{T} > 500 MeV/c","p");
    leg1e->Draw();

    //Comparison 2b
    TCanvas *c1e1 = new TCanvas("c1e1");
    frame_e1->Draw();
    h1a2->Draw("same");
    h1c2->Draw("P same");

    TLegend *leg1e1 = new TLegend(0.15,0.675,0.6,0.875);
    leg1e1->SetBorderSize(0);leg1e1->SetFillStyle(0);
    leg1e1->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1e1->AddEntry("h1a2","Generated charged particles w/ P_{T} > 500 MeV/c","fl");
    leg1e1->AddEntry("h1c2","Truth-seeded tracks w/ P_{T} > 500 MeV/c","p");
    leg1e1->Draw();

    //Comparison 1 ratio
    TCanvas *c1f = new TCanvas("c1f");
    h1rb1->Draw();
    h1rc1->Draw("same");

    TLegend *leg1f = new TLegend(0.575,0.25,0.875,0.45);
    leg1f->SetBorderSize(0);leg1f->SetFillStyle(0);
    leg1f->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1f->AddEntry("h1rb1","Real-seeded tracking","l");
    leg1f->AddEntry("h1rc1","Truth-seeded tracking","l");
    leg1f->Draw();

    //Comparison 2 ratio
    TCanvas *c1g = new TCanvas("c1g");
    h1rb2->Draw();
    h1rc2->Draw("same");

    TLegend *leg1g = new TLegend(0.575,0.25,0.875,0.45);
    leg1g->SetBorderSize(0);leg1g->SetFillStyle(0);
    leg1g->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg1g->AddEntry("h1rb2","Real-seeded tracking","l");
    leg1g->AddEntry("h1rc2","Truth-seeded tracking","l");
    leg1g->Draw();

    //Hit-based associations -- real-seeded tracks
    TCanvas *c2a = new TCanvas("c2a");
    c2a->SetLogy();
    h2a->Draw();

    TLegend *leg2a = new TLegend(0.25,0.6,0.6,0.875);
    leg2a->SetBorderSize(0);leg2a->SetFillStyle(0);
    leg2a->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg2a->Draw();

    //Hit-based associations -- truth-seeded tracks
    TCanvas *c2b = new TCanvas("c2b");
    c2b->SetLogy();
    h2b->Draw();

    TLegend *leg2b = new TLegend(0.25,0.6,0.6,0.875);
    leg2b->SetBorderSize(0);leg2a->SetFillStyle(0);
    leg2b->SetHeader(Form("Pythia8: %dx%d GeV, Q^{2} > %d GeV^{2}",ebeam,pbeam,Q2_min));
    leg2b->Draw();

    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Print plots to pdf file
    c1a->Print(fmt::format("{}.pdf[", output_prefix).c_str());
    c1a->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1b->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1c->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1d->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1e->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1e1->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1f->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c1g->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c2a->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c2b->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c2b->Print(fmt::format("{}.pdf]", output_prefix).c_str());

}