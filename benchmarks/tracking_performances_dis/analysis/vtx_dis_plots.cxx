#include <fstream>
#include <iostream>
#include <string>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

double StudentT(double *x, double *par);

void vtx_dis_plots(const std::string& config_name)
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

    TH2D* hgVr = (TH2D*) file->Get("recoVsMCTracks");
    TH2D* heff = (TH2D*) file->Get("recoVtxEff");
    
    TH2D* hg1 = (TH2D*) file->Get("genVtxYvsXHist");
    TH2D* hg2 = (TH2D*) file->Get("genVtxRvsZHist");
    
    TH2D* hr1 = (TH2D*) file->Get("recoVtxYvsX");
    TH2D* hr2 = (TH2D*) file->Get("recoVtxRvsZ");
    
    TH2D* hres1r = (TH2D*) file->Get("vtxResXvsGenTrk");
    TH2D* hres2r = (TH2D*) file->Get("vtxResYvsGenTrk");
    TH2D* hres3r = (TH2D*) file->Get("vtxResZvsGenTrk");
    
    TH2D* hres1g = (TH2D*) file->Get("vtxResXvsRecoTrk");
    TH2D* hres2g = (TH2D*) file->Get("vtxResYvsRecoTrk");
    TH2D* hres3g = (TH2D*) file->Get("vtxResZvsRecoTrk");
    
    TH1D* hng1 = (TH1D*) file->Get("numGenTracks");
    TH1D* hng2 = (TH1D*) file->Get("numGenTrkswithVtx");
    
    TH1D* hnr1 = (TH1D*) file->Get("numRecoTracks");
    TH1D* hnr2 = (TH1D*) file->Get("numRecoTrkswithVtx");
    
    //--------------------------------------------------------------------------------------------------------------------------------------------
  
    //Fit Function  
    TF1 *myfunction = new TF1("fit", StudentT, -2, 2, 4);
 
    //Fitting res v/s MC histograms
    myfunction->SetParameters(hres1g->GetMaximum(), 0, 0.05, 1);
    hres1g->FitSlicesY(myfunction, 0, -1, 10);

    myfunction->SetParameters(hres2g->GetMaximum(), 0, 0.05, 1);
    hres2g->FitSlicesY(myfunction, 0, -1, 10);
  
    myfunction->SetParameters(hres3g->GetMaximum(), 0, 0.05, 1);
    hres3g->FitSlicesY(myfunction, 0, -1, 10);
    
    //Fitting res v/s RC histograms
    myfunction->SetParameters(hres1r->GetMaximum(), 0, 0.05, 1);
    hres1r->FitSlicesY(myfunction, 0, -1, 10);
  
    myfunction->SetParameters(hres2r->GetMaximum(), 0, 0.05, 1);
    hres2r->FitSlicesY(myfunction, 0, -1, 10);
  
    myfunction->SetParameters(hres3r->GetMaximum(), 0, 0.05, 1);
    hres3r->FitSlicesY(myfunction, 0, -1, 10);
    
    //--------------------------------------------------------------------------------------------------------------------------------------------
    
    // Make new efficiency histograms
    TH1D *vtxEffVsGenTrkHist = new TH1D("vtxEffVsGenTrk","",31,-0.5,30.5);
    TH1D *vtxEffVsRecoTrkHist = new TH1D("vtxEffVsRecoTrk","",31,-0.5,30.5);
    
    // Fill efficiency v/s MC histogram
    for(int i=0; i<=hng1->GetNbinsX(); i++)
      {
	float neventsMC = hng1->GetBinContent(i);
	float nvtxevtsMC = hng2->GetBinContent(i);
		
	if(neventsMC != 0)
	{
		vtxEffVsGenTrkHist->SetBinContent(i, nvtxevtsMC/neventsMC);
		vtxEffVsGenTrkHist->SetBinError(i, sqrt((nvtxevtsMC+1)/(neventsMC+2)*((nvtxevtsMC+2)/(neventsMC+3)-(nvtxevtsMC+1)/(neventsMC+2))));
	}
      }
      
    vtxEffVsGenTrkHist->SetMarkerColor(kRed);
    vtxEffVsGenTrkHist->SetMarkerStyle(8); vtxEffVsGenTrkHist->SetMarkerSize(1.2);
    vtxEffVsGenTrkHist->SetTitle("Vertexing Efficiency vs MC Tracks");
    vtxEffVsGenTrkHist->GetXaxis()->SetTitle("N_{MC}");
    vtxEffVsGenTrkHist->GetYaxis()->SetTitle("Vertexing Efficiency");
    vtxEffVsGenTrkHist->GetYaxis()->SetRangeUser(0, 1.2);
    
    //Fill efficiency v/s RC histogram
    for(int i=0; i<=hnr1->GetNbinsX(); i++)
      {
	float neventsRC = hnr1->GetBinContent(i);
	float nvtxevtsRC = hnr2->GetBinContent(i);
		
	if(neventsRC != 0)
	{
		vtxEffVsRecoTrkHist->SetBinContent(i, nvtxevtsRC/neventsRC);
		vtxEffVsRecoTrkHist->SetBinError(i, sqrt((nvtxevtsRC+1)/(neventsRC+2)*((nvtxevtsRC+2)/(neventsRC+3)-(nvtxevtsRC+1)/(neventsRC+2))));
	}
      }
      
    vtxEffVsRecoTrkHist->SetMarkerColor(kRed);
    vtxEffVsRecoTrkHist->SetMarkerStyle(8); vtxEffVsRecoTrkHist->SetMarkerSize(1.2);
    vtxEffVsRecoTrkHist->SetTitle("Vertexing Efficiency vs RC Tracks");
    vtxEffVsRecoTrkHist->GetXaxis()->SetTitle("N_{RC}");
    vtxEffVsRecoTrkHist->GetYaxis()->SetTitle("Vertexing Efficiency");
    vtxEffVsRecoTrkHist->GetYaxis()->SetRangeUser(0, 1.2);
    
    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Make plots and save to PDF file

    // Update Style
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetOptStat(0);

    std::cout<<"Making plots..."<<std::endl;
  
    // MC vs RC associated particles
    TCanvas *c1 = new TCanvas("c1","MC vs Reco Tracks",800,600);
    c1->cd(1);
    auto func = new TF1("func","x",-1,40);
    func->SetLineStyle(2); func->SetLineColor(1);
    
    hgVr->Draw();
    func->Draw("SAMEL");

    //Vertexing Efficiency
    TCanvas *c2 = new TCanvas("c2","Vertexing Efficiency",800,600);
    int nevents(0);
    for(int i=0; i<=heff->GetNbinsX(); i++) nevents = nevents + heff->GetBinContent(i);
	
    c2->cd(1);
    heff->Draw("p");
    
    // MC vx versus vy Vertex
    TCanvas *c3 = new TCanvas("c3","MC Vertices",800,600);
    c3->cd(1);
    hg1->Draw();
  
    // MC vr versus vz Vertex
    TCanvas *c4 = new TCanvas("c4","MC Vertices",800,600);
    c4->cd(1);
    hg2->Draw();
    
    // Reconstructed Vertex vx versus vy
    TCanvas *c5 = new TCanvas("c5","Reconstructed Vertices",800,600);
    c5->cd(1);
    hr1->Draw();
  
    // Reconstructed Vertex vr versus vz
    TCanvas *c6 = new TCanvas("c6","Reconstructed Vertices",800,600);
    c6->cd(1);
    hr2->Draw();
  
    //Vertex Resolution vs MC Tracks
    TCanvas *c7 = new TCanvas("c7","VtxResX vs MCTrks",800,600);
    c7->cd(1);
    hres1g->Draw("COLZ");
    
    TCanvas *c8 = new TCanvas("c8","VtxResY vs MCTrks",800,600);
    c8->cd(1);
    hres2g->Draw("COLZ");
    
    TCanvas *c9 = new TCanvas("c9","VtxResZ vs MCTrks",800,600);
    c9->cd(1);
    hres3g->Draw("COLZ");
  
    //Vertex Resolution vs RC Tracks
    TCanvas *c10 = new TCanvas("c10","VtxResX vs RCTrks",800,600);
    c10->cd(1);
    hres1r->Draw("COLZ");
    
    TCanvas *c11 = new TCanvas("c11","VtxResY vs RCTrks",800,600);
    c11->cd(1);
    hres2r->Draw("COLZ");
    
    TCanvas *c12 = new TCanvas("c12","VtxResZ vs RCTrks",800,600);
    c12->cd(1);
    hres3r->Draw("COLZ");
  
    // Res Sigma v/s MC tracks
    TCanvas *c13 = new TCanvas("c13","Vertex Resolution vs MC Tracks",800,600);
    c13->cd(1);
  
    TH1D *resXsigmaM = (TH1D*)gDirectory->Get("vtxResXvsGenTrk_2");
    TH1D *resYsigmaM = (TH1D*)gDirectory->Get("vtxResYvsGenTrk_2");
    TH1D *resZsigmaM = (TH1D*)gDirectory->Get("vtxResZvsGenTrk_2");
    
    resXsigmaM->SetMarkerStyle(20); resYsigmaM->SetMarkerStyle(21); resZsigmaM->SetMarkerStyle(22);
    resXsigmaM->SetMarkerSize(1.2); resYsigmaM->SetMarkerSize(1.2); resZsigmaM->SetMarkerSize(1.2);
    resXsigmaM->SetMarkerColor(kBlue); resYsigmaM->SetMarkerColor(kRed); resZsigmaM->SetMarkerColor(kBlack);
    resXsigmaM->SetTitle("Vertex Resolution Sigma vs MC Tracks"); resYsigmaM->SetTitle("Vertex Resolution Sigma vs MC Tracks"); 
    resZsigmaM->SetTitle("Vertex Resolution Sigma vs MC Tracks");
    resZsigmaM->GetXaxis()->SetTitle("N_{MC}");
    resXsigmaM->GetYaxis()->SetTitle("#sigma (mm)"); resYsigmaM->GetYaxis()->SetTitle("#sigma (mm)"); resZsigmaM->GetYaxis()->SetTitle("#sigma (mm)");
    resXsigmaM->GetYaxis()->SetRangeUser(0, 1); resYsigmaM->GetYaxis()->SetRangeUser(0, 1); resZsigmaM->GetYaxis()->SetRangeUser(0, 1);
    
    resXsigmaM->Draw("P");
    resYsigmaM->Draw("PSAME");
    resZsigmaM->Draw("PSAME");
    
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
    legend1->AddEntry(resXsigmaM, "v_{x}", "lep");
    legend1->AddEntry(resYsigmaM, "v_{y}", "lep");
    legend1->AddEntry(resZsigmaM, "v_{z}", "lep");
    legend1->Draw();
  
    // Res Mean v/s MC tracks
    TCanvas *c14 = new TCanvas("c14","Vertex Resolution Mean vs MC Tracks",800,600);
    c14->cd(1);

    TH1D *resXmeanM = (TH1D*)gDirectory->Get("vtxResXvsGenTrk_1");
    TH1D *resYmeanM = (TH1D*)gDirectory->Get("vtxResYvsGenTrk_1");
    TH1D *resZmeanM = (TH1D*)gDirectory->Get("vtxResZvsGenTrk_1");
    
    resXmeanM->SetMarkerStyle(20); resYmeanM->SetMarkerStyle(21); resZmeanM->SetMarkerStyle(22);
    resXmeanM->SetMarkerSize(1.2); resYmeanM->SetMarkerSize(1.2); resZmeanM->SetMarkerSize(1.2);
    resXmeanM->SetMarkerColor(kBlue); resYmeanM->SetMarkerColor(kRed); resZmeanM->SetMarkerColor(kBlack);
    resXmeanM->SetTitle("Vertex Resolution Mean vs MC Tracks"); resYmeanM->SetTitle("Vertex Resolution Mean vs MC Tracks");
    resZmeanM->SetTitle("Vertex Resolution Mean vs MC Tracks");
    resZmeanM->GetXaxis()->SetTitle("N_{MC}");
    resXmeanM->GetYaxis()->SetTitle("#mu (mm)"); resYmeanM->GetYaxis()->SetTitle("#mu (mm)"); resZmeanM->GetYaxis()->SetTitle("#mu (mm)");
    resXmeanM->GetYaxis()->SetRangeUser(-1, 1); resYmeanM->GetYaxis()->SetRangeUser(-1, 1); resZmeanM->GetYaxis()->SetRangeUser(-1, 1);
    
    resXmeanM->Draw("P");
    resYmeanM->Draw("PSAME");
    resZmeanM->Draw("PSAME");
    
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
    legend2->AddEntry(resXmeanM, "v_{x}", "lep");
    legend2->AddEntry(resYmeanM, "v_{y}", "lep");
    legend2->AddEntry(resZmeanM, "v_{z}", "lep");
    legend2->Draw();
  
    //Res Sigma v/s RC tracks
    TCanvas *c15 = new TCanvas("c15","Vertex Resolution vs RC Tracks",800,600);
    c15->cd(1);  
    
    TH1D *resXsigmaR = (TH1D*)gDirectory->Get("vtxResXvsRecoTrk_2");
    TH1D *resYsigmaR = (TH1D*)gDirectory->Get("vtxResYvsRecoTrk_2");
    TH1D *resZsigmaR = (TH1D*)gDirectory->Get("vtxResZvsRecoTrk_2");
  
    resXsigmaR->SetMarkerStyle(20); resYsigmaR->SetMarkerStyle(21); resZsigmaR->SetMarkerStyle(22);
    resXsigmaR->SetMarkerSize(1.2); resYsigmaR->SetMarkerSize(1.2); resZsigmaR->SetMarkerSize(1.2);
    resXsigmaR->SetMarkerColor(kBlue); resYsigmaR->SetMarkerColor(kRed); resZsigmaR->SetMarkerColor(kBlack);
    resXsigmaR->SetTitle("Vertex Resolution Sigma vs RC Tracks"); resYsigmaR->SetTitle("Vertex Resolution Sigma vs RC Tracks"); 
    resZsigmaR->SetTitle("Vertex Resolution Sigma vs RC Tracks");
    resZsigmaR->GetXaxis()->SetTitle("N_{RC}");
    resXsigmaR->GetYaxis()->SetTitle("#sigma (mm)"); resYsigmaR->GetYaxis()->SetTitle("#sigma (mm)"); resZsigmaR->GetYaxis()->SetTitle("#sigma (mm)");
    resXsigmaR->GetYaxis()->SetRangeUser(0, 1); resYsigmaR->GetYaxis()->SetRangeUser(0, 1); resZsigmaR->GetYaxis()->SetRangeUser(0, 1);
    
    resXsigmaR->Draw("P");
    resYsigmaR->Draw("PSAME");
    resZsigmaR->Draw("PSAME");
    
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
    legend3->AddEntry(resXsigmaR, "v_{x}", "lep");
    legend3->AddEntry(resYsigmaR, "v_{y}", "lep");
    legend3->AddEntry(resZsigmaR, "v_{z}", "lep");
    legend3->Draw();
 
    //Res Mean v/s RC tracks
    TCanvas *c16 = new TCanvas("c16","Vertex Resolution Mean vs RC Tracks",800,600);
    c16->cd(1);
  
    TH1D *resXmeanR = (TH1D*)gDirectory->Get("vtxResXvsRecoTrk_1");
    TH1D *resYmeanR = (TH1D*)gDirectory->Get("vtxResYvsRecoTrk_1");
    TH1D *resZmeanR = (TH1D*)gDirectory->Get("vtxResZvsRecoTrk_1");
  
    resXmeanR->SetMarkerStyle(20); resYmeanR->SetMarkerStyle(21); resZmeanR->SetMarkerStyle(22);
    resXmeanR->SetMarkerSize(1.2); resYmeanR->SetMarkerSize(1.2); resZmeanR->SetMarkerSize(1.2);
    resXmeanR->SetMarkerColor(kBlue); resYmeanR->SetMarkerColor(kRed); resZmeanR->SetMarkerColor(kBlack);
    resXmeanR->SetTitle("Vertex Resolution Mean vs RC Tracks"); resYmeanR->SetTitle("Vertex Resolution Mean vs RC Tracks");
    resZmeanR->SetTitle("Vertex Resolution Mean vs RC Tracks");
    resZmeanR->GetXaxis()->SetTitle("N_{RC}");
    resXmeanR->GetYaxis()->SetTitle("#mu (mm)"); resYmeanR->GetYaxis()->SetTitle("#mu (mm)"); resZmeanR->GetYaxis()->SetTitle("#mu (mm)");
    resXmeanR->GetYaxis()->SetRangeUser(-1, 1); resYmeanR->GetYaxis()->SetRangeUser(-1, 1); resZmeanR->GetYaxis()->SetRangeUser(-1, 1);
    
    resXmeanR->Draw("P");
    resYmeanR->Draw("PSAME");
    resZmeanR->Draw("PSAME");
    
    TLegend *legend4 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
    legend4->AddEntry(resXmeanR, "v_{x}", "lep");
    legend4->AddEntry(resYmeanR, "v_{y}", "lep");
    legend4->AddEntry(resZmeanR, "v_{z}", "lep");
    legend4->Draw();
  
    //Vertexing Efficiency vs MC Tracks
    TCanvas *c17 = new TCanvas("c17","Vertexing Efficiency vs MC Tracks",800,600);
    c17->cd(1);

    vtxEffVsGenTrkHist->Draw("P");

    //Vertexing Efficiency vs RC Tracks
    TCanvas *c18 = new TCanvas("c18","Vertexing Efficiency vs RC Tracks",800,600);
    c18->cd(1);
    
    vtxEffVsRecoTrkHist->Draw("P");
    
    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Print plots to pdf file
    c1->Print(fmt::format("{}.pdf[", output_prefix).c_str());
    c1->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c2->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c3->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c4->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c5->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c6->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c7->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c8->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c9->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c10->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c11->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c12->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c13->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c14->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c15->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c16->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c17->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c18->Print(fmt::format("{}.pdf", output_prefix).c_str());
    c18->Print(fmt::format("{}.pdf]", output_prefix).c_str());
    
}

//Defining Fitting Function
double StudentT(double *x, double *par){
	double norm = par[0];
  	double mean = par[1];
  	double sigma = par[2];
  	double nu = par[3];

	double pi = 3.14;
  	double st = norm * (TMath::Gamma((nu+1.0)/2.0)/(TMath::Gamma(nu/2.0)*TMath::Sqrt(pi*nu)*sigma)) * TMath::Power( (1.0+TMath::Power((x[0]-mean)/sigma,2.0)/nu), (-(nu+1.0)/2.0) );
	return st;
}
