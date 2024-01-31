#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <TStyle.h>

void SetStyle() {
    // Set global plot format variables
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.18);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetTitleSize(0.04, "XYZ");

    gStyle->SetTitleOffset(4.0, "Z");
    // gStyle->SetTitleOffset(1.0, "Y");
    
}

//----------------------------------------------------------------------
TH1* AcceptancePlot(TDirectory* inputDir, TString ReconHistName, TString AllHistName, TString Tag="Quasi-Real") {

    // Read in the plots from the input file
    TH1* ReconPlot = (TH1*)inputDir->Get(ReconHistName);
    TH1* AllPlot   = (TH1*)inputDir->Get(AllHistName);

    // Check plots exist
    if (!ReconPlot || !AllPlot) {
        std::cout << "Error: plots "<< ReconHistName <<" and/or "<< AllHistName <<" not found in input file" << std::endl;
        return nullptr;
    }

    //Divide the reconstructed plot by the all plot
    ReconPlot->Divide(AllPlot);

    return ReconPlot;

}

TH2* RatePlot(TDirectory* inputDir, int Module, int Layer, TString Tag="Quasi-Real") {
    
    TString histName = "rate/module"+std::to_string(Module)+"/layer"+std::to_string(Layer)+"/hxIDmodule"+std::to_string(Module)+"layer"+std::to_string(Layer)+"yIDmodule"+std::to_string(Module)+"layer"+std::to_string(Layer);

    // Read in the plots from the input file
    TH2* RatePlot = (TH2*)inputDir->Get(histName);
    // Check plots exist
    if (!RatePlot) {
        std::cout << "Error: plot "<< histName <<" not found in input file" << std::endl;
        return nullptr;
    }

    // Format the plot
    int rebin = 32;
    RatePlot->Rebin2D(rebin,rebin);
    RatePlot->Scale(1.0/(rebin*rebin));
    TString title = "Tagger module "+std::to_string(Module)+", layer "+std::to_string(Layer)+" - Mean "+Tag+" rate per "+std::to_string(rebin)+"x"+std::to_string(rebin)+" pixel group";
    // TString title = "Tagger module 2, layer 0 rate, integrated over "+std::to_string(rebin)+"x"+std::to_string(rebin)+" pixel group";
    RatePlot->SetTitle(title);
    RatePlot->GetXaxis()->SetTitle("x pixel");
    RatePlot->GetYaxis()->SetTitle("y pixel");
    RatePlot->GetZaxis()->SetTitle("Average Rate [Hz/55 #mum pixel]");

    RatePlot->SetStats(0);

    return RatePlot;
    
}


void FormatAcceptancePlots(TDirectory* inputDir, TFile* outputFile, TString Tag="Quasi-Real") {

    //E-Q2 acceptance plots
    TString ReconEQ2Name = "Reconstructed/hReconstructedprimElogQ2";
    TString AllEQ2Name   = "All/hAllprimElogQ2";

    TH1* AcceptEQ2 = AcceptancePlot(inputDir,ReconEQ2Name,AllEQ2Name);

    TCanvas* canvasEQ2 = new TCanvas("AcceptanceEQ2Canvas", "AcceptanceEQ2Canvas", 1600, 1200);

    // Draw the plots on the canvas 
    canvasEQ2->cd();

    TString title = "Acceptance as a function of scattered electron energy and reaction log_{10}(Q^{2})";
    AcceptEQ2->SetTitle(title);
    AcceptEQ2->GetXaxis()->SetTitle("E_{e'} [GeV]");
    AcceptEQ2->GetYaxis()->SetTitle("log_{10}(Q^{2}) [GeV^{2}]");
    AcceptEQ2->GetZaxis()->SetTitle("Acceptance");

    AcceptEQ2->SetStats(0);
    AcceptEQ2->Draw("colz");

    // Save the canvas to output file
    outputFile->WriteTObject(canvasEQ2);

    // Clean up
    delete canvasEQ2;

    //E-Theta acceptance plots
    TString ReconEThetaName = "Reconstructed/hReconstructedprimEtheta";
    TString AllEThetaName   = "All/hAllprimEtheta";

    TH1* AcceptETheta = AcceptancePlot(inputDir,ReconEThetaName,AllEThetaName);

    TCanvas* canvasETheta = new TCanvas("AcceptanceEThetaCanvas", "AcceptanceEThetaCanvas", 1600, 1200);

    // Draw the plots on the canvas
    canvasETheta->cd();

    title = "Acceptance as a function of scattered electron energy and theta";
    AcceptETheta->SetTitle(title);
    AcceptETheta->GetXaxis()->SetTitle("E_{e'} [GeV]");
    AcceptETheta->GetYaxis()->SetTitle("Theta_{e'} [mrad]");
    AcceptETheta->GetZaxis()->SetTitle("Acceptance");

    AcceptETheta->SetStats(0);
    AcceptETheta->Draw("colz");

    // Save the canvas to output file
    outputFile->WriteTObject(canvasETheta);

    // Clean up
    delete canvasETheta;

}

//----------------------------------------------------------------------

void FormatRatePlots(TDirectory* inputDir, TFile* outputFile, TString Tag="Quasi-Real") {

    TString histName = "rate/module2/layer0/hxIDmodule2layer0yIDmodule2layer0";

    TCanvas* canvas = new TCanvas("RateCanvas", "RateCanvas", 3200, 1200);
    canvas->Divide(2,1);

    // Read in the plots from the input file    
    TH2* RatePlot1_0 = RatePlot(inputDir,1,0,Tag);
    TH2* RatePlot2_0 = RatePlot(inputDir,2,0,Tag);
   
    // Draw the plots on the canvas 
    canvas->cd(1);
    gPad->SetLogz();
    RatePlot1_0->Draw("colz");

    canvas->cd(2);
    gPad->SetLogz();
    RatePlot2_0->Draw("colz");

    // Save the canvas to output file
    outputFile->WriteTObject(canvas);

    // Clean up
    delete canvas;

}   

//----------------------------------------------------------------------

void FormatReconstructionPlots(TDirectory* inputDir, TFile* outputFile, TString Tag="Quasi-Real") {
    
    TString Q2HistName = "reconQ2VsPrimQ2";
    TString EHistName  = "reconEVsPrimE";

    TCanvas* canvas = new TCanvas("ReconCanvas", "ReconCanvas", 2400, 1200);

    // Read in the plots from the input file
    TH1* Q2plot = (TH1*)inputDir->Get(Q2HistName);
    TH1* Eplot  = (TH1*)inputDir->Get(EHistName);

    // Check plots exist
    if (!Q2plot || !Eplot) {
        std::cout << "Error: plots "<< Q2HistName <<" and/or "<< EHistName <<" not found in input file" << std::endl;
        return;
    }
    
    // Draw the plots on the canvas 
    canvas->Divide(2,1);
    canvas->cd(1);
    gPad->SetLogz();
    
    Q2plot->SetTitle("Reconstructed Q^{2} vs. primary Q^{2}");
    Q2plot->GetXaxis()->SetTitle("log_{10}(Q^{2}_{prim}) [GeV^{2}]");
    Q2plot->GetYaxis()->SetTitle("log_{10}(Q^{2}_{recon}) [GeV^{2}]");
    Q2plot->GetZaxis()->SetTitle("Counts");

    Q2plot->SetStats(0);
    Q2plot->Draw("colz");

    canvas->cd(2);
    gPad->SetLogz();
    
    Eplot->SetTitle("Reconstructed electron energy vs. primary electron energy");
    Eplot->GetXaxis()->SetTitle("E_{prim} [GeV]");
    Eplot->GetYaxis()->SetTitle("E_{recon} [GeV]");
    Eplot->GetZaxis()->SetTitle("Counts");

    Eplot->SetStats(0);
    Eplot->Draw("colz");

    // Save the canvas to output file
    outputFile->WriteTObject(canvas);

    // Clean up
    delete canvas;

}

//----------------------------------------------------------------------
// This function is called by the benchmarking script
//----------------------------------------------------------------------
void Postprocess(TString inName="LOWQ2QRRates3.root", TString outName="LOWQ2Plots.root", TString Tag="Quasi-Real") {
    
    SetStyle();
    
    // Open the input file
    TFile* inputFile = TFile::Open(inName);
    if (!inputFile) {
        std::cout << "Error opening input file:" << inName << std::endl;
        return;
    }

    // Check the directory LOWQ2 exists and cd into it    
    TDirectory* dir = inputFile->GetDirectory("LOWQ2");
    if (!dir) {
        std::cout << "Error: directory LOWQ2 not found in input file" << std::endl;
        return;
    }

    // Open the output file
    TFile* outputFile = TFile::Open(outName, "RECREATE");

    // Format the plots

    //Check if AcceptanceDistributions directory exists
    TDirectory* dirA = dir->GetDirectory("AcceptanceDistributions");
    if (dirA) {
        FormatAcceptancePlots(dirA, outputFile,Tag);
    }

    //Check if SimDistributions directory exists
    TDirectory* dirS = dir->GetDirectory("SimDistributions");
    if (dirS) {
        FormatRatePlots(dirS, outputFile,Tag);
    }

    //Check if ReconstructedDistributions directory exists
    TDirectory* dirR = dir->GetDirectory("ReconstructedDistributions");
    if (dirR) {
        FormatReconstructionPlots(dirR, outputFile,Tag);
    }
        
    inputFile->Close();
    outputFile->Close();
    
}

//----------------------------------------------------------------------

void PostprocessLOWQ2() {
    Postprocess("LOWQ2QRRates3.root", "plots/LOWQ2QRPlots.root", "Quasi-Real");
    Postprocess("LOWQ2BremsRates3.root", "plots/LOWQ2BremsPlots.root", "Bremsstrahlung");
}