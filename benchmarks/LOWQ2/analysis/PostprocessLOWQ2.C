#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <TStyle.h>

//----------------------------------------------------------------------
// Set global plot format variables
//----------------------------------------------------------------------
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
// Create acceptance plots
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

//----------------------------------------------------------------------
// Create rate plots
//---------------------------------------------------------------------- 
TH2* RatePlot(TDirectory* inputDir, int Module, int Layer, TString Tag="Quasi-Real") {
    
    TString histName = "AllHits/module"+std::to_string(Module)+"/layer"+std::to_string(Layer)+"/hxPixelyPixel";

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

//----------------------------------------------------------------------
// Create formatted acceptance plots on a canvas
//---------------------------------------------------------------------- 
void FormatAcceptancePlots(TDirectory* inputDir, TFile* outputFile, TString Tag="Quasi-Real") {

    //----------------------------------------------------------------------
    // E-Q2 acceptance plot
    //----------------------------------------------------------------------
    TString ReconEQ2Name = "Reconstructed-Track/hprimElogQ2";
    TString AllEQ2Name   = "Generated/hprimElogQ2";

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

    //----------------------------------------------------------------------
    // E-Theta acceptance plot
    //----------------------------------------------------------------------
    TString ReconEThetaName = "Reconstructed-Track/hprimEtheta";
    TString AllEThetaName   = "Generated/hprimEtheta";

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
    
    //----------------------------------------------------------------------
    // Integrated acceptance plot    
    //----------------------------------------------------------------------
    TString ReconIntegratedName = "IntegratedAcceptance";
    
    // Read in the plots from the input file
    TH1* IntegratedAcceptancePlot = (TH1*)inputDir->Get(ReconIntegratedName);
    // Check plots exist    
    if (!IntegratedAcceptancePlot) {
        std::cout << "Error: plot "<< ReconIntegratedName <<" not found in input file" << std::endl;
        return;
    }

    TCanvas* canvasIntegratedAcceptance = new TCanvas("IntegratedAcceptance", "IntegratedAcceptance", 1600, 1200);

    //Get xAxis title
    TString xAxisTitle = IntegratedAcceptancePlot->GetXaxis()->GetTitle();
    //Break up xAxis title into words by "_"
    TObjArray* xAxisTitleWords = xAxisTitle.Tokenize("_");
    //SetBinLabel for each bin
    for (int i=1; i<=xAxisTitleWords->GetEntries(); i++) {
        IntegratedAcceptancePlot->GetXaxis()->SetBinLabel(i, xAxisTitleWords->At(i-1)->GetName());
    }

    // Format the plot
    IntegratedAcceptancePlot->SetTitle("Integrated acceptance");
    IntegratedAcceptancePlot->GetXaxis()->SetTitle("");
    IntegratedAcceptancePlot->GetYaxis()->SetTitle("Acceptance");

    IntegratedAcceptancePlot->SetStats(0);
    IntegratedAcceptancePlot->Draw();

        
    // Get the number of bins in the histogram
    int nBins = IntegratedAcceptancePlot->GetNbinsX();

    // Create a TText object to draw the text
    TText t;
    t.SetTextAlign(22); // Center the text
    t.SetTextSize(0.02); // Set the text size

    // Loop over the bins
    for (int i = 1; i <= nBins; ++i) {
        // Get the bin content
        double binContent = IntegratedAcceptancePlot->GetBinContent(i);

        // Get the bin center
        double binCenter = IntegratedAcceptancePlot->GetBinCenter(i);

        // Draw the bin content at the bin center
        t.DrawText(binCenter, binContent+0.02, Form("%.3f", binContent));
    }

    // Update the canvas to show the changes
    gPad->Update();

    // Save the canvas to output file
    outputFile->WriteTObject(canvasIntegratedAcceptance);

    // Clean up
    delete canvasIntegratedAcceptance;

}

//----------------------------------------------------------------------
// Create formatted rate plots on a canvas
//----------------------------------------------------------------------
void FormatRatePlots(TDirectory* inputDir, TFile* outputFile, TString Tag="Quasi-Real") {

    TString histName = "AllHits/module2/layer0/hxPixelyPixel";

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
// Create formatted reconstruction plots on a canvas
//----------------------------------------------------------------------
void FormatReconstructionPlots(TDirectory* inputDir, TFile* outputFile, TString Tag="Quasi-Real") {
    
    TString EHistName     = "reconEVsPrimE";
    TString ThetaHistName = "reconThetaVsPrimTheta";
    TString PhiHistName   = "thetacut/reconPhiVsPrimPhi";

    TString EResName      = "ERes";
    TString ThetaResName  = "thetaRes";
    TString PhiResName    = "thetacut/phiRes";

    TCanvas* canvas = new TCanvas("ReconCanvas", "ReconCanvas", 2400, 1200);

    // Read in the plots from the input file
    TH1* EPlot      = (TH1*)inputDir->Get(EHistName);
    TH1* ThetaPlot = (TH1*)inputDir->Get(ThetaHistName);
    TH1* PhiPlot    = (TH1*)inputDir->Get(PhiHistName);

    TH1* EResPlot     = (TH1*)inputDir->Get(EResName);
    TH1* ThetaResPlot = (TH1*)inputDir->Get(ThetaResName);
    TH1* PhiResPlot   = (TH1*)inputDir->Get(PhiResName);
    
    // Check plots exist
    if (!ThetaPlot || !EPlot || !PhiPlot || !ThetaResPlot || !EResPlot || !PhiResPlot) {
        std::cout << "Error: plots "<< ThetaHistName <<", "<< EHistName <<", "<< PhiHistName <<", "<< ThetaResName <<", "<< EResName <<", "<< PhiResName <<" not found in input file" << std::endl;
        return;
    }
    
    // Draw the plots on the canvas 
    canvas->Divide(3,2);
    canvas->cd(1);
    gPad->SetLogz();
    
    EPlot->SetTitle("Reconstructed electron energy vs. primary electron energy");
    EPlot->GetXaxis()->SetTitle("E_{prim} [GeV]");
    EPlot->GetYaxis()->SetTitle("E_{recon} [GeV]");
    EPlot->GetZaxis()->SetTitle("Counts");

    EPlot->SetStats(0);
    EPlot->Draw("colz");

    canvas->cd(2);
    gPad->SetLogz();
    
    ThetaPlot->SetTitle("Reconstructed Theta vs. primary Theta");
    ThetaPlot->GetXaxis()->SetTitle("Theta_{prim} [mrad]");
    ThetaPlot->GetYaxis()->SetTitle("Theta_{recon} [mrad]");
    ThetaPlot->GetZaxis()->SetTitle("Counts");

    ThetaPlot->SetStats(0);
    ThetaPlot->Draw("colz");

    canvas->cd(3);
    gPad->SetLogz();

    PhiPlot->SetTitle("Reconstructed Phi vs. primary Phi");
    PhiPlot->GetXaxis()->SetTitle("Phi_{prim} [deg]");
    PhiPlot->GetYaxis()->SetTitle("Phi_{recon} [deg]");
    PhiPlot->GetZaxis()->SetTitle("Counts");

    PhiPlot->SetStats(0);
    PhiPlot->Draw("colz");

    canvas->cd(4);

    EResPlot->SetTitle("Reconstructed electron energy resolution");
    EResPlot->GetXaxis()->SetTitle("E_{recon} - E_{prim} / E_{prim}");
    EResPlot->GetYaxis()->SetTitle("Counts");

    EResPlot->SetStats(0);
    EResPlot->Draw();

    canvas->cd(5);

    ThetaResPlot->SetTitle("Reconstructed Theta resolution");
    ThetaResPlot->GetXaxis()->SetTitle("Theta_{recon} - Theta_{prim} [mrad]");
    ThetaResPlot->GetYaxis()->SetTitle("Counts");

    ThetaResPlot->SetStats(0);
    ThetaResPlot->Draw();

    canvas->cd(6);

    PhiResPlot->SetTitle("Reconstructed Phi resolution");
    PhiResPlot->GetXaxis()->SetTitle("Phi_{recon} - Phi_{prim} [deg]");
    PhiResPlot->GetYaxis()->SetTitle("Counts");

    PhiResPlot->SetStats(0);
    PhiResPlot->Draw();

    // Save the canvas to output file
    outputFile->WriteTObject(canvas);

    // Clean up
    delete canvas;

}

//----------------------------------------------------------------------
// This function is called by the benchmarking script maybe
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
    TDirectory* dirA = dir->GetDirectory("Acceptance");
    if (dirA) {
        FormatAcceptancePlots(dirA, outputFile,Tag);
    }

    //Check if SimDistributions directory exists
    TDirectory* dirS = dir->GetDirectory("Rates");
    if (dirS) {
        FormatRatePlots(dirS, outputFile,Tag);
    }

    //Check if ReconstructedDistributions directory exists
    TDirectory* dirR = dir->GetDirectory("Reconstruction");
    if (dirR) {
        FormatReconstructionPlots(dirR, outputFile,Tag);
    }
        
    inputFile->Close();
    outputFile->Close();
    
}

//----------------------------------------------------------------------
// Main function to create canvases
//----------------------------------------------------------------------
void PostprocessLOWQ2() {
    Postprocess("plots/LOWQ2QRRecon2.root", "plots/LOWQ2QR_FormattedPlots.root", "Quasi-Real");
    Postprocess("plots/LOWQ2BremsRecon2.root", "plots/LOWQ2Brems_FormattedPlots.root", "Bremsstrahlung");
}