#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>

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

    // // Set global plot format variables
    // int font = 42; // 
    // gStyle->SetTextFont(font); // Set the font for text
    // gStyle->SetLabelFont(font, "XYZ"); // Set the font for axis labels
    // gStyle->SetTitleFont(font, "XYZ"); // Set the font for titles

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
TH2* RatePlot(TDirectory* inputDir, int Module, int Layer, TString Tag="Quasi-Real", TString inTag="AllHits") {
    
    TString histName = inTag+"/module"+std::to_string(Module)+"/layer"+std::to_string(Layer)+"/hxPixelyPixelRate";

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
    AcceptETheta->GetYaxis()->SetTitle("#theta_{e'} [mrad]");
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
        t.DrawText(binCenter, binContent+0.02, Form("%.1f %s", binContent*100,"%"));
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

    TCanvas* canvas2 = new TCanvas("RateCanvasOverlay", "RateCanvasOverlay", 3200, 1200);
    canvas2->Divide(2,1);

    // Draw the plots on the canvas 
    canvas2->cd(1);
    gPad->SetLogz();
    RatePlot1_0->Draw("colz");

    //Add a grid ontop of the histogram showing 448x512 pixel chip boundaries
    double xMin = RatePlot1_0->GetXaxis()->GetXmin();
    double xMax = RatePlot1_0->GetXaxis()->GetXmax();
    double yMin = RatePlot1_0->GetYaxis()->GetXmin()+512/4;
    double yMax = RatePlot1_0->GetYaxis()->GetXmax()+512/4;
    std::vector<int> verticalLineWidths   = {3,1,3,1,3,1,3};
    std::vector<int> horizontalLineWidths = {3,1,2,1,2,1,2,1,2,1,2,1,3};
    std::vector<int> horizontalLineStyles = {1,7,1,7,1,7,1,7,1,7,1,7,1};
    std::vector<int> verticalLineColors   = {kRed,kBlue,kRed,kBlue,kRed,kBlue,kRed};
    std::vector<int> horizontalLineColors = {kRed,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kBlue,kRed};
    //Vertical lines
    for (int i=0; i<=6; i++) {
        TLine* line = new TLine(xMin+i*448, yMin, xMin+i*448, yMax);
        line->SetLineColor(verticalLineColors[i]);
        line->SetLineWidth(verticalLineWidths[i]);
        line->Draw();
    }
    //Horizontal lines
    for (int i=0; i<=12; i++) {
        TLine* line = new TLine(xMin, yMin+i*256, xMax, yMin+i*256);
        line->SetLineColor(horizontalLineColors[i]);
        line->SetLineWidth(horizontalLineWidths[i]);
        line->SetLineStyle(horizontalLineStyles[i]);
        line->Draw();
    }

    gPad->Update();

    canvas2->cd(2);
    gPad->SetLogz();
    RatePlot2_0->Draw("colz");

    //Add a grid on top of the histogram showing 448x512 pixel chip boundaries
    xMin = RatePlot2_0->GetXaxis()->GetXmin();
    xMax = RatePlot2_0->GetXaxis()->GetXmax();
    yMin = RatePlot2_0->GetYaxis()->GetXmin()+512/4;
    yMax = RatePlot2_0->GetYaxis()->GetXmax()+512/4;
    //Vertical lines
    for (int i = 0; i <= 6; i++) {
        TLine* line = new TLine(xMin+i*448, yMin, xMin+i*448, yMax);
        line->SetLineColor(verticalLineColors[i]);
        line->SetLineWidth(verticalLineWidths[i]);
        line->Draw();
    }
    //Horizontal lines
    for (int i = 0; i <= 12; i++) {
        TLine* line = new TLine(xMin, yMin+i*256, xMax, yMin+i*256);
        line->SetLineColor(horizontalLineColors[i]);
        line->SetLineWidth(horizontalLineWidths[i]);
        line->SetLineStyle(horizontalLineStyles[i]);
        line->Draw();
    }

    gPad->Update();

    // Save the canvas to output file
    outputFile->WriteTObject(canvas);

    // Clean up
    delete canvas;

    // Canvas showing primary and secondary hits
    // Todo: Neaten up
    TCanvas* canvas3 = new TCanvas("PrimarySecondary-RateCanvas", "PrimarySecondary-RateCanvas", 2400, 2400);
    canvas3->Divide(2,2);

    TH2* RatePlotPrimary1_0 = RatePlot(inputDir,1,0,Tag,"PrimaryHits");
    TH2* RatePlotPrimary2_0 = RatePlot(inputDir,2,0,Tag,"PrimaryHits");

    TH2* RatePlotSecondary1_0 = RatePlot(inputDir,1,0,Tag,"SecondaryHits");
    TH2* RatePlotSecondary2_0 = RatePlot(inputDir,2,0,Tag,"SecondaryHits");

    // Draw the plots on the canvas
    canvas3->cd(1);
    gPad->SetLogz();
    RatePlotPrimary1_0->Draw("colz");

    canvas3->cd(2);
    gPad->SetLogz();
    RatePlotPrimary2_0->Draw("colz");

    canvas3->cd(3);
    gPad->SetLogz();
    RatePlotSecondary1_0->Draw("colz");

    canvas3->cd(4);
    gPad->SetLogz();
    RatePlotSecondary2_0->Draw("colz");

    // Save the canvas to output file
    outputFile->WriteTObject(canvas3);

    // Clean up
    delete canvas3;


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
    TH2* EPlot      = (TH2*)inputDir->Get(EHistName);
    TH2* ThetaPlot  = (TH2*)inputDir->Get(ThetaHistName);
    TH2* PhiPlot    = (TH2*)inputDir->Get(PhiHistName);

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
    
    EPlot->SetTitle("Reconstructed electron energy vs. Primary electron energy");
    EPlot->GetXaxis()->SetTitle("E_{prim} [GeV]");
    EPlot->GetYaxis()->SetTitle("E_{recon} [GeV]");
    EPlot->GetZaxis()->SetTitle("Counts");
    EPlot->Rebin2D(2,2);

    EPlot->SetStats(0);
    EPlot->Draw("colz");

    canvas->cd(2);
    gPad->SetLogz();
    
    ThetaPlot->SetTitle("Reconstructed #theta vs. Primary #theta");
    ThetaPlot->GetXaxis()->SetTitle("#theta_{prim} [mrad]");
    ThetaPlot->GetYaxis()->SetTitle("#theta_{recon} [mrad]");
    ThetaPlot->GetZaxis()->SetTitle("Counts");
    ThetaPlot->Rebin2D(2,2);

    ThetaPlot->SetStats(0);
    ThetaPlot->Draw("colz");

    canvas->cd(3);
    gPad->SetLogz();

    PhiPlot->SetTitle("Reconstructed #varphi vs. Primary #varphi (#theta>1mrad)");
    PhiPlot->GetXaxis()->SetTitle("#phi_{prim} [deg]");
    PhiPlot->GetYaxis()->SetTitle("#phi_{recon} [deg]");
    PhiPlot->GetZaxis()->SetTitle("Counts");
    PhiPlot->Rebin2D(2,2);

    PhiPlot->SetStats(0);
    PhiPlot->Draw("colz");

    canvas->cd(4);

    EResPlot->SetTitle("Electron energy resolution");
    EResPlot->GetXaxis()->SetTitle("(E_{recon} - E_{prim}) / E_{prim}");
    EResPlot->GetYaxis()->SetTitle("Counts");

    EResPlot->SetStats(0);
    EResPlot->Draw();

    // Write fitted Gaussian standard deviation in pad
    //Fit Gaussian to histogram maximum bin +- 10 bins
    int maxBin = EResPlot->GetMaximumBin();
    int fitMin = maxBin-10;
    int fitMax = maxBin+10;    
    EResPlot->Fit("gaus", "Q", "", EResPlot->GetBinCenter(fitMin), EResPlot->GetBinCenter(fitMax));
    double EResStdDev = EResPlot->GetFunction("gaus")->GetParameter(2);    
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.03);
    latex->DrawLatexNDC(0.2, 0.8, Form("#sigma_{E} = %.2f %s", EResStdDev*100, "%"));

    // Remove fitted Gaussian from histogram
    EResPlot->GetListOfFunctions()->Remove(EResPlot->GetFunction("gaus"));

    canvas->cd(5);

    ThetaResPlot->SetTitle("Theta resolution");
    ThetaResPlot->GetXaxis()->SetTitle("#theta_{recon} - #theta_{prim} [mrad]");
    ThetaResPlot->GetYaxis()->SetTitle("Counts");

    ThetaResPlot->SetStats(0);
    ThetaResPlot->Draw();

    // Write fitted Gaussian standard deviation in pad
    //Fit Gaussian to histogram maximum bin +- 10 bins
    maxBin = ThetaResPlot->GetMaximumBin();
    fitMin = maxBin-10;
    fitMax = maxBin+10;
    ThetaResPlot->Fit("gaus", "Q", "", ThetaResPlot->GetBinCenter(fitMin), ThetaResPlot->GetBinCenter(fitMax));
    double ThetaResStdDev = ThetaResPlot->GetFunction("gaus")->GetParameter(2);
    latex->DrawLatexNDC(0.2, 0.8, Form("#sigma_{#theta} = %.2f mrad", ThetaResStdDev));

    // Remove fitted Gaussian from histogram
    ThetaResPlot->GetListOfFunctions()->Remove(ThetaResPlot->GetFunction("gaus"));

    canvas->cd(6);

    PhiResPlot->SetTitle("Phi resolution (#theta>1mrad)");
    PhiResPlot->GetXaxis()->SetTitle("#phi_{recon} - #phi_{prim} [deg]");
    PhiResPlot->GetYaxis()->SetTitle("Counts");

    PhiResPlot->SetStats(0);
    PhiResPlot->Draw();

    // Write fitted Gaussian standard deviation in pad
    //Fit Gaussian to histogram maximum bin +- 10 bins
    maxBin = PhiResPlot->GetMaximumBin();
    fitMin = maxBin-10;
    fitMax = maxBin+10;
    PhiResPlot->Fit("gaus", "Q", "", PhiResPlot->GetBinCenter(fitMin), PhiResPlot->GetBinCenter(fitMax));
    double PhiResStdDev = PhiResPlot->GetFunction("gaus")->GetParameter(2);
    latex->DrawLatexNDC(0.2, 0.8, Form("#sigma_{#phi} = %.1f deg", PhiResStdDev));

    // Remove fitted Gaussian from histogram
    PhiResPlot->GetListOfFunctions()->Remove(PhiResPlot->GetFunction("gaus"));

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
//     Postprocess("plots/LOWQ2QRRecon2.root", "plots/LOWQ2QR_FormattedPlots.root", "Quasi-Real");
     Postprocess("plots/LOWQ2BremsRecon3.root", "plots/LOWQ2Brems_FormattedPlots3.root", "Bremsstrahlung");
}