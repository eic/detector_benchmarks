// Format the log10(Q^2) reconstruction vs MC histogram from reconstructionAnalysis output
// Usage: root -l 'format_Q2_Reconstruction.C("reconstruction_results.root", "log10Q2_Reconstruction.png")'

#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TLine.h>

void format_Q2_Reconstruction(const char* infile = "reconstruction_results.root",
                               const char* outname = "log10Q2_Reconstruction.png") {
    // Style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(60);
    
    // Load histogram
    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        printf("Cannot open %s\n", infile);
        return;
    }
    
    TH2* h = (TH2*)f->Get("log10Q2_vs_log10Q2");
    if (!h) {
        printf("Histogram log10Q2_vs_log10Q2 not found\n");
        f->Close();
        return;
    }
    
    // Canvas
    TCanvas* c = new TCanvas("c", "Q2 Reconstruction", 1200, 700);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.10);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetLogz();
    
    // Draw
    h->SetTitle(";log_{10}(Q^{2})_{MC} [GeV^{2}];log_{10}(Q^{2})_{reco} [GeV^{2}]");
    h->Draw("COLZ");
    
    
    double lm = gPad->GetLeftMargin();
    double rm = gPad->GetRightMargin();
    double tm = gPad->GetTopMargin();

    // Title
    TLatex title;
    title.SetNDC();
    title.SetTextFont(72);
    title.SetTextSize(0.035);
    title.SetTextAlign(12);
    title.DrawLatex(lm, 1.0 - tm + 0.3*tm, "ePIC Simulation 25.10.4");

    TLatex title2;
    title2.SetNDC();
    title2.SetTextFont(72);          
    title2.SetTextAlign(32);         // right-top
    title2.SetTextSize(0.035);
    double x2 = 1.0 - rm;
    double y2 = 1.0 - tm + 0.3*tm;   // adjust 0.3*tm to move up/down
    title2.DrawLatex(x2, y2, "Low Q^{2} Tagger: Q^{2} Reconstruction");
    
    c->Update();

    // Adjust palette (color bar) to match main histogram vertical extent
    TPaletteAxis* pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
    if (pal) {
        double y1 = c->GetBottomMargin();          // lower bound at bottom margin
        double y2 = 1.0 - c->GetTopMargin();       // upper bound at top margin
        pal->SetY1NDC(y1);
        pal->SetY2NDC(y2);

        // Set a consistent width for the palette axis
        double x2_pal = 1.0 - c->GetRightMargin() + 0.05; // slight outward shift
        double x1_pal = x2_pal - 0.046;                        // width of color bar
        pal->SetX1NDC(x1_pal);
        pal->SetX2NDC(x2_pal);

        pal->SetLabelSize(0.035);
        pal->SetTitleSize(0.035);
        c->Modified();
        c->Update();
    }
    c->SaveAs(outname);
    
    f->Close();
}
