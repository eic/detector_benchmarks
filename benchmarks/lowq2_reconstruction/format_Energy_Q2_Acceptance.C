// Format the log10(Q2) vs E acceptance histogram from reconstructionAnalysis output
// Usage: root -l 'format_Energy_Q2_Acceptance.C("reconstruction_results.root", "Energy_Q2_Acceptance.png")'

#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaletteAxis.h>

void format_Energy_Q2_Acceptance(const char* infile = "reconstruction_results.root",
                                 const char* outname = "Energy_Q2_Acceptance.png") {
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
    
    TH2* h = (TH2*)f->Get("hLog10Q2_vs_E_acceptance");
    if (!h) {
        printf("Histogram hLog10Q2_vs_E_acceptance not found\n");
        f->Close();
        return;
    }
    
    // Canvas
    TCanvas* c = new TCanvas("c", "Energy vs Q2 Acceptance", 900, 700);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.16);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    
    // Draw
    h->SetTitle(";E_{e^{-}} [GeV];log_{10}(Q^{2}) [GeV^{2}]");
    h->SetMinimum(0.0);
    h->SetMaximum(1.0);
    h->Draw("COLZ");
    
    double lm = gPad->GetLeftMargin();
    double rm = gPad->GetRightMargin();
    double tm = gPad->GetTopMargin();

    // Title
    TLatex title;
    title.SetNDC();
    title.SetTextFont(72);
    // title.SetTextAlign(33);         // right-top
    title.SetTextSize(0.035);
    title.SetTextAlign(12);
    title.DrawLatex(lm, 1.0 - tm + 0.3*tm, "ePIC Simulation 25.10.4");

    TLatex title2;
    title2.SetNDC();
    title2.SetTextFont(72);          
    title2.SetTextAlign(32);         // right-top
    title2.SetTextSize(0.035);
    double x2 = 1.0 - rm;
    double y2 = 1.0 - tm + 0.3*tm;   // adjust 0.5*tm to move up/down
    title2.DrawLatex(x2, y2, "Low Q^{2} Tagger: Q^{2} vs E Acceptance");
    
    c->Update();

    // Adjust palette (color bar) to match main histogram vertical extent
    TPaletteAxis* pal = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
    if (pal) {
        double y1 = c->GetBottomMargin();          // lower bound at bottom margin
        double y2 = 1.0 - c->GetTopMargin();       // upper bound at top margin
        pal->SetY1NDC(y1);
        pal->SetY2NDC(y2);

        // Set a consistent width for the palette axis
        double x2 = 1.0 - c->GetRightMargin() + 0.05; // slight outward shift
        double x1 = x2 - 0.046;                        // width of color bar
        pal->SetX1NDC(x1);
        pal->SetX2NDC(x2);

        pal->SetLabelSize(0.035);
        pal->SetTitleSize(0.035);
        c->Modified();
        c->Update();
    }
    c->SaveAs(outname);
    
    f->Close();
}
