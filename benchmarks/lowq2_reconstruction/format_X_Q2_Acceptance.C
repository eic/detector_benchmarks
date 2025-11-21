#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <vector>
#include <string>
#include <TLatex.h> 

void X_Q2_Acceptance() {

    bool add_epic_logo = true;

    // Define input files and their labels
    std::vector<std::string> files = {
        "lowq2-5x41/reconstruction_results.root",
        "lowq2-10x100/reconstruction_results.root",
        "lowq2-18x275/reconstruction_results.root"
    };
    
    std::vector<std::string> labels = {
        "5x41 GeV",
        "10x100 GeV",
        "18x275 GeV"
    };
    
    std::vector<int> colors = {kP6Blue, kP6Yellow, kP6Red};
    
    std::string histName = "hLog10Q2_vs_log10x_acceptance";
    
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "x-Q^{2} Acceptance", 800, 600);

    
    TLegend *leg = new TLegend(0.15, 0.7, 0.35, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    
    bool firstDraw = true;
    
    // Loop over files
    for (size_t i = 0; i < files.size(); i++) {
        TFile *f = TFile::Open(files[i].c_str());
        if (!f || f->IsZombie()) {
            printf("Cannot open file: %s\n", files[i].c_str());
            continue;
        }
        
        TH2 *h = (TH2*)f->Get(histName.c_str());
        if (!h) {
            printf("Cannot find histogram %s in %s\n", histName.c_str(), files[i].c_str());
            f->Close();
            continue;
        }
        
        // Create contour at level > 0 to show acceptance region
        h->SetContour(1);
        h->SetContourLevel(0, 0.005);
        h->SetLineColor(colors[i]);
        h->SetLineWidth(2);
        // h->GetXaxis()->SetRangeUser(-10, 1);
        
        if (firstDraw) {
            h->Draw("CONT3");
            h->SetTitle(";log_{10}(x);log_{10}(Q^{2}) [GeV^{2}]");
            firstDraw = false;
        } else {
            h->Draw("CONT3 SAME");
        }
        
        leg->AddEntry(h, labels[i].c_str(), "l");
    }
    
    double lm = gPad->GetLeftMargin();
    double rm = gPad->GetRightMargin();
    double tm = gPad->GetTopMargin();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(72);         
    title.SetTextSize(0.035);       // adjust as desired
    title.SetTextAlign(13);         // left-top alignment
    title.DrawLatex(lm, 1.0 - tm + 0.5*tm, "ePIC Simulation 25.10.4");
    
    TLatex title2;
    title2.SetNDC();
    title2.SetTextFont(72);          // plain (change to 62 for bold)
    title2.SetTextAlign(33);         // right-top
    title2.SetTextSize(0.035);
    double x2 = 1.0 - rm;
    double y2 = 1.0 - tm + 0.5*tm;   // adjust 0.5*tm to move up/down
    title2.DrawLatex(x2, y2, "Low Q^{2} Tagger: x-Q^{2} Acceptance");


    leg->Draw();

//    if(add_epic_logo)
//      {
//        // ===== Add ePIC logo to the figure if desired ======
//        TImage *logo = TImage::Open("EPIC-logo_black_small.png");
//        TPad *pad2 = new TPad("pad2", "Pad 2", 0.8, 0.8, 0.93, 0.93); // Create a new pad and then draw the image in it
//        pad2->Draw();
//        pad2->cd(); // Enter the new pad
//        logo->Draw();
//      }


    c1->Update();

    c1->SaveAs("X_Q2_Acceptance.png");
}