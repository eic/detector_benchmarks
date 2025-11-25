#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <vector>
#include <string>
#include <iostream>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TPad.h>
#include <TImage.h>
#include <cmath> 

void format_X_Q2_Acceptance(TString file5 = "lowq2-5x41/reconstruction_results.root",
                             TString file10 = "lowq2-10x100/reconstruction_results.root",
                             TString file18 = "lowq2-18x275/reconstruction_results.root",
                             TString outname = "X_Q2_Acceptance.png",
                             bool show_unresolvable_region = true) {

    // Global style tweaks
    gStyle->SetOptStat(0);

    bool add_epic_logo = false;

    // Define input files and their labels
    std::vector<TString> files = {
        file5,
        file10,
        file18
    };
    
    std::vector<TString> labels = {
        "5x41 GeV",
        "10x100 GeV",
        "18x275 GeV"
    };
    
    std::vector<int> colors = {kP6Blue, kP6Yellow, kP6Red};
    std::vector<int> fillStyles = {3004, 3005, 3006}; // Different hatched patterns for each energy
    
    TString histName           = "hLog10Q2_vs_log10x_acceptance";
    TString meanHistName       = "Q2_rel_res_vs_log10Q2_1";  // FitSlicesY creates this name from the relative resolution histogram
    TString resolutionHistName = "Q2_rel_res_vs_log10Q2_2";  // FitSlicesY creates this name from the relative resolution histogram
    
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "x-Q^{2} Acceptance with Resolution Regions", 1000, 700);

    TLegend *leg = new TLegend(0.15, 0.65, 0.45, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    
    bool firstDraw = true;
    double resolution_threshold = 1.0;
    
    printf("Looking for histogram: %s\\n", histName.Data());
    printf("Looking for resolution histogram: %s\\n", resolutionHistName.Data());
    
    // Loop over files to draw acceptance regions
    for (size_t i = 0; i < files.size(); i++) {
        TFile *f = TFile::Open(files[i]);
        if (!f || f->IsZombie()) {
            printf("Cannot open file: %s\n", files[i].Data());
            continue;
        }
        
        // Determine Q2 threshold for this specific file
        double q2_threshold = -10.0; // Default fallback
        TH1 *h_mean = (TH1*)f->Get(meanHistName);
        TH1 *h_res  = (TH1*)f->Get(resolutionHistName);
        if (h_res && h_mean) {
            for (int bin = 1; bin <= h_res->GetNbinsX(); bin++) {
                double log10q2 = h_res->GetBinCenter(bin);
                double rel_resolution = h_res->GetBinContent(bin); // This is now relative resolution directly
                double mean_offset = h_mean->GetBinContent(bin);

                double error = h_res->GetBinError(bin);
                
                // Skip bins with insufficient data or invalid fits
                if (error == 0) continue;                    // Require minimum entries
                if (rel_resolution <= 0) continue;           // Skip non-positive resolution
                if (rel_resolution > 2.0) continue;          // Skip unreasonably large relative resolution (>200%)

                // std::cout << "Bin " << bin << ": log10q2 = " << log10q2 << ", rel_resolution = " << rel_resolution << std::endl;
                
                // Look for where relative resolution exceeds threshold (poor resolution region)
                if (rel_resolution < resolution_threshold && std::abs(mean_offset) < resolution_threshold) {
                    q2_threshold = log10q2;
                    break;
                }
            }
        }
        
        printf("%s: Q2 threshold = log10(Q2) = %.2f (Q2 = %.2e GeV^2)\n", 
               labels[i].Data(), q2_threshold, std::pow(10.0, q2_threshold));
        
        TH2 *h = (TH2*)f->Get(histName);
        if (!h) {
            printf("Cannot find histogram %s in %s\n", histName.Data(), files[i].Data());
            f->Close();
            continue;
        }
        
        printf("Found histogram %s with %d entries, max value: %f\n", 
               histName.Data(), (int)h->GetEntries(), h->GetMaximum());
        
        h->SetContour(1);
        h->SetContourLevel(0, 0.005);
        h->SetLineColor(colors[i]);
        h->SetLineWidth(2);
        
        if (firstDraw) {
            h->Draw("CONT3");
            h->SetTitle(";log_{10}(x);log_{10}(Q^{2}) [GeV^{2}]");
            firstDraw = false;
        } else {
            h->Draw("CONT3 SAME");
        }
        // Add to legend
        leg->AddEntry(h, labels[i], "l");
        // break;

        if (show_unresolvable_region){
            TH2 *cut_hist = (TH2*)h->Clone("cut_hist");
            //Set range to below the Q2 cut
            cut_hist->GetYaxis()->SetRangeUser(-10.0, q2_threshold);
            
            // Set fill properties with unique pattern for each energy
            cut_hist->SetFillColorAlpha(colors[i],0.3);
            cut_hist->SetFillStyle(fillStyles[i]); // Use different pattern for each energy
            cut_hist->SetLineWidth(0); // No outline
            
            // Draw as filled contour
            cut_hist->SetContour(1);
            cut_hist->SetContourLevel(0, 0.005);
            cut_hist->Draw("CONT0 F SAME");

            leg->AddEntry(cut_hist, Form("%s Smallest resolvable bin", labels[i].Data()), "f");
        }


        
        // f->Close();
    }

    
    // Ensure we're back on the main canvas before drawing titles
    // c1->cd();
    
    double lm = gPad->GetLeftMargin();
    double rm = gPad->GetRightMargin();
    double tm = gPad->GetTopMargin();

    TLatex title;
    title.SetNDC();
    title.SetTextFont(72);         
    title.SetTextSize(0.035);
    title.SetTextAlign(12);
    title.DrawLatex(lm, 1.0 - tm + 0.3*tm, "ePIC Simulation 25.10.4");
    
    TLatex title2;
    title2.SetNDC();
    title2.SetTextFont(72);          
    title2.SetTextAlign(32);
    title2.SetTextSize(0.035);
    double x2 = 1.0 - rm;
    double y2 = 1.0 - tm + 0.3*tm;
    title2.DrawLatex(x2, y2, "Low Q^{2} Tagger: x-Q^{2} Acceptance");
    
    // Add description if showing unresolvable region
    // if (show_unresolvable_region) {
    //     TLatex desc;
    //     desc.SetNDC();
    //     desc.SetTextFont(42);
    //     desc.SetTextSize(0.025);
    //     desc.SetTextAlign(12);
    //     desc.DrawLatex(0.5, 0.85, "Shaded area: Poor Q^{2} Resolution Region");
    // }

    // Ensure we're on the main canvas
    // c1->cd();
    
    // Draw the legend
    leg->Draw();

   if(add_epic_logo)
     {
       // ===== Add ePIC logo to the figure if desired ======
       TImage *logo = TImage::Open("EPIC-logo_black_small.png");
       TPad *pad2 = new TPad("pad2", "Pad 2", 0.8, 0.8, 0.93, 0.93); // Create a new pad and then draw the image in it
       pad2->Draw();
       pad2->cd(); // Enter the new pad
       logo->Draw();
     }


    c1->Update();

    c1->SaveAs(outname);
}