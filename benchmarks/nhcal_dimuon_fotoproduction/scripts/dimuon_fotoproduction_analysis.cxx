#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "TLorentzVector.h"

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>
#include <TCanvas.h>

#include "TROOT.h"
#include "TRandom.h"
#include "TH3.h"

#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"

#include <podio/Frame.h>
#include <podio/CollectionBase.h>
#include "podio/ROOTReader.h"
#include "podio/CollectionIDTable.h"
#include "podio/ObjectID.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticleCollectionData.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleData.h"

#include "edm4hep/SimCalorimeterHitCollectionData.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitData.h"
#include "edm4hep/SimCalorimeterHit.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollectionData.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitObj.h"

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/Cluster.h"
#include "edm4eic/ClusterData.h"

#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitCollectionData.h"
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHitData.h"
#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitObj.h"

#include "edm4eic/InclusiveKinematicsCollection.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4hep/utils/vector_utils.h"

#include <edm4eic/vector_utils_legacy.h>
#include <edm4hep/Vector3f.h>

using namespace std;
using namespace ROOT;
using namespace TMath;
using namespace edm4hep;

constexpr double ETA_MIN = -4.14, ETA_MAX = -1.16;

inline bool inNHCal(double eta) {return (eta >= ETA_MIN && eta <= ETA_MAX);}

int dimuon_fotoproduction_analysis(const string& filename, string outname_pdf, string outname_png) {
    gStyle->SetOptStat(0);
    podio::ROOTReader reader;
    reader.openFile(filename);
    unsigned nEvents = reader.getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    TH2D* hEtaPt = new TH2D("hEtaPt", "Muon #eta vs p_{T};#eta;p_{T} [GeV]", 120, -6., 6., 100, 0., 5.);

    const int nb = 60;
    TH1D* hQ2_all = new TH1D("hQ2_all", "Q^{2}", nb, 1e-6, 1.0);
    TH1D* hX_all  = new TH1D("hX_all", "x",     nb, 1e-12, 1e-3);
    TH1D* hY_all  = new TH1D("hY_all", "y",     nb, 0., 1.);
    TH1D* hW_all  = new TH1D("hW_all", "W",     nb, 0., 160.);

    TH1D *hQ2_in[3], *hX_in[3], *hY_in[3], *hW_in[3];
    for (int i = 0; i < 3; ++i) {
        hQ2_in[i] = new TH1D(Form("hQ2_in_%d", i), "", nb, 1e-6, 1.0);
        hX_in[i]  = new TH1D(Form("hX_in_%d",  i), "", nb, 1e-12, 1e-3);
        hY_in[i]  = new TH1D(Form("hY_in_%d",  i), "", nb, 0., 1.);
        hW_in[i]  = new TH1D(Form("hW_in_%d",  i), "", nb, 0., 160.);
    }

    for (unsigned ev = 0; ev < nEvents; ev++) {
        auto frameData = reader.readNextEntry(podio::Category::Event);
        if (!frameData) continue;
        podio::Frame frame(std::move(frameData));

        auto& kinCol = frame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsTruth");
        auto& mcCol  = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        if (!kinCol.isValid() || !mcCol.isValid()) continue;
        if (kinCol.empty()) continue;

        const auto& kin = kinCol.at(0);
        double Q2 = kin.getQ2(), x = kin.getX(), y = kin.getY(), W = kin.getW();

        edm4hep::MCParticle m1, m2;
        for (const auto& p : mcCol) {
            if (abs(p.getPDG()) != 13 || p.getGeneratorStatus() == 0) continue;
            if (p.getPDG() == 13)  m1 = p;
            if (p.getPDG() == -13) m2 = p;
        }

        TLorentzVector v1(m1.getMomentum().x, m1.getMomentum().y, m1.getMomentum().z, m1.getEnergy());
        TLorentzVector v2(m2.getMomentum().x, m2.getMomentum().y, m2.getMomentum().z, m2.getEnergy());

        hEtaPt->Fill(v1.Eta(), v1.Pt());
        hEtaPt->Fill(v2.Eta(), v2.Pt());

        hQ2_all->Fill(Q2);
        hX_all->Fill(x);
        hY_all->Fill(y);
        hW_all->Fill(W);

        int nInNH = inNHCal(v1.Eta()) + inNHCal(v2.Eta());
        if (nInNH >= 0 && nInNH <= 2) {
            hQ2_in[nInNH]->Fill(Q2);
            hX_in[nInNH]->Fill(x);
            hY_in[nInNH]->Fill(y);
            hW_in[nInNH]->Fill(W);
        }
    }

    auto makeEffMultiGraph = [](TH1D* h_all, TH1D** h_in, const char* title, const char* xlabel, bool logx) {
        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.63, 0.7, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);

        Color_t colors[3] = {kBlue, kRed, kBlack};
        Style_t markers[3] = {20, 21, 22};

        for (int i = 0; i < 3; ++i) {
            std::vector<double> x_vals, y_vals;
            for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
                double all = h_all->GetBinContent(b);
                double sel = h_in[i]->GetBinContent(b);
                if (all < 5) continue;
                x_vals.push_back(h_all->GetBinCenter(b));
                y_vals.push_back(100. * sel / all);
            }
            TGraph* g = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
            g->SetMarkerColor(colors[i]);
            g->SetMarkerStyle(markers[i]);
            g->SetMarkerSize(1.0);
            g->SetLineColor(0);
            mg->Add(g, "P");
            leg->AddEntry(g, Form("%d mu-in nHCAL", i), "p");
        }

        std::vector<double> x_vals, y_vals;
        for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
            double all = h_all->GetBinContent(b);
            double sel = h_in[1]->GetBinContent(b) + h_in[2]->GetBinContent(b);
            if (all < 5) continue;
            x_vals.push_back(h_all->GetBinCenter(b));
            y_vals.push_back(100. * sel / all);
        }
        TGraph* gsum = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
        gsum->SetMarkerStyle(25); 
        gsum->SetMarkerColor(kMagenta + 1);
        gsum->SetFillColor(kMagenta+1);  
        gsum->SetLineColor(0);
        gsum->SetMarkerSize(1.0);
        mg->Add(gsum, "P");
        leg->AddEntry(gsum, "All eligible", "p");

        mg->SetTitle(Form("%s;%s;geom. acc. [%%]", title, xlabel));
        double xmin = mg->GetXaxis()->GetXmin();
        double xmax = mg->GetXaxis()->GetXmax();
        TLine* line = new TLine(xmin, 100.0, xmax, 100.0);
        line->SetLineStyle(2);           
        line->SetLineColor(kGray+1);     
        line->SetLineWidth(2);
        mg->GetListOfFunctions()->Add(line);
        return std::make_tuple(mg, leg, logx);
    };

    TCanvas* canvas = new TCanvas("canvas", "muon analysis", 1600, 800);
    canvas->Divide(3, 2);
    canvas->cd(1); hEtaPt->Draw("COLZ");

    auto [mg_q2, leg_q2, log_q2] = makeEffMultiGraph(hQ2_all, hQ2_in, "Geom. acc. vs Q^{2}", "Q^{2} [GeV^{2}]", true);
    canvas->cd(2); if (log_q2) gPad->SetLogx(); gPad->SetGrid(); mg_q2->Draw("A"); leg_q2->Draw();

    auto [mg_x, leg_x, log_x] = makeEffMultiGraph(hX_all, hX_in, "Geom. acc. vs x", "x", true);
    canvas->cd(3); if (log_x) gPad->SetLogx(); gPad->SetGrid(); mg_x->Draw("A"); leg_x->Draw();

    auto [mg_y, leg_y, log_y] = makeEffMultiGraph(hY_all, hY_in, "Geom. acc. vs y", "y", false);
    canvas->cd(4); gPad->SetGrid(); mg_y->Draw("A"); leg_y->Draw();

    auto [mg_w, leg_w, log_w] = makeEffMultiGraph(hW_all, hW_in, "Geom. acc. vs W", "W [GeV]", false);
    canvas->cd(5); gPad->SetGrid(); mg_w->Draw("A"); leg_w->Draw();

    canvas->SaveAs(outname_pdf.c_str());
    canvas->SaveAs(outname_png.c_str());

    return 0;
}