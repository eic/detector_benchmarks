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

constexpr double ETA_MIN = -5, ETA_MAX = 5;

inline bool inNHCal(double eta) {
  double a = std::abs(eta);
  return (a >= ETA_MIN && a <= ETA_MAX);
}


int light_by_light_analysis(const string& filename, string outname_pdf, string outname_png) {
    podio::ROOTReader reader;
    reader.openFile(filename);
    unsigned nEvents = reader.getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    TH2D* hEtaPt = new TH2D("hEtaPt", "Muon #eta vs p_{T};#eta;p_{T} [GeV]", 120, -6., 6., 100, 0., 5.);
    TH1D *hQ2[3], *hX[3], *hY[3], *hW[3];

    for (int i = 0; i < 3; i++) {
        hQ2[i] = new TH1D(Form("hQ2_%dmu", i), Form("Q^{2}, %d mu in nHCal;Q^{2} [GeV^{2}]", i), 100, 0., 100.);
        hX [i] = new TH1D(Form("hx_%dmu",  i), Form("x,    %d mu in nHCal;x", i), 100, 0., 1.);
        hY [i] = new TH1D(Form("hy_%dmu",  i), Form("y,    %d mu in nHCal;y", i), 100, 0., 1.);
        hW [i] = new TH1D(Form("hW_%dmu",  i), Form("W,    %d mu in nHCal;W [GeV]", i), 100, 0., 200.);
    }

    for (unsigned ev = 0; ev < nEvents; ev++) {
        auto frameData = reader.readNextEntry(podio::Category::Event);
        if (!frameData) continue;
        podio::Frame frame(std::move(frameData));

        auto& kinCol = frame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsTruth");
        auto& mcCol  = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        if (!kinCol.isValid() || !mcCol.isValid()) continue;

        if (kinCol.empty()) {
            std::cerr << "Warning: InclusiveKinematicsTruth is empty in event " << ev << std::endl;
            continue;
        }
        const auto& kin = kinCol.at(0);
        double Q2 = kin.getQ2(), x = kin.getX(), y = kin.getY(), W = kin.getW();
        //cout << "Event " << ev << ": MCParticles = " << mcCol.size() << ", Kin = " << kinCol.size() << endl;

        for (const auto& p : mcCol) {
            //if (p.getPDG() != 22 || p.getGeneratorStatus() != 2) continue;
            //if (p.getPDG() == 22) std::cout << "Foton znaleziony, genStat=" << p.getGeneratorStatus() << std::endl;
            auto daughters = p.getDaughters();
            if (daughters.size() != 2) continue;

            edm4hep::MCParticle m1, m2;
            for (const auto& d : mcCol) {
                if (d.getObjectID() == daughters[0].getObjectID()) m1 = d;
                if (d.getObjectID() == daughters[1].getObjectID()) m2 = d;
            }
            if (abs(m1.getPDG()) != 13 || abs(m2.getPDG()) != 13) continue;

            TLorentzVector v1(m1.getMomentum().x, m1.getMomentum().y, m1.getMomentum().z, m1.getEnergy());
            TLorentzVector v2(m2.getMomentum().x, m2.getMomentum().y, m2.getMomentum().z, m2.getEnergy());
            
            hEtaPt->Fill(v1.Eta(), v1.Pt());
            hEtaPt->Fill(v2.Eta(), v2.Pt());

            int nInNH = inNHCal(v1.Eta()) + inNHCal(v2.Eta());
            //std::cout << "Muon eta: " << v1.Eta() << ", " << v2.Eta() << " --> " << nInNH << " in nHCal" << std::endl;

            hQ2[nInNH]->Fill(Q2);
            hX [nInNH]->Fill(x);
            hY [nInNH]->Fill(y);
            hW [nInNH]->Fill(W);
        }
    }

    auto normalize = [](TH1D* h[3]) {
        int nBins = h[0]->GetNbinsX();
        array<TH1D*, 3> out;
        for (int k = 0; k < 3; ++k) {
            out[k] = (TH1D*)h[k]->Clone(Form("f_%s", h[k]->GetName()));
            for (int b = 1; b <= nBins; ++b) {
                double sum = h[0]->GetBinContent(b) + h[1]->GetBinContent(b) + h[2]->GetBinContent(b);
                if (sum > 0)
                    out[k]->SetBinContent(b, (h[k]->GetBinContent(b) / sum) * 100.);
            }
            int col = (k == 0 ? kBlack : k == 1 ? kBlue : kRed);
            out[k]->SetLineColor(col);
            out[k]->SetMarkerColor(col);
        }
        return out;
    };

    auto fQ2 = normalize(hQ2), fX = normalize(hX), fY = normalize(hY), fW = normalize(hW);

    TCanvas* canvas = new TCanvas("canvas", "muon analysis", 1600, 800);
    canvas->Divide(3, 2);
    canvas->cd(1); hEtaPt->Draw("COLZ");
    canvas->cd(2); for (int k = 0; k < 3; ++k) fQ2[k]->Draw(k ? "HIST SAME" : "HIST");
    canvas->cd(3); for (int k = 0; k < 3; ++k) fX[k]->Draw(k ? "HIST SAME" : "HIST");
    canvas->cd(4); for (int k = 0; k < 3; ++k) fY[k]->Draw(k ? "HIST SAME" : "HIST");
    canvas->cd(5); for (int k = 0; k < 3; ++k) fW[k]->Draw(k ? "HIST SAME" : "HIST");

    canvas->SaveAs(outname_pdf.c_str());
    canvas->SaveAs(outname_png.c_str());

    TCanvas* c2 = new TCanvas("geom_acceptance", "Geometric Acceptance vs Kinematic Variables", 1600, 800);
    c2->Divide(2, 2);


    vector<tuple<TH1D**, const char*>> toPlot = {
        {hQ2, "Q^{2} [GeV^{2}]"},
        {hX,  "x"},
        {hY,  "y"},
        {hW,  "W [GeV]"}
    };

    for (int i = 0; i < 4; ++i) {
        c2->cd(i+1);
        auto hset = get<0>(toPlot[i]);
        const char* xaxis = get<1>(toPlot[i]);

        TH1D* hFrac[4];
        for (int k = 0; k < 4; ++k) {
            hFrac[k] = (TH1D*)hset[0]->Clone(Form("frac_%d_%s", k, xaxis));
            hFrac[k]->Reset();
        }

        int nBins = hset[0]->GetNbinsX();
        for (int b = 1; b <= nBins; ++b) {
            double c0 = hset[0]->GetBinContent(b);
            double c1 = hset[1]->GetBinContent(b);
            double c2 = hset[2]->GetBinContent(b);
            double sum = c0 + c1 + c2;

            if (sum == 0) continue;

            hFrac[0]->SetBinContent(b, c0 / sum * 100.);
            hFrac[1]->SetBinContent(b, c1 / sum * 100.);
            hFrac[2]->SetBinContent(b, c2 / sum * 100.);
            hFrac[3]->SetBinContent(b, 100.);
        }

        const int colors[4] = {kBlue+2, kBlack, kRed+1, kMagenta+1};
        const char* labels[4] = {
            "0-in nHCal", "1-in nHCal", "2-in nHCal", "All calo."
        };

        for (int k = 0; k < 4; ++k) {
            hFrac[k]->SetLineColor(colors[k]);
            hFrac[k]->SetMarkerColor(colors[k]);
            hFrac[k]->SetLineWidth(k == 3 ? 3 : 2); 
            hFrac[k]->SetTitle(Form("geom. acc. vs %s", xaxis));
            hFrac[k]->GetXaxis()->SetTitle(xaxis);
            hFrac[k]->GetYaxis()->SetTitle("geom. acc. (%)");
            hFrac[k]->SetMinimum(0); hFrac[k]->SetMaximum(105);
            hFrac[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
        }

        TLine* l = new TLine(hFrac[0]->GetXaxis()->GetXmin(), 100.0, hFrac[0]->GetXaxis()->GetXmax(), 100.0);
        l->SetLineStyle(2);
        l->Draw("SAME");
    }

    c2->SaveAs("results/nhcal_light_by_light/geom_acceptance.pdf");
    c2->SaveAs("results/nhcal_light_by_light/geom_acceptance.png");

    return 0;
}