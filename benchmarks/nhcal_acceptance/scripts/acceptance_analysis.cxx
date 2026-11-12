#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <set>
#include <iostream>
#include <TString.h>
#include <TStyle.h>

#include <podio/Frame.h>
#include <podio/Reader.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/CaloHitContributionCollection.h>

R__LOAD_LIBRARY(podioIO)
R__LOAD_LIBRARY(edm4hep)

using namespace std;

int acceptance_analysis(TString filename, string outname_pdf, string outname_png) 
{
    auto reader = podio::makeReader(filename.Data());
    const unsigned nEvents = reader.getEntries("events");

    int nEtaBins = 100;
    int nPhiBins = 100;
    double etaMin = -5, etaMax = 0;

    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    TH2D* hEtaPhiAll = new TH2D("hEtaPhiAll", "All #pi- (status==1); #eta; #phi[rad]",
                                nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

    TH2D* hEtaPhiDetected = new TH2D("hEtaPhiDetected", "#pi- detected in nHCal; #eta; #phi[rad]",
                                     nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

    for (unsigned ev = 0; ev < nEvents; ++ev)
    {
        const auto frame = reader.readNextEvent();
        const auto& mcParticles = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        const auto& contributions = frame.get<edm4hep::CaloHitContributionCollection>("HcalEndcapNHitsContributions");

        map<int, pair<double, double>> pi_minus_eta_phi; 
        set<int> detected;
        for (size_t i = 0; i < mcParticles.size(); ++i) 
        {
            const auto mc = mcParticles[i];
            if (mc.getPDG() == -211 && mc.getGeneratorStatus() == 1) 
            {
                float px = mc.getMomentum().x;
                float py = mc.getMomentum().y;
                float pz = mc.getMomentum().z;

                float p = sqrt(px * px + py * py + pz * pz);
                float eta = 0.5 * log((p + pz) / (p - pz + 1e-8));
                float phi = atan2(py, px);

                hEtaPhiAll->Fill(eta, phi);
                pi_minus_eta_phi[i] = make_pair(eta, phi);
            }
        }

        for (const auto& contrib : contributions) {
            const auto particle = contrib.getParticle();
            if (!particle.isAvailable() || particle.getObjectID().collectionID != mcParticles.getID()) continue;
            int idx = particle.getObjectID().index;
            if (pi_minus_eta_phi.count(idx)) {
                detected.insert(idx);
            }
        }

        for (auto idx : detected) {
            auto [eta, phi] = pi_minus_eta_phi[idx];
            hEtaPhiDetected->Fill(eta, phi);
        }
    }

    TH2D* hAcceptance = (TH2D*)hEtaPhiAll->Clone("hAcceptance");
    hAcceptance->Divide(hEtaPhiDetected);
    hAcceptance->SetTitle("#pi- detected/All");
    hAcceptance->SetMinimum(0);
    hAcceptance->SetMaximum(1);

    TCanvas *canvas = new TCanvas("canvas", "pi- All", 1600, 600);
    canvas->Divide(3,1);
    canvas->cd(1);
    hEtaPhiAll->Draw("COLZ");
    canvas->cd(2);
    hEtaPhiDetected->Draw("COLZ");
    canvas->cd(3);
    hAcceptance->Draw("COLZ");

    canvas->SaveAs(outname_pdf.c_str());
    canvas->SaveAs(outname_png.c_str());

    return 0;
}
