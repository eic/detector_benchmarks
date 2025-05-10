#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <set>
#include <iostream>
#include <TString.h>

using namespace std;

int nhcal_acceptance_analysis(TString filename, TString outname) 
{
    TChain *chain = new TChain("events");
    chain->Add(filename);
    

    TTreeReader reader(chain);

    
    TTreeReaderArray<int> mc_pdg(reader, "MCParticles.PDG");
    TTreeReaderArray<int> mc_genStatus(reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> mc_px(reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> mc_py(reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> mc_pz(reader, "MCParticles.momentum.z");

    TTreeReaderArray<int> contrib_particle_idx(reader, "_HcalEndcapNHitsContributions_particle.index");
    TTreeReaderArray<unsigned int> contrib_particle_cid(reader, "_HcalEndcapNHitsContributions_particle.collectionID");


    // Histograms eta-phi
    int nEtaBins = 200;
    int nPhiBins = 128;
    double etaMin = -5, etaMax = 0;

    TH2D* hEtaPhiAll = new TH2D("hEtaPhiAll", "All #pi- (status==1); #eta[1]; #phi[rad]",
                                nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

    TH2D* hEtaPhiDetected = new TH2D("hEtaPhiDetected", "#pi- detected in nHCal; #eta[1]; #phi[rad]",
                                     nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

    // TH2D* hAcceptance = new TH2D("hAcceptance", "Acceptance: pi- in nHCal / all; #eta[1]; #phi[rad]",
    //                              nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

    while (reader.Next()) 
    {
        map<int, pair<double, double>> pi_minus_eta_phi; 
        set<int> detected;
        for (size_t i = 0; i < mc_pdg.GetSize(); ++i) 
        {
            if (mc_pdg[i] == -211 && mc_genStatus[i] == 1) 
            {
                float px = mc_px[i];
                float py = mc_py[i];
                float pz = mc_pz[i];

                float p = sqrt(px * px + py * py + pz * pz);
                float eta = 0.5 * log((p + pz) / (p - pz + 1e-8));
                float phi = atan2(py, px);

                hEtaPhiAll->Fill(eta, phi);
                pi_minus_eta_phi[i] = make_pair(eta, phi);
            }
        }

        // Check, if pi- have contributions nHCal
        for (size_t i = 0; i < contrib_particle_idx.GetSize(); ++i) {
            int idx = contrib_particle_idx[i];
            if (pi_minus_eta_phi.count(idx)) {
                detected.insert(idx);
            }
        }

        // Filling histogram with detected pi-
        for (auto idx : detected) {
            auto [eta, phi] = pi_minus_eta_phi[idx];
            hEtaPhiDetected->Fill(eta, phi);
        }
    }

    // Calculating acceptance eta-phi
    // for (int etaBin = 1; etaBin <= hAcceptance->GetNbinsX(); ++etaBin) {
    //     for (int phiBin = 1; phiBin <= hAcceptance->GetNbinsY(); ++phiBin) {
    //         double all = hEtaPhiAll->GetBinContent(etaBin, phiBin);
    //         double detected = hEtaPhiDetected->GetBinContent(etaBin, phiBin);
    //         if (all > 0) {
    //             hAcceptance->SetBinContent(etaBin, phiBin, detected / all);
    //         } else {
    //             hAcceptance->SetBinContent(etaBin, phiBin, 0);
    //         }
    //     }
    // }

    TH2D* hAcceptance = (TH2D*)hEtaPhiAll->Clone("hAcceptance");
    hAcceptance->Divide(hEtaPhiDetected);
    hAcceptance->SetTitle("#pi- detected/All");

    TCanvas *canvas = new TCanvas("canvas", "pi- All", 1600, 600);
    canvas->Divide(3,1);
    canvas->cd(1);
    hEtaPhiAll->Draw("COLZ");
    canvas->cd(2);
    hEtaPhiDetected->Draw("COLZ");
    canvas->cd(3);
    hAcceptance->Draw("COLZ");

    canvas->Print(outname);
    return 0;
}


