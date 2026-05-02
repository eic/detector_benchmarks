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

int acceptance_analysis(TString filename, string outname_pdf, string outname_png) 
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

    int nEtaBins = 100;
    int nPhiBins = 100;
    double etaMin = -5, etaMax = 0;

    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    TH2D* hEtaPhiAll = new TH2D("hEtaPhiAll", "All #pi- (status==1); #eta; #phi[rad]",
                                nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

    TH2D* hEtaPhiDetected = new TH2D("hEtaPhiDetected", "#pi- detected in nHCal; #eta; #phi[rad]",
                                     nEtaBins, etaMin, etaMax, nPhiBins, -TMath::Pi(), TMath::Pi());

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

        for (size_t i = 0; i < contrib_particle_idx.GetSize(); i++) {
            int idx = contrib_particle_idx[i];
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


