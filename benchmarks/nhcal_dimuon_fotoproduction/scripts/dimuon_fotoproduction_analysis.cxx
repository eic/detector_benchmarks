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
#include <unordered_map> 
#include <unordered_set>

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
#include "edm4eic/InclusiveKinematics.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4hep/utils/vector_utils.h"

#include <edm4eic/vector_utils_legacy.h>
#include <edm4hep/Vector3f.h>

#include "edm4eic/Track.h"
#include "edm4eic/TrackSegment.h"
#include "edm4eic/TrackSegmentCollectionData.h"
#include "edm4eic/TrackPoint.h"

#include "edm4eic/TrackSegmentCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/ReconstructedParticle.h"

#include "edm4eic/MCRecoCalorimeterHitAssociationCollection.h"
#include "edm4eic/MCRecoCalorimeterHitAssociation.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"
#include "edm4eic/MCRecoParticleAssociation.h"

using namespace std;
using namespace ROOT;
using namespace TMath;
using namespace edm4hep;

constexpr double ETA_MIN = -4.14, ETA_MAX = -1.16;

inline bool inNHCal(double eta) {return (eta >= ETA_MIN && eta <= ETA_MAX);}

auto addPrefixAfterSlash = [](const string& path,
                              const string& prefix) {
    const auto slash = path.find_last_of('/');
    if (slash == string::npos) return prefix + path;
    return path.substr(0, slash+1) + prefix + path.substr(slash+1);
};

TLine* makeLine(double x1, double y1, double x2, double y2) {
        TLine* l = new TLine(x1, y1, x2, y2);
        l->SetLineColor(kRed);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw("same");
        return l;
    }

inline double dist3(const edm4hep::Vector3f& a, const edm4hep::Vector3f& b) {
    const double dx = double(a.x) - double(b.x);
    const double dy = double(a.y) - double(b.y);
    const double dz = double(a.z) - double(b.z);
    return sqrt(dx*dx + dy*dy + dz*dz);
}

inline tuple<TMultiGraph*, TLegend*, bool>
makeEffMultiGraph(TH1D* h_all,
                  TH1D* const h_in[3],
                  const char* title,
                  const char* xlabel,
                  bool logx)
{
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg   = new TLegend(0.63, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    Color_t colors[3]  = {kBlue, kRed, kBlack};
    Style_t markers[3] = {20, 21, 22};

    int added = 0;

    for (int i = 0; i < 3; ++i) {
        if (!h_all || !h_in[i]) continue;

        vector<double> x_vals, y_vals;
        x_vals.reserve(h_all->GetNbinsX());
        y_vals.reserve(h_all->GetNbinsX());

        for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
            double all = h_all->GetBinContent(b);
            if (all < 2) continue;
            double sel = h_in[i]->GetBinContent(b);

            x_vals.push_back(h_all->GetBinCenter(b));
            y_vals.push_back(100. * sel / all);
        }

        if (!x_vals.empty()) {
            TGraph* g = new TGraph((int)x_vals.size(), x_vals.data(), y_vals.data());
            g->SetMarkerColor(colors[i]);
            g->SetMarkerStyle(markers[i]);
            g->SetMarkerSize(1.0);
            g->SetLineColor(0);
            mg->Add(g, "P");
            leg->AddEntry(g, Form("%d mu-in nHCAL", i), "p");
            ++added;
        }
    }

    if (h_all && h_in[0] && h_in[1] && h_in[2]) {
        vector<double> x_vals, y_vals;
        x_vals.reserve(h_all->GetNbinsX());
        y_vals.reserve(h_all->GetNbinsX());

        for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
            double all = h_all->GetBinContent(b);
            if (all < 2) continue;

            double sel = h_in[0]->GetBinContent(b)
                       + h_in[1]->GetBinContent(b)
                       + h_in[2]->GetBinContent(b);

            x_vals.push_back(h_all->GetBinCenter(b));
            y_vals.push_back(100. * sel / all);
        }

        if (!x_vals.empty()) {
            TGraph* gsum = new TGraph((int)x_vals.size(), x_vals.data(), y_vals.data());
            gsum->SetMarkerStyle(25);
            gsum->SetMarkerColor(kMagenta + 1);
            gsum->SetFillColor(kMagenta + 1);
            gsum->SetLineColor(0);
            gsum->SetMarkerSize(1.0);
            mg->Add(gsum, "P");
            leg->AddEntry(gsum, "All eligible", "p");
            ++added;
        }
    }

    mg->SetTitle(Form("%s;%s;geom. acc. [%%]", title ? title : "", xlabel ? xlabel : ""));

    if (added == 0 && h_all) {
        TH1D* hframe = new TH1D("hframe_tmp", "", 10,
                                h_all->GetXaxis()->GetXmin(),
                                h_all->GetXaxis()->GetXmax());
        hframe->SetMinimum(0.0);
        hframe->SetMaximum(110.0);
        hframe->GetXaxis()->SetTitle(xlabel);
        hframe->GetYaxis()->SetTitle("geom. acc. [%]");
        hframe->Draw();
    }

    return make_tuple(mg, leg, logx);
}

int dimuon_fotoproduction_analysis(const string& filename, string outname_pdf, string outname_png) {

    gStyle->SetOptStat(0);
    podio::ROOTReader reader;
    reader.openFile(filename);
    unsigned nEvents = reader.getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    TH2D* hEtaPt = new TH2D("hEtaPt", "Muon #eta vs p_{T};#eta;p_{T} [GeV]", 100, -6., 6., 100, 0., 7.);

    TH2D* hx_Q2 = new TH2D("hx_Q2", "Muon x vs Q^{2}; x; Q^{2}[GeV^{2}]", 100, 1e-6, 1e-3, 100, 1e-2, 1.0);
    TH2D* hEta1_Eta2 = new TH2D("hEta1_Eta2", "Muon #eta_{+} vs #eta_{-}; #eta_{+} (PDG=-13); #eta_{-} (PDG=13)", 100, -6., 6., 100, -6., 6.);

    constexpr int NBINS = 60; 

    TH1D* hQ2_all = new TH1D("hQ2_all", "Q^{2}", NBINS, 1e-2, 1.0);
    TH1D* hX_all  = new TH1D("hX_all", "x",     NBINS, 1e-6, 1e-3);
    TH1D* hY_all  = new TH1D("hY_all", "y",     NBINS, 0., 1.);
    TH1D* hW_all  = new TH1D("hW_all", "W",     NBINS, 0., 160.);

    TH1D* hQ2_in[3], *hX_in[3], *hY_in[3], *hW_in[3];
    for (int i = 0; i < 3; ++i) {
        hQ2_in[i] = new TH1D(Form("hQ2_in_%d", i), "", NBINS, 1e-2, 1.0);
        hX_in[i]  = new TH1D(Form("hX_in_%d",  i), "", NBINS, 1e-6, 1e-3);
        hY_in[i]  = new TH1D(Form("hY_in_%d",  i), "", NBINS, 0., 1.);
        hW_in[i]  = new TH1D(Form("hW_in_%d",  i), "", NBINS, 0., 160.);
    }

    TH1D* hQ2_rec_all = new TH1D("hQ2_rec_all", "Q^{2} (rec)", NBINS, 1e-2, 1.0);
    TH1D* hX_rec_all  = new TH1D("hX_rec_all", "x (rec)",      NBINS, 1e-6, 1e-3);
    TH1D* hY_rec_all  = new TH1D("hY_rec_all", "y (rec)",      NBINS, 0., 1.);
    TH1D* hW_rec_all  = new TH1D("hW_rec_all", "W (rec)",      NBINS, 0., 160.);

    TH1D* hQ2_rec_in[3], *hX_rec_in[3], *hY_rec_in[3], *hW_rec_in[3];
    for (int i = 0; i < 3; ++i) {
        hQ2_rec_in[i] = new TH1D(Form("hQ2_rec_in_%d", i), "", NBINS, 1e-2, 1.0);
        hX_rec_in[i]  = new TH1D(Form("hX_rec_in_%d",  i), "", NBINS, 1e-6, 1e-3);
        hY_rec_in[i]  = new TH1D(Form("hY_rec_in_%d",  i), "", NBINS, 0., 1.);
        hW_rec_in[i]  = new TH1D(Form("hW_rec_in_%d",  i), "", NBINS, 0., 160.);
    }

    constexpr int    Z_NBINS = 100;
    constexpr double Z_MIN_MM = -4500.0;
    constexpr double Z_MAX_MM = -3600.0;

    constexpr int    DXY_NBINS  = 120;
    constexpr double DXY_MIN_MM = -1200.0; 
    constexpr double DXY_MAX_MM =  1200.0;

    constexpr int    DR_NBINS  = 120;
    constexpr double DR_MIN_MM = 0.0;   
    constexpr double DR_MAX_MM = 1200.0;

    constexpr double P_MIN = 0.0;     // GeV
    constexpr double P_MAX = 25.0;    // GeV
    constexpr int    P_NB  = 50;

    constexpr double E_MIN = 0.0;
    constexpr double E_MAX = 10.0;

    TH1D* hZ_proj = new TH1D("hZ_proj", "z proj. ; z [mm]; N", NBINS, Z_MIN_MM, Z_MAX_MM);

    TH1D* hZ_hits = new TH1D("hZ_hits", "z rec. hits nHCal; z [mm]; N", NBINS, Z_MIN_MM, Z_MAX_MM);
    // To do
    TH3D* hDxDyZ_layer = new TH3D("hDxDyZ_layer","diff rec. and proj.; dx [mm]; dy [mm]; z [mm]", NBINS, DXY_MIN_MM, DXY_MAX_MM, NBINS, DXY_MIN_MM, DXY_MAX_MM, NBINS, Z_MIN_MM, Z_MAX_MM); 
    // To to
    TH2D* hDxDy_all = new TH2D("hDxDy_all",
                           "dx, dy ; dx [mm]; dy [mm]",
                           DXY_NBINS, DXY_MIN_MM, DXY_MAX_MM, DXY_NBINS, DXY_MIN_MM, DXY_MAX_MM);
    // To do
    TH2D* hDrZ_layer = new TH2D("hDrZ_layer","diff rec. vs proj.; dr [mm]; z [mm]", NBINS, DR_MIN_MM, DR_MAX_MM, NBINS, Z_MIN_MM, Z_MAX_MM); 
    // To do
    TH1D* hDr_all = new TH1D("hDr_all",
                         "dr = #sqrt{dx^{2}+dy^{2}} sum by layers; dr [mm]; N",
                         DR_NBINS, DR_MIN_MM, DR_MAX_MM);

    TH2D* hE_z = new TH2D("hE_z", "Energy of hits vs z; z [mm]; E [GeV]", NBINS, Z_MIN_MM, Z_MAX_MM, NBINS, E_MIN, E_MAX);

    TH2D* hEsum_z = new TH2D("hEsum_z", "Sum of energy in layer vs z_{layer}; z_{layer} [mm]; E_{sum} [GeV]",
                            NBINS, Z_MIN_MM, Z_MAX_MM, NBINS, E_MIN, E_MAX);

    TH1D* hP_all_mu = new TH1D("hP_all_mu", "All muons MC; p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX);

    constexpr double DR_CUTS_CM[3] = {5.0, 7.0, 10.0};
    TH1D* hP_pass_dr[3] = {
        new TH1D("hP_pass_dr5cm",  "Accepted (dr<5cm);  p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX),
        new TH1D("hP_pass_dr7cm",  "Accepted (dr<7cm);  p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX),
        new TH1D("hP_pass_dr10cm", "Accepted (dr<10cm); p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX)
    };

    // To do
    const double MIP_ENERGY_GEV = 0.002; 
    const double E_CUT_FACTORS[3] = {1.0, 1.3, 1.7}; 
    TH1D* hP_pass_Ecut[3] = {
        new TH1D("hP_pass_E<1.0MIP", "Accepted (E<1.0 MIP); p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX),
        new TH1D("hP_pass_E<1.3MIP", "Accepted (E<1.3 MIP); p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX),
        new TH1D("hP_pass_E<1.7MIP", "Accepted (E<1.7 MIP); p_{MC} [GeV]; N", P_NB, P_MIN, P_MAX)
    };

    // To do
    TH1D* hP_pass_combo[3][3]; // [idr][ie]
    for (int idr=0; idr<3; ++idr) {
        for (int ie=0; ie<3; ++ie) {
            hP_pass_combo[idr][ie] = new TH1D(
                Form("hP_pass_dr%.0fcm_E<%.1fMIP", DR_CUTS_CM[idr], E_CUT_FACTORS[ie]),
                Form("Accepted (dr<%.0f cm & E<%.1f MIP); p_{MC} [GeV]; N", DR_CUTS_CM[idr], E_CUT_FACTORS[ie]),
                P_NB, P_MIN, P_MAX
            );
        }
    }

    const double DR_CUTS_MM[3] = {50.0, 70.0, 100.0};  

    auto tf = [](bool v){ return v ? "true" : "false"; };

    for (unsigned ev = 0; ev < nEvents; ev++) {
        auto frameData = reader.readNextEntry(podio::Category::Event);
        if (!frameData) continue;
        podio::Frame frame(std::move(frameData));

        auto& kinCol = frame.get<edm4eic::InclusiveKinematicsCollection>("InclusiveKinematicsTruth");
        auto& mcCol  = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        auto& recParts = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
        auto& projSegs = frame.get<edm4eic::TrackSegmentCollection>("CalorimeterTrackProjections");
        auto& hcalRec  = frame.get<edm4eic::CalorimeterHitCollection>("HcalEndcapNRecHits");
        auto& assocCol = frame.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

        if (!kinCol.isValid() || !mcCol.isValid()) continue;
        if (kinCol.empty()) continue;

        const auto& kin = kinCol.at(0);
        double Q2 = kin.getQ2(), x = kin.getX(), y = kin.getY(), W = kin.getW();

        edm4hep::MCParticle m1, m2;
        for (const auto& p : mcCol) {
            if (abs(p.getPDG()) != 13 || p.getGeneratorStatus() == 0) continue;
            if (p.getPDG() == 13)  m1 = p;
            else if (p.getPDG() == -13) m2 = p;
        }

        TLorentzVector v1(m1.getMomentum().x, m1.getMomentum().y, m1.getMomentum().z, m1.getEnergy());
        TLorentzVector v2(m2.getMomentum().x, m2.getMomentum().y, m2.getMomentum().z, m2.getEnergy());

        hEtaPt->Fill(v1.Eta(), v1.Pt());
        hEtaPt->Fill(v2.Eta(), v2.Pt());

        hx_Q2->Fill(x,Q2); 
        hEta1_Eta2->Fill(v2.Eta(), v1.Eta());

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

        if(nInNH >= 1) {
            hP_all_mu->Fill(v1.P());
            hP_all_mu->Fill(v2.P());
        }

        constexpr double THRESH_MM = 100.0; // 10 cm w mm

        bool m1_has_rec = false;
        bool m2_has_rec = false;


        struct TaggedReco { edm4eic::ReconstructedParticle reco; int muTag; /*1=m1, 2=m2*/ };
        vector<TaggedReco> matchedRecos;

        auto find_associated_reco = [&](const edm4hep::MCParticle& mc, int muTag)->bool {
            if (!mc.isAvailable()) return false;
            try {
                if (!assocCol.isValid() || assocCol.empty()) {
                    return false;
                }
                auto simIDs = assocCol.simID();
                auto recIDs = assocCol.recID();
                const uint32_t mc_idx = mc.getObjectID().index;
                bool found = false;
                for (size_t i=0; i<assocCol.size() && i<simIDs.size() && i<recIDs.size(); ++i) {
                    if (simIDs[i] == mc_idx) {
                        uint32_t ridx = recIDs[i];
                        if (!recParts.isValid() || ridx >= recParts.size()) {continue;}
                        auto reco = recParts.at(ridx);
                        if (reco.isAvailable()) {
                            matchedRecos.push_back({reco, muTag});
                            found = true;
                        }
                    }
                }
                return found;
            } catch (...) {
                return false;
            }
        };

        int assocCount = 0;
        if (m1.isAvailable() && abs(m1.getPDG())==13) assocCount += find_associated_reco(m1, 1) ? 1 : 0;
        if (m2.isAvailable() && abs(m2.getPDG())==13) assocCount += find_associated_reco(m2, 2) ? 1 : 0;

        if (!recParts.isValid() || !projSegs.isValid() || !hcalRec.isValid() || recParts.empty() || projSegs.empty() || hcalRec.empty()) continue;
        
        map<double, double> layerData;
        for (const auto& hit : hcalRec) {    
            double z = hit.getPosition().z;
            double zBin = round(z);        
            layerData[zBin] += hit.getEnergy(); 
            hZ_hits->Fill(z);                
            hE_z->Fill(z, hit.getEnergy());   
        }

        for (const auto& [zValue, sumEnergy] : layerData) {hEsum_z->Fill(zValue, sumEnergy);}

        struct TaggedTrack { edm4eic::Track tr; int muTag; };
        vector<TaggedTrack> allTracks;
        for (const auto& R : matchedRecos) {
            for (const auto& tr : R.reco.getTracks()) {
                if (tr.isAvailable()) allTracks.push_back({tr, R.muTag});
            }
        }
        struct SegTag { edm4eic::TrackSegment seg; int muTag; };
        vector<SegTag> segsTagged;
        for (const auto& seg : projSegs) {
            auto linkedTr = seg.getTrack();
            if (!linkedTr.isAvailable()) continue;
            for (const auto& TT : allTracks) {
                if (linkedTr.getObjectID() == TT.tr.getObjectID()) {
                    segsTagged.push_back({seg, TT.muTag});
                    break;
                }
            }
        }

        for (const auto& ST : segsTagged) {
            auto points = ST.seg.getPoints();
            double segBestMin = 1e12;
            double segHitEnergy = 1e12;

            for (const auto& pt : points) {
                const auto& ptPosition = pt.position;  

                if(pt.system != 113) continue; 
                //if(pt.surface != 1) continue;

                double localHitEnergy = 1e12;
                double localMin = 1e12;
                for (const auto& hit : hcalRec) {
                    const auto& hpos = hit.getPosition();                
                    const double d = dist3(ptPosition, hpos);
                    if (d < localMin) { localMin = d;  }
                }
                if (localMin < segBestMin) {
                    segBestMin = localMin;
                }

                hZ_proj->Fill(ptPosition.z);
            } 

            if (segBestMin <= THRESH_MM) {
                cout << "[MATCH] muTag=" << ST.muTag
                    << " d=" << segBestMin << " mm (<= " << THRESH_MM << ")\n";
                if (ST.muTag == 1) {
                    m1_has_rec = true; 
                    for (int idr=0; idr<3; ++idr) if (segBestMin < DR_CUTS_MM[idr]) hP_pass_dr[idr]->Fill(v1.P());
                }
                if (ST.muTag == 2) {
                    m2_has_rec = true; 
                    for (int idr=0; idr<3; ++idr) if (segBestMin < DR_CUTS_MM[idr]) hP_pass_dr[idr]->Fill(v2.P());
                }
            } else {
                cout << "[INFO] no match <= " << THRESH_MM
                    << " mm; best distance (this seg) = " << segBestMin
                    << " mm, for muTag=" << ST.muTag << "\n";
            }
        } 
        

        int nInNH_rec_local = 0;
        if (m1_has_rec) ++nInNH_rec_local;
        if (m2_has_rec) ++nInNH_rec_local;

        cout << "[RES] m1_has_rec=" << (m1_has_rec ? "true" : "false")
                << ", m2_has_rec=" << (m2_has_rec ? "true" : "false")
                << " | nInNH_rec_local=" << nInNH_rec_local << "\n";

        hQ2_rec_all->Fill(Q2);
        hX_rec_all->Fill(x);
        hY_rec_all->Fill(y);
        hW_rec_all->Fill(W);

        if (nInNH_rec_local >= 0 && nInNH_rec_local <= 2) {
            hQ2_rec_in[nInNH_rec_local]->Fill(Q2);
            hX_rec_in[nInNH_rec_local]->Fill(x);
            hY_rec_in[nInNH_rec_local]->Fill(y);
            hW_rec_in[nInNH_rec_local]->Fill(W);
        }
    }

    TCanvas* canvas_sim = new TCanvas("canvas_sim", "muon analysis", 1600, 800);
    canvas_sim->Divide(4, 2);
    canvas_sim->cd(1); hEtaPt->Draw("COLZ");

    auto [mg_q2, leg_q2, log_q2] = makeEffMultiGraph(hQ2_all, hQ2_in, "Geom. acc. vs Q^{2}", "Q^{2} [GeV^{2}]", true);
    canvas_sim->cd(2); if (log_q2) gPad->SetLogx(); gPad->SetGrid(); mg_q2->Draw("A"); leg_q2->Draw();

    auto [mg_x, leg_x, log_x] = makeEffMultiGraph(hX_all, hX_in, "Geom. acc. vs x", "x", true);
    canvas_sim->cd(3); if (log_x) gPad->SetLogx(); gPad->SetGrid(); mg_x->Draw("A"); leg_x->Draw();

    auto [mg_y, leg_y, log_y] = makeEffMultiGraph(hY_all, hY_in, "Geom. acc. vs y", "y", false);
    canvas_sim->cd(4); gPad->SetGrid(); mg_y->Draw("A"); leg_y->Draw();

    auto [mg_w, leg_w, log_w] = makeEffMultiGraph(hW_all, hW_in, "Geom. acc. vs W", "W [GeV]", false);
    canvas_sim->cd(5); gPad->SetGrid(); mg_w->Draw("A"); leg_w->Draw();

    canvas_sim->cd(6); gPad->SetLogx(); gPad->SetGrid(); gPad->SetLogy(); hx_Q2->Draw("COLZ");
    canvas_sim->cd(7); hEta1_Eta2->Draw("COLZ"); 
    
    

    auto* lv1 = makeLine(ETA_MIN, -6, ETA_MIN, 6);
    auto* lv2 = makeLine(ETA_MAX, -6, ETA_MAX, 6);
    auto* lh1 = makeLine(hEta1_Eta2->GetXaxis()->GetXmin(), ETA_MIN,
                        hEta1_Eta2->GetXaxis()->GetXmax(), ETA_MIN);
    auto* lh2 = makeLine(hEta1_Eta2->GetXaxis()->GetXmin(), ETA_MAX,
                        hEta1_Eta2->GetXaxis()->GetXmax(), ETA_MAX);

    canvas_sim->SaveAs(outname_pdf.c_str());
    canvas_sim->SaveAs(outname_png.c_str());

    TCanvas* canvas_rec = new TCanvas("canvas_rec", "muon analysis rec", 1600, 800);
    canvas_rec->Divide(4, 2);
    canvas_rec->cd(1); hEtaPt->Draw("COLZ");

    auto [mg_rec_q2, leg_rec_q2, log_rec_q2] = makeEffMultiGraph(hQ2_rec_all, hQ2_rec_in, "Geom. acc. vs Q^{2}", "Q^{2} [GeV^{2}]", true);
    canvas_rec->cd(2); if (log_rec_q2) gPad->SetLogx(); gPad->SetGrid(); mg_rec_q2->Draw("A"); leg_rec_q2->Draw();

    auto [mg_rec_x, leg_rec_x, log_rec_x] = makeEffMultiGraph(hX_rec_all, hX_rec_in, "Geom. acc. vs x", "x", true);
    canvas_rec->cd(3); if (log_rec_x) gPad->SetLogx(); gPad->SetGrid(); mg_rec_x->Draw("A"); leg_rec_x->Draw();

    auto [mg_rec_y, leg_rec_y, log_rec_y] = makeEffMultiGraph(hY_rec_all, hY_rec_in, "Geom. acc. vs y", "y", false);
    canvas_rec->cd(4); gPad->SetGrid(); mg_rec_y->Draw("A"); leg_rec_y->Draw();

    auto [mg_rec_w, leg_rec_w, log_rec_w] = makeEffMultiGraph(hW_rec_all, hW_rec_in, "Geom. acc. vs W", "W [GeV]", false);
    canvas_rec->cd(5); gPad->SetGrid(); mg_rec_w->Draw("A"); leg_rec_w->Draw();

    canvas_rec->cd(6); gPad->SetLogx(); gPad->SetGrid(); hx_Q2->Draw("COLZ");
    canvas_rec->cd(7); hEta1_Eta2->Draw("COLZ");

    string rec_outname_pdf = addPrefixAfterSlash(outname_pdf, "rec_");
    string rec_outname_png = addPrefixAfterSlash(outname_png, "rec_");
    canvas_rec->SaveAs(rec_outname_png.c_str());
    canvas_rec->SaveAs(rec_outname_pdf.c_str());

    TCanvas* canvas = new TCanvas("canvas", "Extra histograms (ALL-IN-ONE)", 1800, 1200); 
    canvas->Divide(3,3);                                                                 
    canvas->cd(1); gPad->SetGrid(); hZ_proj->Draw();                                     
    canvas->cd(2); gPad->SetGrid(); hZ_hits->Draw();                                     
    canvas->cd(3); gPad->SetGrid(); hDxDy_all->Draw("COLZ");                              
    canvas->cd(4); gPad->SetGrid(); hDr_all->Draw();                                      
    canvas->cd(5); gPad->SetGrid(); hE_z->Draw("COLZ");                                  
    canvas->cd(6); gPad->SetGrid(); hEsum_z->Draw("COLZ");                                
    canvas->cd(7); gPad->SetGrid(); hP_all_mu->SetLineColor(kBlack); hP_all_mu->Draw();  
    canvas->cd(8); gPad->SetGrid();                                                       
    TH1D* hEff_dr[3]; 
    for (int idr=0; idr<3; ++idr) { 
        hEff_dr[idr] = (TH1D*)hP_pass_dr[idr]->Clone( 
            Form("hEff_dr%.0fcm", DR_CUTS_CM[idr]) 
        ); 
        hEff_dr[idr]->SetTitle("Matching efficiency (progi dr); p_{MC} [GeV]; eff"); 
        hEff_dr[idr]->Divide(hP_all_mu); 
    } 

    hEff_dr[0]->SetLineColor(kBlue);     hEff_dr[0]->SetMinimum(0.0); hEff_dr[0]->SetMaximum(1.1); 
    hEff_dr[1]->SetLineColor(kRed);   
    hEff_dr[2]->SetLineColor(kGreen+2); 

    hEff_dr[0]->Draw("HIST");            
    hEff_dr[1]->Draw("HIST SAME");       
    hEff_dr[2]->Draw("HIST SAME");       
    { auto leg8 = new TLegend(0.55,0.65,0.88,0.88);
    leg8->AddEntry(hEff_dr[0],"dr < 5 cm","l");   
    leg8->AddEntry(hEff_dr[1],"dr < 7 cm","l");   
    leg8->AddEntry(hEff_dr[2],"dr < 10 cm","l");  
    leg8->Draw(); }  
                                                                   
    canvas->cd(9); gPad->SetGrid();                                                      
    hP_pass_Ecut[0]->SetLineColor(kBlue);   hP_pass_Ecut[0]->Draw();                   
    hP_pass_Ecut[1]->SetLineColor(kRed);    hP_pass_Ecut[1]->Draw("SAME");             
    hP_pass_Ecut[2]->SetLineColor(kGreen+2);hP_pass_Ecut[2]->Draw("SAME");              
    { auto leg9 = new TLegend(0.55,0.65,0.88,0.88);                                     
      leg9->AddEntry(hP_pass_Ecut[0],"E < 1.0 MIP","l");                               
      leg9->AddEntry(hP_pass_Ecut[1],"E < 1.3 MIP","l");                              
      leg9->AddEntry(hP_pass_Ecut[2],"E < 1.7 MIP","l");                               
      leg9->Draw(); }                                                                    
    canvas->SaveAs(addPrefixAfterSlash(outname_png, "extra_ALL_").c_str());              
    canvas->SaveAs(addPrefixAfterSlash(outname_pdf, "extra_ALL_").c_str());               

    return 0;
}