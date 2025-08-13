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

using namespace std;
using namespace ROOT;
using namespace TMath;
using namespace edm4hep;

constexpr double ETA_MIN = -4.14, ETA_MAX = -1.16;

inline bool inNHCal(double eta) {return (eta >= ETA_MIN && eta <= ETA_MAX);}

int dimuon_fotoproduction_analysis(const string& filename, string outname_pdf, string outname_png) {

    int debug = 0;

    gStyle->SetOptStat(0);
    podio::ROOTReader reader;
    reader.openFile(filename);
    unsigned nEvents = reader.getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    TH2D* hEtaPt = new TH2D("hEtaPt", "Muon #eta vs p_{T};#eta;p_{T} [GeV]", 120, -6., 6., 100, 0., 5.);

    const int nb = 60;
    TH1D* hQ2_all = new TH1D("hQ2_all", "Q^{2}", nb, 1e-2, 1.0);
    TH1D* hX_all  = new TH1D("hX_all", "x",     nb, 1e-6, 1e-3);
    TH1D* hY_all  = new TH1D("hY_all", "y",     nb, 0., 1.);
    TH1D* hW_all  = new TH1D("hW_all", "W",     nb, 0., 160.);

    TH1D *hQ2_in[3], *hX_in[3], *hY_in[3], *hW_in[3];
    for (int i = 0; i < 3; ++i) {
        hQ2_in[i] = new TH1D(Form("hQ2_in_%d", i), "", nb, 1e-2, 1.0);
        hX_in[i]  = new TH1D(Form("hX_in_%d",  i), "", nb, 1e-6, 1e-3);
        hY_in[i]  = new TH1D(Form("hY_in_%d",  i), "", nb, 0., 1.);
        hW_in[i]  = new TH1D(Form("hW_in_%d",  i), "", nb, 0., 160.);
    }

    TH1D* hQ2_rec_all = new TH1D("hQ2_rec_all", "Q^{2} (rec)", nb, 1e-2, 1.0);
    TH1D* hX_rec_all  = new TH1D("hX_rec_all", "x (rec)",      nb, 1e-6, 1e-3);
    TH1D* hY_rec_all  = new TH1D("hY_rec_all", "y (rec)",      nb, 0., 1.);
    TH1D* hW_rec_all  = new TH1D("hW_rec_all", "W (rec)",      nb, 0., 160.);

    TH1D *hQ2_rec_in[3], *hX_rec_in[3], *hY_rec_in[3], *hW_rec_in[3];
    for (int i = 0; i < 3; ++i) {
        hQ2_rec_in[i] = new TH1D(Form("hQ2_rec_in_%d", i), "", nb, 1e-2, 1.0);
        hX_rec_in[i]  = new TH1D(Form("hX_rec_in_%d",  i), "", nb, 1e-6, 1e-3);
        hY_rec_in[i]  = new TH1D(Form("hY_rec_in_%d",  i), "", nb, 0., 1.);
        hW_rec_in[i]  = new TH1D(Form("hW_rec_in_%d",  i), "", nb, 0., 160.);
    }

    for (unsigned ev = 0; ev < nEvents; ev++) {
        auto frameData = reader.readNextEntry(podio::Category::Event);
        if (!frameData) continue;
        podio::Frame frame(std::move(frameData));
        
        if (debug && ev == 0) {
            std::cout << "==== Available collections in event 0 ====\n";
            for (const auto& name : frame.getAvailableCollections()) {
                std::cout << "  - " << name << "\n";
            }
            std::cout << "==========================================\n";
        }

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
            else if (p.getPDG() == -13) m2 = p;
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

        //================================ New module ================================//

        auto pickKinematicsName = [&](const podio::Frame& f)->std::string {
            const char* cands[] = {
                "InclusiveKinematicsElectron",
                "InclusiveKinematicsDA",
                "InclusiveKinematicsJB",
                "InclusiveKinematicsML",
                "InclusiveKinematicsSigma",
                "InclusiveKinematicsESigma",
                "InclusiveKinematicseSigma"
            };
            for (auto nm : cands) {
                auto& c = f.get<edm4eic::InclusiveKinematicsCollection>(nm);
                if (c.isValid() && !c.empty()) return nm;
            }
            return "";
        };

        std::string kinName = pickKinematicsName(frame);
        bool haveKinRec = !kinName.empty();
        double Q2r=0, xr=0, yr=0, Wr=0;

        if (debug) std::cout << "  [DEBUG] REC kinematics collection: " << (haveKinRec ? kinName : "<none>") << "\n";
        if (haveKinRec) {
            auto& kinRecCol = frame.get<edm4eic::InclusiveKinematicsCollection>(kinName);
            const auto& kinRec = kinRecCol.at(0);
            Q2r = kinRec.getQ2(); xr = kinRec.getX(); yr = kinRec.getY(); Wr = kinRec.getW();
            if (debug) std::cout << "  RECO Q2=" << Q2r << " x=" << xr << " y=" << yr << " W=" << Wr << "\n";
        }

        auto& segCol = frame.get<edm4eic::TrackSegmentCollection>("CentralTrackSegments");
        auto& rpCol  = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
        auto& rcCol  = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
        auto& rcPID  = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedRealPIDParticles");

        constexpr const char* HCAL_NAME = "HcalEndcapNRecHits";
        auto& hcalN = frame.get<edm4eic::CalorimeterHitCollection>(HCAL_NAME);

        if (debug) {
        std::cout << " RECO collections:\n";
        std::cout << " ReconstructedParticles valid=" << rpCol.isValid()  << " size=" << rpCol.size()  << "\n";
        std::cout << " ReconstructedChargedParticles valid=" << rcCol.isValid() << " size=" << rcCol.size() << "\n";
        std::cout << " ReconstructedChargedRealPIDParticles valid=" << rcPID.isValid() << " size=" << rcPID.size() << "\n";
        std::cout << " [DEBUG] HCAL used: " << HCAL_NAME
                    << " valid=" << hcalN.isValid() << " size=" << hcalN.size() << "\n";
        }

        bool haveSeg   = segCol.isValid() && !segCol.empty();
        bool haveHits  = hcalN.isValid() && !hcalN.empty();     
        bool haveAnyRP = (rpCol.isValid()  && !rpCol.empty())
                    || (rcCol.isValid()  && !rcCol.empty())
                    || (rcPID.isValid()  && !rcPID.empty());

        if (!haveSeg || !haveHits || !haveAnyRP) {
        if (debug) {
            if (!haveSeg)   std::cout << " [WARN] brak CentralTrackSegments\n";
            if (!haveHits)  std::cout << " [WARN] brak HCAL hitów w " << HCAL_NAME << "\n";
            if (!haveAnyRP) std::cout << " [WARN] brak RP kolekcji do znaku ładunku\n";
        }
        } else {

        const float zBin = 10.f; // mm
        std::unordered_map<int,int> zcount;
        zcount.reserve(hcalN.size());
        for (const auto& h : hcalN) {
            int bz = int(std::floor(h.getPosition().z / zBin));
            ++zcount[bz];
        }
        int bestBz = 0, bestCnt = -1;
        for (auto& kv : zcount) if (kv.second > bestCnt) { bestCnt = kv.second; bestBz = kv.first; }
        const float zTarget = bestBz * zBin;

        const float zWin   = 80.f; // mm
        const float voxel  = 50.f; // mm

        struct K2D { int ix, iy; };
        struct H2DHash { size_t operator()(const K2D& k) const noexcept {
            return (k.ix*73856093) ^ (k.iy*19349663);
        }};
        struct H2DEq { bool operator()(const K2D& a, const K2D& b) const noexcept {
            return a.ix==b.ix && a.iy==b.iy;
        }};
        auto k2 = [&](float x, float y)->K2D {
            return { int(std::floor(x/voxel)), int(std::floor(y/voxel)) };
        };

        std::unordered_map<K2D, std::vector<edm4eic::CalorimeterHit>, H2DHash, H2DEq> grid;
        grid.reserve(hcalN.size());
        for (const auto& h : hcalN) {
            const auto& p = h.getPosition();
            if (std::fabs(p.z - zTarget) > zWin) continue;
            grid[k2(p.x, p.y)].push_back(h);
        }

        if (debug) {
            size_t kept = 0;
            for (auto& kv : grid) kept += kv.second.size();
            std::cout << " [DEBUG] zTarget=" << zTarget << " mm, kept NHCal hits ~layer: "
                    << kept << " cells=" << grid.size() << "\n";
        }
        auto oid = [](const auto& obj){ return obj.getObjectID(); };
        struct OIDHash {
            size_t operator()(const podio::ObjectID& id) const noexcept {
                return (static_cast<size_t>(static_cast<uint32_t>(id.collectionID)) << 32)
                    ^ static_cast<uint32_t>(id.index);
            }
        };
        struct OIDEq {
            bool operator()(const podio::ObjectID& a, const podio::ObjectID& b) const noexcept {
                return a.collectionID == b.collectionID && a.index == b.index;
            }
        };

        std::unordered_set<podio::ObjectID, OIDHash, OIDEq> muonTracks;
        std::unordered_map<podio::ObjectID, int, OIDHash, OIDEq> muonTrackSign; 

        auto isMuonRP = [&](const edm4eic::ReconstructedParticle& rp)->bool {
            if constexpr (requires { rp.getPDG(); }) {
                int pdg = rp.getPDG();
                return std::abs(pdg) == 13; 
            } else {
                return false;
            }
        };

        auto add_muons_from = [&](const edm4eic::ReconstructedParticleCollection& col){
            if (!col.isValid() || col.empty()) return;
            for (const auto& rp : col) {
                if (!isMuonRP(rp)) continue;
                double q = rp.getCharge();
                int sgn = (q > 0) ? +1 : (q < 0 ? -1 : 0);
                if (sgn == 0) continue; 
                for (const auto& trh : rp.getTracks()) {
                    if (!trh.isAvailable()) continue;
                    auto id = oid(trh);
                    muonTracks.insert(id);
                    muonTrackSign[id] = sgn;
                }
            }
        };

        add_muons_from(rcPID);
        if (muonTracks.empty()) add_muons_from(rcCol);
        if (muonTracks.empty()) add_muons_from(rpCol);

        if (debug) {
            std::cout << "  [DEBUG] muonTracks=" << muonTracks.size()
                    << " muonTrackSign=" << muonTrackSign.size() << "\n";
        }
            auto projectToZ = [&](const edm4eic::TrackSegment& ts, float zt, edm4hep::Vector3f& out)->bool {
                auto pts = ts.getPoints();
                if (pts.size() < 2) return false;
                const edm4eic::TrackPoint *pA=nullptr, *pB=nullptr;
                for (const auto& tp : pts) {
                    if (!pA || tp.pathlength > pA->pathlength) { pB = pA; pA = &tp; }
                    else if (!pB || tp.pathlength > pB->pathlength) { pB = &tp; }
                }
                if (!pA || !pB) return false;
                auto A = pA->position; auto B = pB->position;
                float dz = (A.z - B.z);
                if (std::fabs(dz) < 1e-3f) return false;
                float t = (zt - A.z) / ((A.z - B.z));
                out.x = A.x + t * (A.x - B.x);
                out.y = A.y + t * (A.y - B.y);
                out.z = zt;
                return true;
            };

            auto anyHitNear = [&](const edm4hep::Vector3f& p, float R)->bool {
                const float R2 = R*R;
                auto key = k2(p.x, p.y);
                for (int dx=-1; dx<=1; ++dx)
                for (int dy=-1; dy<=1; ++dy) {
                    auto it = grid.find( K2D{key.ix+dx, key.iy+dy} );
                    if (it == grid.end()) continue;
                    for (const auto& h : it->second) {
                        const auto& hp = h.getPosition();
                        float dxp = p.x - hp.x, dyp = p.y - hp.y;
                        if (dxp*dxp + dyp*dyp <= R2) return true;
                    }
                }
                return false;
            };

            const float Rmatch = 50.f; // 5 cm

            bool matchedPlus=false, matchedMinus=false;
            int segTried=0, segProj=0, segHit=0;

            for (size_t i=0; i<segCol.size(); ++i) {
                const auto& ts = segCol.at(i);
                auto trh = ts.getTrack();
                if (!trh.isAvailable()) continue;
                ++segTried;

                edm4hep::Vector3f pproj{};
                if (!projectToZ(ts, zTarget, pproj)) continue;
                ++segProj;

                if (anyHitNear(pproj, Rmatch)) {
                    ++segHit;

                    auto id = oid(trh);
                    if (muonTracks.find(id) == muonTracks.end()) {
                        continue;
                    }

                    int sgn = muonTrackSign.count(id) ? muonTrackSign[id] : 0;
                    if (sgn == +1)      matchedPlus  = true;
                    else if (sgn == -1) matchedMinus = true;

                    if (matchedPlus && matchedMinus) break;
                }
            }

            if (debug) {
                std::cout << "  [DEBUG] projZ=" << zTarget << " mm | segments tried=" << segTried
                        << " projected=" << segProj << " hit-near=" << segHit
                        << " → plus=" << matchedPlus << " minus=" << matchedMinus << "\n";
            }

            int matchedMuons = (matchedMinus ? 1 : 0) + (matchedPlus ? 1 : 0);
            if (debug) std::cout << "  [DEBUG] matchedMuons=" << matchedMuons << "\n";

            if (haveKinRec) {
                hQ2_rec_all->Fill(Q2r);
                hX_rec_all->Fill(xr);
                hY_rec_all->Fill(yr);
                hW_rec_all->Fill(Wr);

                if (matchedMuons >= 0 && matchedMuons <= 2) {
                    hQ2_rec_in[matchedMuons]->Fill(Q2r);
                    hX_rec_in[matchedMuons]->Fill(xr);
                    hY_rec_in[matchedMuons]->Fill(yr);
                    hW_rec_in[matchedMuons]->Fill(Wr);
                }
            } else if (debug) {
                std::cout << "  [DEBUG] brak REC kinematyki — nie wypełniam h*_rec_* w tym evencie\n";
            }
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

        int added = 0;

        for (int i = 0; i < 3; ++i) {
            std::vector<double> x_vals, y_vals;
            for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
                double all = h_all->GetBinContent(b);
                double sel = h_in[i]->GetBinContent(b);
                if (all < 5) continue;
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

        {
            std::vector<double> x_vals, y_vals;
            for (int b = 1; b <= h_all->GetNbinsX(); ++b) {
                double all = h_all->GetBinContent(b);
                double sel = h_in[1]->GetBinContent(b) + h_in[2]->GetBinContent(b);
                if (all < 5) continue;
                x_vals.push_back(h_all->GetBinCenter(b));
                y_vals.push_back(100. * sel / all);
            }
            if (!x_vals.empty()) {
                TGraph* gsum = new TGraph((int)x_vals.size(), x_vals.data(), y_vals.data());
                gsum->SetMarkerStyle(25);
                gsum->SetMarkerColor(kMagenta + 1);
                gsum->SetFillColor(kMagenta+1);
                gsum->SetLineColor(0);
                gsum->SetMarkerSize(1.0);
                mg->Add(gsum, "P");
                leg->AddEntry(gsum, "All eligible", "p");
                ++added;
            }
        }

        mg->SetTitle(Form("%s;%s;geom. acc. [%%]", title, xlabel));

        if (added == 0) {
            TH1D* hframe = new TH1D("hframe_tmp", "", 10, h_all->GetXaxis()->GetXmin(), h_all->GetXaxis()->GetXmax());
            hframe->SetMinimum(0.0);
            hframe->SetMaximum(110.0);
            hframe->GetXaxis()->SetTitle(xlabel);
            hframe->GetYaxis()->SetTitle("geom. acc. [%]");
            hframe->Draw();
        }

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

    TCanvas* canvas_rec = new TCanvas("canvas_rec", "muon analysis rec", 1600, 800);
    canvas_rec->Divide(3, 2);
    //canvas_rec->cd(1); hEtaPt->Draw("COLZ");

    auto [mg_rec_q2, leg_rec_q2, log_rec_q2] = makeEffMultiGraph(hQ2_rec_all, hQ2_rec_in, "Geom. acc. vs Q^{2}", "Q^{2} [GeV^{2}]", true);
    canvas_rec->cd(2); if (log_rec_q2) gPad->SetLogx(); gPad->SetGrid(); mg_rec_q2->Draw("A"); leg_rec_q2->Draw();

    auto [mg_rec_x, leg_rec_x, log_rec_x] = makeEffMultiGraph(hX_rec_all, hX_rec_in, "Geom. acc. vs x", "x", true);
    canvas_rec->cd(3); if (log_rec_x) gPad->SetLogx(); gPad->SetGrid(); mg_rec_x->Draw("A"); leg_rec_x->Draw();

    auto [mg_rec_y, leg_rec_y, log_rec_y] = makeEffMultiGraph(hY_rec_all, hY_rec_in, "Geom. acc. vs y", "y", false);
    canvas_rec->cd(4); gPad->SetGrid(); mg_rec_y->Draw("A"); leg_rec_y->Draw();

    auto [mg_rec_w, leg_rec_w, log_rec_w] = makeEffMultiGraph(hW_rec_all, hW_rec_in, "Geom. acc. vs W", "W [GeV]", false);
    canvas_rec->cd(5); gPad->SetGrid(); mg_rec_w->Draw("A"); leg_rec_w->Draw();

    canvas_rec->SaveAs("results/nhcal_dimuon_fotoproduction/rec_analysis_epic_full.png");
    canvas_rec->SaveAs("results/nhcal_dimuon_fotoproduction/rec_analysis_epic_full.pdf");

    return 0;
}