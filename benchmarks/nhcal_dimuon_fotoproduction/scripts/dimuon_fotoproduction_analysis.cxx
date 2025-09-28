#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
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
#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Shapes.h"


#include "podio/Frame.h"
#include "podio/CollectionBase.h"
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

#include "edm4eic/vector_utils_legacy.h"
#include "edm4hep/Vector3f.h"

#include "edm4eic/Track.h"
#include "edm4eic/TrackSegment.h"
#include "edm4eic/TrackSegmentCollectionData.h"
#include "edm4eic/TrackPoint.h"
#include "edm4eic/TrackParameters.h"

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

dd4hep::Detector* det = NULL;

constexpr double ETA_MIN = -4.14, ETA_MAX = -1.16;

inline bool inNHCal(double eta) {return (eta >= ETA_MIN && eta <= ETA_MAX);}

inline string addPrefixAfterSlash(const string& path,
                                       const string& prefix) {
    const auto slash = path.find_last_of('/');
    if (slash == string::npos)
        return prefix + path;
    return path.substr(0, slash + 1) + prefix + path.substr(slash + 1);
}

TLine* makeLine(double x1, double y1, double x2, double y2) {
    TLine* l = new TLine(x1, y1, x2, y2);
    l->SetLineColor(kRed);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw("same");
    return l;
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

inline double getPlasticThicknessMM(dd4hep::Detector& det,
                                    const dd4hep::DDSegmentation::BitFieldCoder* dec,
                                    dd4hep::DDSegmentation::CellID cid,
                                    int slice_idx,
                                    int plastic_slice_value = 3)
{
    const double NaN = numeric_limits<double>::quiet_NaN();
    try {
        if (!dec) throw runtime_error("decoder==nullptr");
        if (dec->get(cid, slice_idx) != plastic_slice_value)
            throw runtime_error("cell is not plastic");

        auto de = det.volumeManager().lookupDetElement(cid);
        if (!de.isValid()) throw runtime_error("lookupDetElement failed");

        auto pv = de.placement();
        if (!pv.isValid()) throw runtime_error("DetElement has no placement");

        auto vol = pv.volume();
        if (!vol.isValid()) throw runtime_error("Invalid Volume");

        auto* box = dynamic_cast<TGeoBBox*>(vol.solid().ptr());
        if (!box) throw runtime_error("Solid is not TGeoBBox");

        const double dz_cm = box->GetDZ();                     
        const double thickness_cm = 2.0 * dz_cm;
        return thickness_cm;
    } catch (const exception& e) {
        cerr << "[WARN] getPlasticThicknessMM: " << e.what() << " (cellID=" << cid << ")\n";
        return NaN;
    }
}

struct GeomState {
    double x, y, z;    
    double px, py, pz;  
};

// --- tiny helpers ---
inline double hypot2(double a, double b){ return sqrt(a*a + b*b); }
inline void unitXY(double px, double py, double& tx, double& ty){
    double pt = max(hypot2(px,py), 1e-16);
    tx = px/pt; ty = py/pt;
}
inline void rot90(double vx, double vy, double& rx, double& ry){
    rx = -vy; ry =  vx; // +90 degree
}

// geometry-only helix (or line) projection from exactly two states
inline pair<double,double>
trackXYatZ_fromTwoStates(double x1,double y1,double z1,
                         double px1,double py1,double pz1,
                         double x2,double y2,double z2,
                         double px2,double py2,double pz2,
                         double zTarget)
{
    constexpr double EPS    = 1e-12;
    constexpr double EPSANG = 1e-6;   // radians

    // unit tangents in XY
    double t1x,t1y,t2x,t2y; unitXY(px1,py1,t1x,t1y); unitXY(px2,py2,t2x,t2y);

    // average dip angle tan(lambda)
    const double pt1 = max(hypot2(px1,py1), 1e-16);
    const double pt2 = max(hypot2(px2,py2), 1e-16);
    const double tanL1 = pz1/pt1, tanL2 = pz2/pt2, tanL = 0.5*(tanL1+tanL2);

    if (abs(tanL) < EPS) {
        // nearly horizontal track -> z-param is unstable; return nearest state's XY
        return (abs(zTarget - z1) <= abs(zTarget - z2)) ? pair{x1,y1} : pair{x2,y2};
    }

    // signed angle between tangents (curvature sign & fallback)
    const double cross = t1x*t2y - t1y*t2x;
    const double dot   = clamp(t1x*t2x + t1y*t2y, -1.0, 1.0);
    const double dPhi  = atan2(cross, dot);      // (-pi, pi]
    const int    sgn   = (dPhi >= 0) ? +1 : -1;       // CCW/CW by geometry

    // normals (left normals wrt tangent)
    double n1x,n1y,n2x,n2y; rot90(t1x,t1y,n1x,n1y); rot90(t2x,t2y,n2x,n2y);

    // radius from shifted-normals intersection (least-squares)
    const double dx = x2 - x1, dy = y2 - y1;
    const double dn_x = n1x - n2x, dn_y = n1y - n2y;
    const double denom = dn_x*dn_x + dn_y*dn_y;

    const bool near_parallel_normals = (denom < 1e-10);
    if (near_parallel_normals || abs(dPhi) < EPSANG) {
        // almost straight -> linear extrapolation in XY by z
        const bool anchorAt1 = (abs(zTarget - z1) <= abs(zTarget - z2));
        double xa = anchorAt1 ? x1 : x2;
        double ya = anchorAt1 ? y1 : y2;
        double za = anchorAt1 ? z1 : z2;
        double tax= anchorAt1 ? t1x: t2x;
        double tay= anchorAt1 ? t1y: t2y;
        double s = (zTarget - za) / tanL;
        return { xa + s*tax, ya + s*tay };
    }

    double R = (dx*dn_x + dy*dn_y) / denom;
    if (R < 0) R = -R;  // keep radius positive; orientation via 'sgn'

    // circle center (average of two per-state centers)
    const double C1x = x1 + sgn * R * n1x;
    const double C1y = y1 + sgn * R * n1y;
    const double C2x = x2 + sgn * R * n2x;
    const double C2y = y2 + sgn * R * n2y;
    const double Cx  = 0.5*(C1x + C2x);
    const double Cy  = 0.5*(C1y + C2y);

    // anchor phase at state closer in z to zTarget
    const bool anchorAt1 = (abs(zTarget - z1) <= abs(zTarget - z2));
    const double xa = anchorAt1 ? x1 : x2;
    const double ya = anchorAt1 ? y1 : y2;
    const double za = anchorAt1 ? z1 : z2;

    const double phiA = atan2(ya - Cy, xa - Cx);
    const double phi  = phiA + sgn * (zTarget - za) / (R * tanL);

    return { Cx + R*cos(phi), Cy + R*sin(phi) };
}

// convenient wrapper for two GeomState objects
inline std::pair<double,double>
trackXYatZ(const GeomState& A, const GeomState& B, double zTarget){
    return trackXYatZ_fromTwoStates(
        A.x,A.y,A.z, A.px,A.py,A.pz,
        B.x,B.y,B.z, B.px,B.py,B.pz,
        zTarget
    );
}

int dimuon_fotoproduction_analysis(const string& filename, string outname_pdf, string outname_png, TString compact_file) {

    gStyle->SetOptStat(0);
    podio::ROOTReader reader;
    reader.openFile(filename);
    unsigned nEvents = reader.getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    det = &(dd4hep::Detector::getInstance());
    det->fromCompact(compact_file.Data());
    det->volumeManager();
    det->apply("DD4hepVolumeManager", 0, 0);

    dd4hep::rec::CellIDPositionConverter cellid_converter(*det);
    auto idSpec = det->readout("HcalEndcapNHits").idSpec();
    auto decoder = idSpec.decoder();
    const int slice_index = decoder->index("slice");
    if (slice_index < 0) {
        cerr << "ERROR: 'slice' field not found in cell ID spec!" << endl;
        return 1;
    }

    //double t_cm = det->constantAsDouble("HcalEndcapNPolystyreneThickness")* 0.1;

    constexpr int NBINS = 60; 
    constexpr int    Z_NBINS = 100;
    constexpr double Z_MIN_MM = -4500.0;
    constexpr double Z_MAX_MM = -3600.0;

    constexpr int    DXY_NBINS  = 120;
    constexpr double DXY_MIN_MM = -1000.0; 
    constexpr double DXY_MAX_MM =  1000.0;

    constexpr int    DR_NBINS  = 120;
    constexpr double DR_MIN_MM = 0.0;   
    constexpr double DR_MAX_MM = 400.0;

    constexpr double P_MIN_GEV = 0.0;  
    constexpr double P_MAX_GEV = 25.0;   
    constexpr int    P_NB  = 50;

    constexpr double E_MIN_GEV = 0.0;
    constexpr double E_MAX_GEV = 8.0;

    constexpr double DR_CUTS_CM[3] = {1.0, 3.0, 5.0};
    constexpr double MIP_ENERGY_GEV = 0.002; 
    constexpr double E_CUTS_GEV[3] = {0.5, 0.7, 1.0}; 
    constexpr double THRESH_MM = DR_CUTS_CM[2]*10;
    constexpr double E_THRESH_GEV = E_CUTS_GEV[2];   

    TH2D* hEtaPt = new TH2D("hEtaPt", "Muon #eta vs p_{T};#eta;p_{T} [GeV]", 100, -6., 6., 100, 0., 7.);

    TH2D* hx_Q2 = new TH2D("hx_Q2", "Muon x vs Q^{2}; x; Q^{2}[GeV^{2}]", 100, 1e-6, 1e-3, 100, 1e-2, 1.0);
    TH2D* hEta1_Eta2 = new TH2D("hEta1_Eta2", "Muon #eta_{+} vs #eta_{-}; #eta_{+} (PDG=-13); #eta_{-} (PDG=13)", 100, -6., 6., 100, -6., 6.);

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

    TH1D* hZ_proj = new TH1D("hZ_proj", "Muon track projections in nHCal; z [mm]; N", NBINS, Z_MIN_MM, Z_MAX_MM);

    TH1D* hZ_hits = new TH1D("hZ_hits", "Reconstructed hit z in nHCal; z [mm]; N", NBINS, Z_MIN_MM, Z_MAX_MM);
    
    TH3D* hDxDyZ_layer = new TH3D("hDxDyZ_layer","3D residuals (rec - proj): dx, dy vs z; dx [mm]; dy [mm]; z [mm]", NBINS, DXY_MIN_MM, DXY_MAX_MM, NBINS, DXY_MIN_MM, DXY_MAX_MM, NBINS, Z_MIN_MM, Z_MAX_MM); 
    
    TH2D* hDxDy_all = new TH2D("hDxDy_all",
                           "Residuals (rec - proj): dx vs dy (all layers); dx [mm]; dy [mm]",
                           DXY_NBINS, DXY_MIN_MM, DXY_MAX_MM, DXY_NBINS, DXY_MIN_MM, DXY_MAX_MM);
    
    TH2D* hDrZ_layer = new TH2D("hDrZ_layer","Radial residual (rec - proj) vs z; dr [mm]; z [mm]", NBINS, DR_MIN_MM, DR_MAX_MM, NBINS, Z_MIN_MM, Z_MAX_MM); 
    
    TH1D* hDr_all = new TH1D("hDr_all",
                         "Radial residual dr = #sqrt{dx^{2}+dy^{2}} (all layers); dr [mm]; N",
                         DR_NBINS, DR_MIN_MM, DR_MAX_MM);

    TH2D* hE_z = new TH2D("hE_z", "Reconstructed hit energy vs layer z; z_{layer} [mm]; E [GeV]", NBINS, Z_MIN_MM, Z_MAX_MM, NBINS, E_MIN_GEV, E_MAX_GEV);

    TH2D* hEsum_z = new TH2D("hEsum_z", "Layer energy sum vs z (reconstructed); z_{layer} [mm]; E_{sum} [GeV]",
                            NBINS, Z_MIN_MM, Z_MAX_MM, NBINS, E_MIN_GEV, E_MAX_GEV);

    TH1D* hP_all_mu = new TH1D("hP_all_mu", "MC muon momentum (muons in nHCal acceptance); p_{MC} [GeV]; N", P_NB, P_MIN_GEV, P_MAX_GEV);

    TH1D* hP_pass_dr[3] = {
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[0]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[0]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[1]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[1]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[2]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[2]), P_NB, P_MIN_GEV, P_MAX_GEV),
    };

    TH1D* hP_pass_Ecut[3] = {
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS_GEV[0]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS_GEV[0]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS_GEV[1]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS_GEV[1]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS_GEV[2]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS_GEV[2]), P_NB, P_MIN_GEV, P_MAX_GEV),
    };

    // To do
    TH1D* hP_pass_combo[3][3]; // [idr][ie]
    for (int idr=0; idr<3; ++idr) {
        for (int ie=0; ie<3; ++ie) {
            hP_pass_combo[idr][ie] = new TH1D(
                Form("hP_pass_dr%.0fcm_E<%.1fMIP", DR_CUTS_CM[idr], E_CUTS_GEV[ie]),
                Form("Accepted (dr<%.0f cm & E<%.1f MIP); p_{MC} [GeV]; N", DR_CUTS_CM[idr], E_CUTS_GEV[ie]),
                P_NB, P_MIN_GEV, P_MAX_GEV
            );
        }
    }

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

            double segMinDistance = 1e12;
            double segHitEnergy = 1e12;
            double t_cm;
            dd4hep::DDSegmentation::CellID segHitCellID = 0;

            GeomState A{}, B{};
            bool haveA = false, haveB = false;

            for (const auto& pt : points) {
                if (pt.system != 113) continue;
                hZ_proj->Fill(pt.position.z);                
                if (!haveA) {
                    A.x = pt.position.x; A.y = pt.position.y; A.z = pt.position.z;
                    A.px = pt.momentum.x; A.py = pt.momentum.y; A.pz = pt.momentum.z;
                    haveA = true;
                } else if (!haveB) {
                    B.x = pt.position.x; B.y = pt.position.y; B.z = pt.position.z;
                    B.px = pt.momentum.x; B.py = pt.momentum.y; B.pz = pt.momentum.z;
                    haveB = true;
                    break;
                }
            }

            if (!haveA || !haveB) {continue;}

            for (const auto& hit : hcalRec) {
                const auto& hpos = hit.getPosition();
                auto [X, Y] = trackXYatZ(A, B, hpos.z);          
                double dx = X - hpos.x;
                double dy = Y - hpos.y;
                const double dr = sqrt(dx*dx + dy*dy);
                hDxDyZ_layer->Fill(dx, dy, hpos.z);
                hDxDy_all->Fill(dx,dy);
                hDrZ_layer->Fill(dr,hpos.z);
                hDr_all->Fill(dr);

                if (dr < segMinDistance) {
                    segMinDistance = dr;
                    segHitEnergy = hit.getEnergy();
                    segHitCellID = hit.getCellID();
                }
            }

            t_cm = getPlasticThicknessMM(*det, decoder, segHitCellID, slice_index, 3);

            if (segMinDistance <= THRESH_MM) {
                cout << "[MATCH] muTag=" << ST.muTag
                    << " dr=" << segMinDistance << " mm (<= " << THRESH_MM << ")\n";
                if (ST.muTag == 1) {
                    m1_has_rec = true; 
                    for (int idr=0; idr<3; ++idr) if (segMinDistance < DR_CUTS_CM[idr] * 10) hP_pass_dr[idr]->Fill(v1.P());
                    for (int ie=0; ie<3; ++ie) if ( segHitEnergy < (E_CUTS_GEV[ie] * t_cm) ) hP_pass_Ecut[ie]->Fill(v1.P());
                }
                if (ST.muTag == 2) {
                    m2_has_rec = true; 
                    for (int idr=0; idr<3; ++idr) if (segMinDistance < DR_CUTS_CM[idr] * 10) hP_pass_dr[idr]->Fill(v2.P());
                    for (int ie=0; ie<3; ++ie) if ( segHitEnergy < (E_CUTS_GEV[ie] * t_cm) ) hP_pass_Ecut[ie]->Fill(v2.P());
                }
            } else {
                cout << "[INFO] no match <= " << THRESH_MM
                    << " mm; best distance (this seg) = " << segMinDistance
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

    TCanvas* canvas_sim = new TCanvas("canvas_sim", "Muon analysis", 1600, 800);
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

    TCanvas* canvas_rec = new TCanvas("canvas_rec", "Muon analysis rec", 1600, 800);
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

    canvas_rec->SaveAs(addPrefixAfterSlash(outname_png, "rec_").c_str());
    canvas_rec->SaveAs(addPrefixAfterSlash(outname_pdf, "rec_").c_str());

    TCanvas* canvas = new TCanvas("canvas", "Muon analysis extra", 1800, 1200); 
    canvas->Divide(3,4);                                                                 
    canvas->cd(1); gPad->SetGrid(); hZ_proj->Draw();                                     
    canvas->cd(2); gPad->SetGrid(); hZ_hits->Draw();                                     
    canvas->cd(3); gPad->SetGrid(); hDxDyZ_layer->Draw("COLZ");     
    canvas->cd(4); gPad->SetGrid(); hDxDy_all->Draw("COLZ"); 
    canvas->cd(5); gPad->SetGrid(); hDrZ_layer->Draw("COLZ");                             
    canvas->cd(6); gPad->SetGrid(); hDr_all->Draw("COLZ");                                      
    canvas->cd(7); gPad->SetGrid(); hE_z->Draw("COLZ");                                  
    canvas->cd(8); gPad->SetGrid(); hEsum_z->Draw("COLZ");                                
    canvas->cd(9); gPad->SetGrid(); hP_all_mu->SetLineColor(kBlack); hP_all_mu->Draw();  
    canvas->cd(10); gPad->SetGrid();                                                       
    TH1D* hEff_dr[3]; 
    for (int idr=0; idr<3; ++idr) { 
        hEff_dr[idr] = (TH1D*)hP_pass_dr[idr]->Clone( 
            Form("hEff_dr%.1fcm", DR_CUTS_CM[idr]) 
        ); 
        hEff_dr[idr]->SetTitle("Matching efficiency vs p_{MC} (dr < d cm); p_{MC} [GeV]; efficiency"); 
        hEff_dr[idr]->Divide(hP_all_mu); 
    } 

    hEff_dr[0]->SetLineColor(kBlue);     hEff_dr[0]->SetMinimum(0.0); hEff_dr[0]->SetMaximum(1.1); 
    hEff_dr[1]->SetLineColor(kRed);   
    hEff_dr[2]->SetLineColor(kGreen+2); 

    hEff_dr[0]->Draw("HIST");            
    hEff_dr[1]->Draw("HIST SAME");       
    hEff_dr[2]->Draw("HIST SAME");       
    { auto leg_dr = new TLegend(0.55,0.65,0.88,0.88);
        for(int idr=0; idr<3; ++idr)
        {leg_dr->AddEntry(hEff_dr[idr],Form("dr < %.1f cm", DR_CUTS_CM[idr]),"l");} 
        leg_dr->Draw(); 
    }  
                                                                   
    canvas->cd(11); gPad->SetGrid();                                                        
    TH1D* hEff_E[3]; 
    for (int ie=0; ie<3; ++ie) { 
        hEff_E[ie] = (TH1D*)hP_pass_Ecut[ie]->Clone( 
            Form("hEff_E%.1fcm", E_CUTS_GEV[ie]) 
        ); 
        hEff_E[ie]->SetTitle("Matching efficiency vs p_{MC} (E < k MIP); p_{MC} [GeV]; efficiency"); 
        hEff_E[ie]->Divide(hP_all_mu); 
    } 

    hEff_E[0]->SetLineColor(kBlue);     hEff_E[0]->SetMinimum(0.0); hEff_E[0]->SetMaximum(1.1); 
    hEff_E[1]->SetLineColor(kRed);   
    hEff_E[2]->SetLineColor(kGreen+2); 

    hEff_E[0]->Draw("HIST");            
    hEff_E[1]->Draw("HIST SAME");       
    hEff_E[2]->Draw("HIST SAME");       
    { auto leg_E = new TLegend(0.55,0.65,0.88,0.88);
        for(int ie=0; ie<3; ++ie)
        {leg_E->AddEntry(hEff_E[ie],Form("E < %.1f MIP", E_CUTS_GEV[ie]),"l");}  
        leg_E->Draw(); 
    }  
    canvas->SaveAs(addPrefixAfterSlash(outname_png, "extra_").c_str());              
    canvas->SaveAs(addPrefixAfterSlash(outname_pdf, "extra_").c_str());              

    return 0;
}