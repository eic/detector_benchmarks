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
#include <TLorentzVector.h>  
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>   

#include "TROOT.h"
#include "TRandom.h"
#include "TH3.h"

#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Shapes.h"
#include <DD4hep/Detector.h>
#include <DD4hep/VolumeManager.h>
#include <DD4hep/Printout.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>

#include "podio/Frame.h"
#include "podio/CollectionBase.h"
#include "podio/ROOTReader.h"
#include "podio/CollectionIDTable.h"
#include "podio/ObjectID.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHit.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollection.h"

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/Cluster.h"
#include "edm4eic/ClusterData.h"

#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitCollection.h"

#include "edm4eic/InclusiveKinematicsCollection.h"
#include "edm4eic/InclusiveKinematics.h"
#include "edm4hep/utils/kinematics.h"
#include "edm4hep/utils/vector_utils.h"
#include "edm4eic/vector_utils_legacy.h"

#include "edm4eic/Track.h"
#include "edm4eic/TrackSegment.h"
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

void MakeLogBins(double *array, int nbins, double binLo, double binHi)
{
    double logMin = log10(binLo);
    double logMax = log10(binHi);
    double binWidth = (logMax - logMin) / nbins;

    for (int i = 0; i <= nbins; ++i) {
    	array[i] = pow(10, logMin + i * binWidth);
    }
}

TLine* makeLine(double x1, double y1, double x2, double y2) {
    TLine* l = new TLine(x1, y1, x2, y2);
    l->SetLineColor(kRed);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw("same");
    return l;
}

static TGraph* makeEffGraph_B(TH1D* h_sel, TH1D* h_all,
                              const char* name,
                              int minAll = 2,
                              bool logx = false)
{
    if (!h_sel || !h_all) return nullptr;

    std::unique_ptr<TH1D> h_eff((TH1D*)h_sel->Clone(Form("tmp_%s", name ? name : "eff")));
    h_eff->SetDirectory(nullptr);
    h_eff->Sumw2();
    h_eff->Divide(h_sel, h_all, 1.0, 1.0, "B");

    const int nb = h_all->GetNbinsX();

    std::unique_ptr<TGraphAsymmErrors> g_log;
    std::unique_ptr<TGraphErrors>      g_lin;

    if (logx) g_log = std::make_unique<TGraphAsymmErrors>();
    else      g_lin = std::make_unique<TGraphErrors>();

    int ip = 0;
    for (int b = 1; b <= nb; ++b) {
        const double all = h_all->GetBinContent(b);
        if (all < minAll) continue;

        const double xl = h_all->GetXaxis()->GetBinLowEdge(b);
        const double xh = h_all->GetXaxis()->GetBinUpEdge(b);

        const double y  = 100.0 * h_eff->GetBinContent(b);
        const double ey = 100.0 * h_eff->GetBinError(b);

        if (logx) {
            if (xl <= 0.0 || xh <= 0.0) continue; 
            const double x  = sqrt(xl * xh);
            const double exl = x - xl;
            const double exh = xh - x;

            g_log->SetPoint(ip, x, y);
            g_log->SetPointError(ip, exl, exh, ey, ey);
        } else {
            const double x  = h_all->GetBinCenter(b);
            const double ex = 0.5 * (xh - xl);

            g_lin->SetPoint(ip, x, y);
            g_lin->SetPointError(ip, ex, ey);
        }
        ++ip;
    }

    TGraph* g = logx ? (TGraph*)g_log.release() : (TGraph*)g_lin.release();
    g->SetName(Form("g_%s", name ? name : "eff"));
    return g;
}

inline void drawEffPanel(TH1D* h_all,
                         TH1D* const h_in[3],
                         const char* title,
                         const char* xlabel,
                         bool logx)
{
    gPad->SetGrid();
    if (logx) gPad->SetLogx();

    auto* mg = new TMultiGraph();
    auto* leg = new TLegend(0.63, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    Color_t colors[3]  = {kBlue, kRed, kBlack};
    Style_t markers[3] = {20, 21, 22};

    int added = 0;

    for (int i = 0; i < 3; ++i) {
        if (!h_all || !h_in[i]) continue;
        auto* g = makeEffGraph_B(h_in[i], h_all, Form("eff_%d", i), 2, logx);
        if (!g || g->GetN() == 0) continue;

        g->SetMarkerColor(colors[i]);
        g->SetMarkerStyle(markers[i]);
        g->SetMarkerSize(1.0);
        g->SetLineColor(kBlack);

        mg->Add(g, "PE");
        leg->AddEntry(g, Form("%d mu-in nHCAL", i), "p");
        ++added;
    }

    if (h_all && h_in[0] && h_in[1] && h_in[2]) {
        std::unique_ptr<TH1D> h_sum((TH1D*)h_in[0]->Clone("h_sum"));
        h_sum->SetDirectory(nullptr);
        h_sum->Add(h_in[1]);
        h_sum->Add(h_in[2]);

        auto* gsum = makeEffGraph_B(h_sum.get(), h_all, "eff_sum", 2, logx);
        if (gsum && gsum->GetN() > 0) {
            gsum->SetMarkerStyle(25);
            gsum->SetMarkerColor(kMagenta + 1);
            gsum->SetFillColor(kMagenta + 1);
            gsum->SetLineColor(kBlack);
            gsum->SetMarkerSize(1.0);

            mg->Add(gsum, "PE");
            leg->AddEntry(gsum, "All eligible", "p");
            ++added;
        }
    }

    mg->SetTitle(Form("%s;%s;geom. acc. [%%]", title ? title : "", xlabel ? xlabel : ""));
    mg->Draw("A");

    gPad->Update();
    if (mg->GetHistogram()) {
        double ymax = mg->GetHistogram()->GetMaximum();
        mg->GetHistogram()->SetMaximum(ymax * 1.35); 
    }
    gPad->Modified();
    gPad->Update();

    if (logx) {
        auto* ax = mg->GetXaxis();
        if (ax) {
            double xmin = ax->GetXmin(), xmax = ax->GetXmax();
            if (xmin <= 0) xmin = 1e-12;
            mg->GetXaxis()->SetLimits(xmin, xmax);
        }
    }

    leg->Draw();
}

inline TVector3 getPlasticDimensionsCM(dd4hep::Detector& det,
                                    const dd4hep::DDSegmentation::BitFieldCoder* dec,
                                    dd4hep::DDSegmentation::CellID cid,
                                    int slice_idx,
                                    int plastic_slice_value = 3)
{
    const double NaN = numeric_limits<double>::quiet_NaN();
    TVector3 dims(NaN, NaN, NaN);
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

        dims.SetXYZ(2.0 * box->GetDX(),
                    2.0 * box->GetDY(),
                    2.0 * box->GetDZ());
        return dims;
    } catch (const exception& e) {
        cerr << "[WARN] getPlasticThicknessMM: " << e.what() << " (cellID=" << cid << ")\n";
        return dims;
    }
}

inline double getPlasticCenterZ_cm(
    const dd4hep::DDSegmentation::BitFieldCoder* decoder,
    const dd4hep::rec::CellIDPositionConverter& cellid_converter,
    dd4hep::DDSegmentation::CellID cid,
    int slice_index,
    int plastic_slice_value = 3)
{
    try {
        if (!decoder) throw std::runtime_error("decoder==nullptr");
        if (decoder->get(cid, slice_index) != plastic_slice_value)
            throw std::runtime_error("cell is not plastic (slice mismatch)");
        return cellid_converter.position(cid).z();
    } catch (const std::exception& e) {
        cerr << "[WARN] getPlasticThicknessCM: " << e.what()
                << " (cellID=" << cid << ")\n";
        return std::numeric_limits<double>::quiet_NaN();;
    }
}

struct GeomState {
    double x, y, z;     
};

inline double hypot2(double a, double b){ return sqrt(a*a + b*b); }
inline int sgnD(double v){ return (v>0) - (v<0); }

inline pair<double,double>
trackXYatZ_fromTwoStates(double x1,double y1,double z1,
                         double x2,double y2,double z2,
                         double zTarget)
{
    constexpr double EPSZ = 1e-12;

    const double dz = z2 - z1;
    if (abs(dz) < EPSZ) {
        const double d1 = abs(zTarget - z1);
        const double d2 = abs(zTarget - z2);
        if (d1 < d2) return {x1, y1};
        if (d2 < d1) return {x2, y2};
        return {0.5*(x1+x2), 0.5*(y1+y2)};
    }

    const double t = (zTarget - z1) / dz;
    const double x = x1 + t * (x2 - x1);
    const double y = y1 + t * (y2 - y1);
    return {x, y};
}

inline pair<double,double>
trackXYatZ(const GeomState& A, const GeomState& B, double zTarget){
    return trackXYatZ_fromTwoStates(
        A.x, A.y, A.z, B.x, B.y, B.z,
        zTarget
    );
}


int pion_rejection_analysis(const string& filename, string outname_pdf, string outname_png, TString compact_file) {

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

    //double thickness_cm = det->constantAsDouble("HcalEndcapNPolystyreneThickness")* 0.1;

    constexpr int NBINS = 60; 

    constexpr double Q2_MIN = 1e-2;
    constexpr double Q2_MAX = 1.0;
    constexpr double X_MIN = 1e-6;
    constexpr double X_MAX = 1e-3;
    constexpr double Y_MIN = 0.0;
    constexpr double Y_MAX = 1.0;
    constexpr double W_MIN = 0.0;
    constexpr double W_MAX = 160.0;

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
    constexpr int    P_NBINS  = 50;

    constexpr double E_MIN_GEV = 2e-2;
    constexpr double E_MAX_GEV = 10.0;

    constexpr int    SIZE = 3;
    constexpr double MIP_ENERGY_GEV = 0.002; 
    constexpr double E_CUTS[SIZE] = {1.5, 1.7, 2.0}; 
    constexpr double E_THRESH = E_CUTS[2]; 
    constexpr double LAYER_PROC[SIZE] = {0.6, 0.7, 0.8};

    double t_cm;
    double DR_CUTS_CM[SIZE] = {7.0, 10.0, 13.0};
    double DR_THRESH_MM = 10*DR_CUTS_CM[2];
    int LAYER_MAX = 10;
    int LAYER_CUTS[SIZE] = {5, 6, 7}; 
    int LAYER_THRESH = LAYER_CUTS[0];

    vector<double> Xbins(NBINS+1), Q2bins(NBINS+1), Ebins(NBINS+1);
    MakeLogBins(Q2bins.data(), NBINS, Q2_MIN, Q2_MAX);
    MakeLogBins(Xbins.data(), NBINS, X_MIN, X_MAX);
    MakeLogBins(Ebins.data(), NBINS, E_MIN_GEV, E_MAX_GEV);

    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetLabelSize(0.06, "XYZ");
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);
    gROOT->ForceStyle();

    TH2D* hEtaPt = new TH2D("hEtaPt", "Pion #eta vs p_{T};#eta;p_{T} [GeV]", 100, -6., 6., 100, 0., 7.);

    TH1D* hZ_proj = new TH1D("hZ_proj", "Pion track projections in nHCal; z [mm]; N_{proj}", NBINS, Z_MIN_MM, Z_MAX_MM);

    TH1D* hZ_hits = new TH1D("hZ_hits", "Reconstructed hit z in nHCal; z [mm]; N_{hits}", NBINS, Z_MIN_MM, Z_MAX_MM);
    
    TH3D* hDxDyZ_layer = new TH3D("hDxDyZ_layer","3D residuals (rec - proj): dx, dy vs z; dx [mm]; dy [mm]; z [mm]", NBINS, DXY_MIN_MM, DXY_MAX_MM, NBINS, DXY_MIN_MM, DXY_MAX_MM, NBINS, Z_MIN_MM, Z_MAX_MM); 
    
    TH2D* hDxDy_all = new TH2D("hDxDy_all",
                           "Residuals (rec - proj): dx vs dy (all layers); dx [mm]; dy [mm]",
                           DXY_NBINS, DXY_MIN_MM, DXY_MAX_MM, DXY_NBINS, DXY_MIN_MM, DXY_MAX_MM);
    
    TH2D* hDrZ_layer = new TH2D("hDrZ_layer","Radial residual (rec - proj) vs z; dr [mm]; z [mm]", NBINS, DR_MIN_MM, DR_MAX_MM, NBINS, Z_MIN_MM, Z_MAX_MM); 
    
    TH1D* hDr_all = new TH1D("hDr_all",
                         "Radial residual dr = #sqrt{dx^{2}+dy^{2}} (all layers); dr [mm]; N",
                         DR_NBINS, DR_MIN_MM, DR_MAX_MM);

    TH2D* hE_z = new TH2D("hE_z", "Reconstructed hit energy vs layer z; z_{layer} [mm]; E [GeV]", NBINS, Z_MIN_MM, Z_MAX_MM, NBINS, Ebins.data());

    TH2D* hEsum_z = new TH2D("hEsum_z", "Layer energy sum vs z (reconstructed); z_{layer} [mm]; E_{sum} [GeV]",
                            NBINS, Z_MIN_MM, Z_MAX_MM, NBINS, Ebins.data());

    TH1D* hE = new TH1D("hE", "Reconstructed hit energy; E [GeV]; N_{hits}", 10000, 0.0, 1.0);

    TH1D* hEsum = new TH1D("hEsum", "Layer energy sum (reconstructed); E_{sum} [GeV]; N_{energy}", 10000, 0.0, 1.0);

    TH1D* hP_all_pion = new TH1D("hP_all_pion", "MC pion momentum (pions in nHCal acceptance); p_{MC} [GeV]; N_{pions}", P_NBINS, P_MIN_GEV, P_MAX_GEV);

    TH1D* hP_pass_dr[3] = {
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[0]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[0]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[1]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[1]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[2]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[2]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
    };

    TH1D* hP_pass_ECut[3] = {
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS[0]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS[0]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS[1]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS[1]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS[2]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS[2]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
    };

    TH1D* hP_pass_LayerCut[3] = {
        new TH1D(Form("hP_pass_Layer< %d layers", LAYER_CUTS[0]), Form("Accepted (Layer < %d layers); p_{MC} [GeV]; N",LAYER_CUTS[0]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_Layer< %d layers", LAYER_CUTS[1]), Form("Accepted (Layer < %d layers); p_{MC} [GeV]; N",LAYER_CUTS[1]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_Layer< %d layers", LAYER_CUTS[2]), Form("Accepted (Layer < %d layers); p_{MC} [GeV]; N",LAYER_CUTS[2]), P_NBINS, P_MIN_GEV, P_MAX_GEV),
    };

    for (unsigned ev = 0; ev < nEvents; ev++) {
        auto frameData = reader.readNextEntry(podio::Category::Event);
        if (!frameData) continue;
        podio::Frame frame(std::move(frameData));

        auto& mcCol  = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        auto& recParts = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedParticles");
        auto& projSegs = frame.get<edm4eic::TrackSegmentCollection>("CalorimeterTrackProjections");
        auto& hcalRec  = frame.get<edm4eic::CalorimeterHitCollection>("HcalEndcapNRecHits");
        auto& assocCol = frame.get<edm4eic::MCRecoParticleAssociationCollection>("ReconstructedParticleAssociations");

        if (!mcCol.isValid()) continue;

        vector<edm4hep::MCParticle> vPions;
        vector<TLorentzVector> vLorentzPions;
        for (const auto& p : mcCol) {
            if (p.getPDG() != -211 || p.getGeneratorStatus() == 0) continue;
            vPions.push_back(p); 
            TLorentzVector v(p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy());
            vLorentzPions.push_back(v); 
            hEtaPt->Fill(v.Eta(), v.Pt());
            if(inNHCal(v.Eta())) {hP_all_pion->Fill(v.P());}
        }

        vector<edm4eic::ReconstructedParticle> matchedRecos;
        auto find_associated_reco = [&](const edm4hep::MCParticle& mc)->void {
            try {
                if (!assocCol.isValid() || assocCol.empty()) return;
                
                const uint32_t mc_idx = mc.getObjectID().index;
                
                for (const auto& assoc : assocCol) {
                    if (assoc.getSimID() == mc_idx) {
                        uint32_t ridx = assoc.getRecID();
                        if (!recParts.isValid() || ridx >= recParts.size()) continue;
                        auto reco = recParts.at(ridx);
                        if (reco.isAvailable()) {
                            matchedRecos.push_back(reco);
                        }
                    }
                }
            } catch (...) {}
        };

        if (!assocCol.isValid() || assocCol.empty()) continue;
        for(const auto&  p: vPions) if (p.getPDG() == -211) find_associated_reco(p);

        if (!recParts.isValid() || !projSegs.isValid() || !hcalRec.isValid() || recParts.empty() || projSegs.empty() || hcalRec.empty()) continue;
        
        set<int> uniqueCentersZ10;
        map<double, double> layerData;
        for (const auto& hit : hcalRec) {    
            double z = hit.getPosition().z;
            double zBin = round(z);
            uint64_t cid = hit.getCellID();
            double zCenter_mm = 10 * getPlasticCenterZ_cm(decoder, cellid_converter, cid, slice_index, /*plastic_slice_value=*/3);
            if (!std::isnan(zCenter_mm)) {
                int z10 = lround(zCenter_mm * 10.0); 
                uniqueCentersZ10.insert(z10); 
            }
            layerData[zBin] += hit.getEnergy(); 
            hZ_hits->Fill(z);                
            hE_z->Fill(z, hit.getEnergy()); 
            hE->Fill(hit.getEnergy());
        }

        vector<double> layerCentersZ;
        layerCentersZ.reserve(uniqueCentersZ10.size());
        for (int z10 : uniqueCentersZ10) layerCentersZ.push_back(z10 / 10.0);
        if(layerCentersZ.size() > LAYER_MAX) LAYER_MAX = layerCentersZ.size();

        for(size_t n = 0; n < SIZE; n++) LAYER_CUTS[n] = static_cast<int>(LAYER_MAX*LAYER_PROC[n]);
        LAYER_THRESH = LAYER_CUTS[0];

        for (const auto& [zValue, sumEnergy] : layerData) {hEsum_z->Fill(zValue, sumEnergy); hEsum->Fill(sumEnergy);}

        vector<edm4eic::Track> allTracks;
        for (const auto& R : matchedRecos) {
            for (const auto& tr : R.getTracks()) {
                if (tr.isAvailable()) allTracks.push_back(tr);
            }
        }
        vector<edm4eic::TrackSegment> segsTagged;
        for (const auto& seg : projSegs) {
            auto linkedTr = seg.getTrack();
            if (!linkedTr.isAvailable()) continue;
            for (const auto& TT : allTracks) {
                if (linkedTr.getObjectID() == TT.getObjectID()) {
                    segsTagged.push_back(seg);
                    break;
                }
            }
        }

        for (size_t s = 0; s < segsTagged.size(); ++s) {

            vector<double> segMinDistance(LAYER_MAX, std::numeric_limits<double>::infinity());
            vector<double> segHitEnergy(LAYER_MAX, std::numeric_limits<double>::quiet_NaN());
            vector<double> thickness_cm(LAYER_MAX, std::numeric_limits<double>::quiet_NaN());
            vector<int> count_DrCuts(SIZE, 0);
            vector<int> count_ECuts(SIZE, 0);

            GeomState A{}, B{};
            bool haveA = false, haveB = false;

            auto points = segsTagged[s].getPoints();

            for (const auto& pt : points) {
                if (pt.system != 113) continue;
                hZ_proj->Fill(pt.position.z);                
                if (!haveA) {
                    A.x = pt.position.x; A.y = pt.position.y; A.z = pt.position.z;
                    haveA = true;
                } else if (!haveB) {
                    B.x = pt.position.x; B.y = pt.position.y; B.z = pt.position.z;
                    haveB = true; 
                    break;
                }
            }

            if (!haveA || !haveB) {continue;}

            for (size_t i = 0; i < LAYER_MAX; ++i) {
                double best_dr_in_layer = 1e10;
                double best_E_in_layer;
                double partLayerEnergySum = 0;
                double ratio_HitPartLayerEnergy = 0;
                dd4hep::DDSegmentation::CellID best_cid_in_layer;

                double zc = layerCentersZ[i];
                auto [X, Y] = trackXYatZ(A, B, zc);

                for (const auto& hit : hcalRec) {
                    const auto& hp = hit.getPosition();

                    if (fabs(hp.z - zc) > 5.0 /*mm*/) continue;

                    const double dx = X - hp.x;
                    const double dy = Y - hp.y;
                    const double dr = sqrt(dx*dx + dy*dy);
                    if(dr < 30*10 /*mm*/) partLayerEnergySum += hit.getEnergy();

                    hDxDyZ_layer->Fill(dx, dy, hp.z);
                    hDxDy_all->Fill(dx, dy);
                    hDrZ_layer->Fill(dr, hp.z);
                    hDr_all->Fill(dr);

                    if (dr < best_dr_in_layer) {
                        best_dr_in_layer = dr;
                        best_E_in_layer  = hit.getEnergy();
                        best_cid_in_layer = hit.getCellID();
                    }
                }

                if (best_dr_in_layer < DR_THRESH_MM) {
                    segMinDistance[i] = best_dr_in_layer;
                    segHitEnergy[i]   = best_E_in_layer;
                    thickness_cm[i]   = getPlasticDimensionsCM(*det, decoder, best_cid_in_layer, slice_index, 3).z();
                    t_cm = thickness_cm[0];
                }
                ratio_HitPartLayerEnergy += segHitEnergy[i]/partLayerEnergySum;
                if (std::isnan(ratio_HitPartLayerEnergy) || fabs(ratio_HitPartLayerEnergy - 1) >= 0.2) continue; 
                for(size_t j = 0; j < SIZE; ++j){
                    if (segMinDistance[i] < DR_CUTS_CM[j] * 10 /*mm*/){++count_DrCuts[j];}
                    if (segHitEnergy[i] < (E_CUTS[j] * MIP_ENERGY_GEV * thickness_cm[i] * 100)){++count_ECuts[j];}
                }
            }
            for (int j = 0; j < SIZE; ++j){
                if (count_DrCuts[j] > LAYER_THRESH) hP_pass_dr[j]->Fill(vLorentzPions[s].P());
                if (count_ECuts[j] > LAYER_THRESH) hP_pass_ECut[j]->Fill(vLorentzPions[s].P());
                if (count_DrCuts[SIZE-1] > LAYER_CUTS[j]) hP_pass_LayerCut[j]->Fill(vLorentzPions[s].P());
            }
        }   
    }

    TCanvas* c = new TCanvas("c", "Pion analysis", 1600, 1000);
    c->Divide(2,2);

    c->cd(1); hEtaPt->Draw("COLZ");
  
    c->SaveAs(outname_pdf.c_str());
    c->SaveAs(outname_png.c_str());

    TCanvas* c_hitTrack = new TCanvas("c_hitTrack", "Pion hit and track analysis", 1600, 1000);
    c_hitTrack->Divide(2,2);                                                                 
    c_hitTrack->cd(1); gPad->SetGrid(); hZ_proj->Draw();                                     
    c_hitTrack->cd(2); gPad->SetGrid(); hZ_hits->Draw();
    c_hitTrack->SaveAs(addPrefixAfterSlash(outname_png, "hitTrack_").c_str());              
    c_hitTrack->SaveAs(addPrefixAfterSlash(outname_pdf, "hitTrack_").c_str()); 

    TCanvas* c_hitTrackDistance = new TCanvas("c_hitTrackDistance", "Pion hit-track distance analysis", 1600, 1000);
    c_hitTrackDistance->Divide(2,2); 
    c_hitTrackDistance->cd(1); gPad->SetGrid(); hDxDyZ_layer->Draw("COLZ");     
    c_hitTrackDistance->cd(2); gPad->SetGrid(); hDxDy_all->Draw("COLZ"); 
    c_hitTrackDistance->cd(3); gPad->SetGrid(); hDrZ_layer->Draw("COLZ");                             
    c_hitTrackDistance->cd(4); gPad->SetGrid(); hDr_all->Draw("COLZ");
    c_hitTrackDistance->SaveAs(addPrefixAfterSlash(outname_png, "hitTrackDistance_").c_str());              
    c_hitTrackDistance->SaveAs(addPrefixAfterSlash(outname_pdf, "hitTrackDistance_").c_str()); 

    TCanvas* c_Eff = new TCanvas("c_Eff", "Pion efficiency analysis", 1600, 1000);
    c_Eff->Divide(2,2);                        
    c_Eff->cd(1); gPad->SetGrid(); hP_all_pion->SetLineColor(kBlack); hP_all_pion->Draw();  
    c_Eff->cd(2); gPad->SetGrid();                                                       
    TH1D* hEff_dr[3]; 
    for (int idr=0; idr<3; ++idr) { 
        hEff_dr[idr] = (TH1D*)hP_pass_dr[idr]->Clone( 
            Form("hEff_dr%.1fcm", DR_CUTS_CM[idr]) 
        ); 
        hEff_dr[idr]->SetDirectory(nullptr);
        hEff_dr[idr]->SetTitle(Form("Efficiency vs p_{MC} (Layers > %d); p_{MC} [GeV]; efficiency", LAYER_THRESH)); 
        hEff_dr[idr]->Divide(hP_pass_dr[idr], hP_all_pion, 1, 1, "B");
    } 

    hEff_dr[0]->SetLineColor(kBlue); hEff_dr[0]->SetMinimum(0.0); hEff_dr[0]->SetMaximum(0.3); 
    hEff_dr[1]->SetLineColor(kRed);   
    hEff_dr[2]->SetLineColor(kGreen+2); 

    hEff_dr[0]->Draw("HIST E");            
    hEff_dr[1]->Draw("HIST E SAME");       
    hEff_dr[2]->Draw("HIST E SAME");       
    { auto leg_dr = new TLegend(0.60,0.70,0.88,0.88);
        for(int idr=0; idr<3; ++idr)
        {leg_dr->AddEntry(hEff_dr[idr],Form("dr < %.1f cm", DR_CUTS_CM[idr]),"l");} 
        leg_dr->Draw(); 
    }  
                                                                   
    c_Eff->cd(3); gPad->SetGrid();                                                        
    TH1D* hEff_E[3]; 
    for (int ie=0; ie<3; ++ie) { 
        hEff_E[ie] = (TH1D*)hP_pass_ECut[ie]->Clone( 
            Form("hEff_E%.1fcm", E_CUTS[ie]) 
        ); 
        hEff_E[ie]->SetDirectory(nullptr);
        hEff_E[ie]->SetTitle(Form("Efficiency vs p_{MC} (Layers > %d, dr < %.1f cm); p_{MC} [GeV]; efficiency",LAYER_THRESH, DR_THRESH_MM/10)); 
        hEff_E[ie]->Divide(hP_pass_ECut[ie], hP_all_pion, 1, 1, "B");
    } 

    hEff_E[0]->SetLineColor(kBlue); hEff_E[0]->SetMinimum(0.0); hEff_E[0]->SetMaximum(0.3); 
    hEff_E[1]->SetLineColor(kRed);   
    hEff_E[2]->SetLineColor(kGreen+2); 

    hEff_E[0]->Draw("HIST E");            
    hEff_E[1]->Draw("HIST E SAME");       
    hEff_E[2]->Draw("HIST E SAME");       
    { auto leg_E = new TLegend(0.60,0.70,0.88,0.88);
        for(int ie=0; ie<3; ++ie)
        {leg_E->AddEntry(hEff_E[ie],Form("E < %.1f MIP", E_CUTS[ie]),"l");}  
        leg_E->Draw(); 
    }
    
    c_Eff->cd(4); gPad->SetGrid();
    TH1D* hEff_Layer[3]; 
    for (int il=0; il<3; ++il) { 
        hEff_Layer[il] = (TH1D*)hP_pass_LayerCut[il]->Clone( 
            Form("hEff_Layer%d", LAYER_CUTS[il]) 
        ); 
        hEff_Layer[il]->SetDirectory(nullptr);
        hEff_Layer[il]->SetTitle(Form("Efficiency vs p_{MC} (dr < %.1f cm); p_{MC} [GeV]; efficiency", DR_THRESH_MM/10)); 
        hEff_Layer[il]->Divide(hP_pass_LayerCut[il], hP_all_pion, 1, 1, "B");
    } 

    hEff_Layer[0]->SetLineColor(kBlue); hEff_Layer[0]->SetMinimum(0.0); hEff_Layer[0]->SetMaximum(0.3); 
    hEff_Layer[1]->SetLineColor(kRed);   
    hEff_Layer[2]->SetLineColor(kGreen+2); 

    hEff_Layer[0]->Draw("HIST E");            
    hEff_Layer[1]->Draw("HIST E SAME");       
    hEff_Layer[2]->Draw("HIST E SAME");       
    { auto leg_Layer = new TLegend(0.60,0.70,0.88,0.88);
        for(int il=0; il<3; ++il)
        {leg_Layer->AddEntry(hEff_Layer[il],Form("L > %d Layers", LAYER_CUTS[il]),"l");}  
        leg_Layer->Draw(); 
    }  
    c_Eff->SaveAs(addPrefixAfterSlash(outname_png, "matching_efficiency_").c_str());              
    c_Eff->SaveAs(addPrefixAfterSlash(outname_pdf, "matching_efficiency_").c_str()); 
    
    TCanvas* c_hEnergy = new TCanvas("c_hEnergy", "Pion hit energy analysis", 1600, 1000);
    c_hEnergy->Divide(2,2); 

    c_hEnergy->cd(1); gPad->SetGrid(); gPad->SetLogy(); hE_z->GetYaxis()->SetMoreLogLabels(); hE_z->Draw("COLZ");
    TProfile* pE_z = hE_z->ProfileX("pE_z"); pE_z->SetLineWidth(3); pE_z->SetLineColor(kRed); pE_z->SetMarkerColor(kRed);
                    pE_z->SetDirectory(nullptr); pE_z->Draw("SAME");

    const double Ecut_GeV = MIP_ENERGY_GEV * t_cm * 100;
    TLine* cut = new TLine(hE_z->GetXaxis()->GetXmin(), Ecut_GeV, hE_z->GetXaxis()->GetXmax(), Ecut_GeV);
                cut->SetLineColor(kRed+1); cut->SetLineStyle(2); cut->SetLineWidth(2);
                cut->Draw("SAME"); 
    auto* leg = new TLegend(0.58, 0.75, 0.84, 0.88); leg->SetTextSize(0.04);leg->AddEntry(cut, Form("E_{cut} = %.3f GeV", Ecut_GeV), "l"); 
                leg->Draw();    

    c_hEnergy->cd(2); gPad->SetGrid(); gPad->SetLogy(); hEsum_z->GetYaxis()->SetMoreLogLabels(); hEsum_z->Draw("COLZ");
    TProfile* pEsum_z = hEsum_z->ProfileX("pEsum_z"); pEsum_z->SetLineWidth(3); pEsum_z->SetLineColor(kRed); pEsum_z->SetMarkerColor(kRed);
                    pEsum_z->SetDirectory(nullptr); pEsum_z->Draw("SAME"); 
    c_hEnergy->cd(3); gPad->SetGrid(); hE->Draw();
    c_hEnergy->cd(4); gPad->SetGrid(); hEsum->Draw();
    
    c_hEnergy->SaveAs(addPrefixAfterSlash(outname_png, "hit_energy_").c_str());              
    c_hEnergy->SaveAs(addPrefixAfterSlash(outname_pdf, "hit_energy_").c_str()); 

    delete c;
    delete c_hitTrack;
    delete c_hitTrackDistance;
    delete c_Eff;
    delete c_hEnergy;

    return 0;
}