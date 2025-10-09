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

inline void setPadMargins(double left = 0.15,
                          double right = 0.05,
                          double bottom = 0.15,
                          double top = 0.10)
{
    if (!gPad) return;
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bottom);
    gPad->SetTopMargin(top);
}

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

static TGraphErrors* makeEffGraph_B(TH1D* h_sel, TH1D* h_all,
                                    const char* name,
                                    int minAll = 2)
{
    if (!h_sel || !h_all) return nullptr;

    std::unique_ptr<TH1D> h_eff((TH1D*)h_sel->Clone(Form("tmp_%s", name ? name : "eff")));
    h_eff->SetDirectory(nullptr);
    h_eff->Sumw2();                         
    h_eff->Divide(h_sel, h_all, 1.0, 1.0, "B");

    const int nb = h_all->GetNbinsX();
    vector<double> xs, ys, exs, eys;
    xs.reserve(nb); ys.reserve(nb); exs.reserve(nb); eys.reserve(nb);

    for (int b = 1; b <= nb; ++b) {
        const double all = h_all->GetBinContent(b);
        if (all < minAll) continue;

        const double x   = h_all->GetBinCenter(b);
        const double ex  = 0.5 * h_all->GetBinWidth(b);      
        const double y   = 100.0 * h_eff->GetBinContent(b);  
        const double ey  = 100.0 * h_eff->GetBinError(b);   

        xs.push_back(x);   ys.push_back(y);
        exs.push_back(ex); eys.push_back(ey);
    }

    auto* g = new TGraphErrors((int)xs.size(), xs.data(), ys.data(), exs.data(), eys.data());
    g->SetName(Form("g_%s", name ? name : "eff"));
    return g;
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

        auto* g = makeEffGraph_B(h_in[i], h_all, Form("eff_%d", i), 2);
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

        auto* gsum = makeEffGraph_B(h_sum.get(), h_all, "eff_sum", 2);
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

    if (added == 0 && h_all) {
        TH1D* hframe = new TH1D("hframe_tmp", "", 10,
                                h_all->GetXaxis()->GetXmin(),
                                h_all->GetXaxis()->GetXmax());
        hframe->SetMinimum(0.0);
        hframe->SetMaximum(110.0);
        hframe->GetXaxis()->SetTitle(xlabel ? xlabel : "");
        hframe->GetYaxis()->SetTitle("geom. acc. [%]");
        hframe->Draw();
    }

    return make_tuple(mg, leg, logx);
}

inline double getPlasticThicknessCM(dd4hep::Detector& det,
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

    //double thickness_cm = det->constantAsDouble("HcalEndcapNPolystyreneThickness")* 0.1;

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
    constexpr double E_MAX_GEV = 6.0;

    constexpr int    SIZE = 3;
    constexpr double DR_CUTS_CM[SIZE] = {7.0, 10.0, 13.0};
    constexpr double DR_THRESH_MM = DR_CUTS_CM[2]*10;
    constexpr double MIP_ENERGY_GEV = 0.002; 
    constexpr double E_CUTS_GEV[SIZE] = {1.5, 1.7, 2.0}; 
    constexpr double E_THRESH_GEV = E_CUTS_GEV[2]; 
    constexpr int    LAYER_CUTS[SIZE] = {5, 6, 7};
    constexpr int    LAYER_THRESH = LAYER_CUTS[2];

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

    TH1D* hE = new TH1D("hE", "Reconstructed hit energy; E [GeV]; N", 10000, 0.0, 1.0);

    TH1D* hEsum = new TH1D("hEsum", "Layer energy sum (reconstructed); E_{sum} [GeV]; N", 10000, 0.0, 1.0);

    TH1D* hP_all_mu = new TH1D("hP_all_mu", "MC muon momentum (muons in nHCal acceptance); p_{MC} [GeV]; N", P_NB, P_MIN_GEV, P_MAX_GEV);

    TH1D* hP_pass_dr[3] = {
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[0]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[0]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[1]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[1]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_dr%.1fcm",DR_CUTS_CM[2]),  Form("Accepted (dr<%.1fcm);  p_{MC} [GeV]; N",DR_CUTS_CM[2]), P_NB, P_MIN_GEV, P_MAX_GEV),
    };

    TH1D* hP_pass_ECut[3] = {
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS_GEV[0]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS_GEV[0]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS_GEV[1]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS_GEV[1]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_E< %.1f MIP", E_CUTS_GEV[2]), Form("Accepted (E< %.1f MIP); p_{MC} [GeV]; N",E_CUTS_GEV[2]), P_NB, P_MIN_GEV, P_MAX_GEV),
    };

    TH1D* hP_pass_LayerCut[3] = {
        new TH1D(Form("hP_pass_Layer< %d layers", LAYER_CUTS[0]), Form("Accepted (Layer < %d layers); p_{MC} [GeV]; N",LAYER_CUTS[0]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_Layer< %d layers", LAYER_CUTS[1]), Form("Accepted (Layer < %d layers); p_{MC} [GeV]; N",LAYER_CUTS[1]), P_NB, P_MIN_GEV, P_MAX_GEV),
        new TH1D(Form("hP_pass_Layer< %d layers", LAYER_CUTS[2]), Form("Accepted (Layer < %d layers); p_{MC} [GeV]; N",LAYER_CUTS[2]), P_NB, P_MIN_GEV, P_MAX_GEV),
    };

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

        bool m1_has_rec = false;
        bool m2_has_rec = false;

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

        if(inNHCal(v1.Eta())) {hP_all_mu->Fill(v1.P());}
        if(inNHCal(v2.Eta())){hP_all_mu->Fill(v2.P());}

        struct TaggedReco { edm4eic::ReconstructedParticle reco; int muTag; /*0=m1, 1=m2*/ };
        vector<TaggedReco> matchedRecos;

        auto find_associated_reco = [&](const edm4hep::MCParticle& mc, int muTag)->void {
            try {
                if (!assocCol.isValid() || assocCol.empty()) {return;}
                auto simIDs = assocCol.simID();
                auto recIDs = assocCol.recID();
                const uint32_t mc_idx = mc.getObjectID().index;
                for (size_t i=0; i<assocCol.size() && i<simIDs.size() && i<recIDs.size(); ++i) {
                    if (simIDs[i] == mc_idx) {
                        uint32_t ridx = recIDs[i];
                        if (!recParts.isValid() || ridx >= recParts.size()) {continue;}
                        auto reco = recParts.at(ridx);
                        if (reco.isAvailable()) {
                            matchedRecos.push_back({reco, muTag});
                        }
                    }
                }
            } catch (...) {}
        };

        if (!assocCol.isValid() || assocCol.empty()) continue;
        if (m1.isAvailable() && abs(m1.getPDG())==13) find_associated_reco(m1, 0);
        if (m2.isAvailable() && abs(m2.getPDG())==13) find_associated_reco(m2, 1);

        if (!recParts.isValid() || !projSegs.isValid() || !hcalRec.isValid() || recParts.empty() || projSegs.empty() || hcalRec.empty()) continue;
        
        set<double> uniqueCentersZ;
        map<double, double> layerData;
        for (const auto& hit : hcalRec) {    
            double z = hit.getPosition().z;
            double zBin = round(z);
            uint64_t cid = hit.getCellID();
            double zCenter_mm = 10 * getPlasticCenterZ_cm(decoder, cellid_converter, cid, slice_index, /*plastic_slice_value=*/3);
            if (!std::isnan(zCenter_mm)) {uniqueCentersZ.insert(zCenter_mm);}

            layerData[zBin] += hit.getEnergy(); 
            hZ_hits->Fill(z);                
            hE_z->Fill(z, hit.getEnergy()); 
            hE->Fill(hit.getEnergy());
        }

        vector<double> layerCentersZ(uniqueCentersZ.begin(), uniqueCentersZ.end());

        for (const auto& [zValue, sumEnergy] : layerData) {hEsum_z->Fill(zValue, sumEnergy); hEsum->Fill(sumEnergy);}

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

            size_t layerCentersZSize = layerCentersZ.size();
            vector<double> segMinDistance(layerCentersZSize, std::numeric_limits<double>::infinity());
            vector<double> segHitEnergy(layerCentersZSize, std::numeric_limits<double>::quiet_NaN());
            vector<double> thickness_cm(layerCentersZSize, std::numeric_limits<double>::quiet_NaN());
            vector<int> count_DrCuts(SIZE, 0);
            vector<int> count_ECuts(SIZE, 0);

            GeomState A{}, B{};
            bool haveA = false, haveB = false;

            auto points = ST.seg.getPoints();

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

            for (int i = 0; i < layerCentersZSize; ++i) {
                double best_dr_in_layer = 1e10;
                double best_E_in_layer;
                dd4hep::DDSegmentation::CellID best_cid_in_layer;

                double zc = layerCentersZ[i];
                auto [X, Y] = trackXYatZ(A, B, zc);

                for (const auto& hit : hcalRec) {
                    const auto& hp = hit.getPosition();

                    if (fabs(hp.z - zc) > 5.0 /*mm*/) continue;

                    const double dx = X - hp.x;
                    const double dy = Y - hp.y;
                    const double dr = sqrt(dx*dx + dy*dy);

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
                    thickness_cm[i]   = getPlasticThicknessCM(*det, decoder, best_cid_in_layer, slice_index, 3);
                }
                for(int j = 0; j < SIZE; ++j){
                    if (segMinDistance[i] < DR_CUTS_CM[j] * 10){++count_DrCuts[j];}
                    if (segHitEnergy[i] < (E_CUTS_GEV[j] * thickness_cm[i] * 100)){++count_ECuts[j];}
                }
                
            }

            if (ST.muTag == 0 && count_DrCuts[SIZE-1] > LAYER_THRESH) m1_has_rec = true;
            if (ST.muTag == 1 && count_DrCuts[SIZE-1] > LAYER_THRESH) m2_has_rec = true;

            for (int j = 0; j < SIZE; ++j){
                const double p = ST.muTag ? v2.P() : v1.P();
                if (count_DrCuts[j] > LAYER_THRESH) hP_pass_dr[j]->Fill(p);
                if (count_ECuts[j] > LAYER_THRESH) hP_pass_ECut[j]->Fill(p);
                if (count_DrCuts[SIZE-1] > LAYER_CUTS[j]) hP_pass_LayerCut[j]->Fill(p);
            }
        }   
        

        int nInNH_rec_local = 0;
        if (m1_has_rec) ++nInNH_rec_local;
        if (m2_has_rec) ++nInNH_rec_local;

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
    canvas_sim->cd(1); setPadMargins(); hEtaPt->Draw("COLZ");

    auto [mg_q2, leg_q2, log_q2] = makeEffMultiGraph(hQ2_all, hQ2_in, "Geom. acc. vs Q^{2}", "Q^{2} [GeV^{2}]", true);
    canvas_sim->cd(2); setPadMargins(); if (log_q2) gPad->SetLogx(); gPad->SetGrid(); mg_q2->Draw("A"); leg_q2->Draw();

    auto [mg_x, leg_x, log_x] = makeEffMultiGraph(hX_all, hX_in, "Geom. acc. vs x", "x", true);
    canvas_sim->cd(3); setPadMargins(); if (log_x) gPad->SetLogx(); gPad->SetGrid(); mg_x->Draw("A"); leg_x->Draw();

    auto [mg_y, leg_y, log_y] = makeEffMultiGraph(hY_all, hY_in, "Geom. acc. vs y", "y", false);
    canvas_sim->cd(4); setPadMargins(); gPad->SetGrid(); mg_y->Draw("A"); leg_y->Draw();

    auto [mg_w, leg_w, log_w] = makeEffMultiGraph(hW_all, hW_in, "Geom. acc. vs W", "W [GeV]", false);
    canvas_sim->cd(5); setPadMargins(); gPad->SetGrid(); mg_w->Draw("A"); leg_w->Draw();

    canvas_sim->cd(6); setPadMargins(); gPad->SetLogx(); gPad->SetGrid(); gPad->SetLogy(); hx_Q2->Draw("COLZ");
    canvas_sim->cd(7); setPadMargins(); hEta1_Eta2->Draw("COLZ"); 
    
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
    canvas_rec->cd(1); setPadMargins(); hEtaPt->Draw("COLZ");

    auto [mg_rec_q2, leg_rec_q2, log_rec_q2] = makeEffMultiGraph(hQ2_rec_all, hQ2_rec_in, "Geom. acc. vs Q^{2}", "Q^{2} [GeV^{2}]", true);
    canvas_rec->cd(2); setPadMargins(); if (log_rec_q2) gPad->SetLogx(); gPad->SetGrid(); mg_rec_q2->Draw("A"); leg_rec_q2->Draw();

    auto [mg_rec_x, leg_rec_x, log_rec_x] = makeEffMultiGraph(hX_rec_all, hX_rec_in, "Geom. acc. vs x", "x", true);
    canvas_rec->cd(3); setPadMargins(); if (log_rec_x) gPad->SetLogx(); gPad->SetGrid(); mg_rec_x->Draw("A"); leg_rec_x->Draw();

    auto [mg_rec_y, leg_rec_y, log_rec_y] = makeEffMultiGraph(hY_rec_all, hY_rec_in, "Geom. acc. vs y", "y", false);
    canvas_rec->cd(4); setPadMargins(); gPad->SetGrid(); mg_rec_y->Draw("A"); leg_rec_y->Draw();

    auto [mg_rec_w, leg_rec_w, log_rec_w] = makeEffMultiGraph(hW_rec_all, hW_rec_in, "Geom. acc. vs W", "W [GeV]", false);
    canvas_rec->cd(5); setPadMargins(); gPad->SetGrid(); mg_rec_w->Draw("A"); leg_rec_w->Draw();

    canvas_rec->cd(6); setPadMargins(); gPad->SetLogx(); gPad->SetGrid(); gPad->SetLogy(); hx_Q2->Draw("COLZ");
    canvas_rec->cd(7); setPadMargins(); hEta1_Eta2->Draw("COLZ");

    canvas_rec->SaveAs(addPrefixAfterSlash(outname_png, "rec_").c_str());
    canvas_rec->SaveAs(addPrefixAfterSlash(outname_pdf, "rec_").c_str());

    double letterSize = 0.045;
    TCanvas* canvas1 = new TCanvas("canvas1", "Muon analysis extra 1", 1800, 1200);
    canvas1->Divide(2,2);                                                                 
    canvas1->cd(1); setPadMargins(); gPad->SetGrid(); hZ_proj->SetTitleSize(letterSize, "XYZ"); hZ_proj->Draw();                                     
    canvas1->cd(2); setPadMargins(); gPad->SetGrid(); hZ_hits->SetTitleSize(letterSize, "XYZ"); hZ_hits->Draw();
    canvas1->Update();
    canvas1->SaveAs(addPrefixAfterSlash(outname_png, "extra1_").c_str());              
    canvas1->SaveAs(addPrefixAfterSlash(outname_pdf, "extra1_").c_str()); 

    TCanvas* canvas2 = new TCanvas("canvas2", "Muon analysis extra 2", 1800, 1200);
    canvas2->Divide(2,2); 
    canvas2->cd(1); setPadMargins(); gPad->SetGrid(); hDxDyZ_layer->SetTitleSize(letterSize, "XYZ"); hDxDyZ_layer->Draw("COLZ");     
    canvas2->cd(2); setPadMargins(); gPad->SetGrid(); hDxDy_all->SetTitleSize(letterSize, "XYZ"); hDxDy_all->Draw("COLZ"); 
    canvas2->cd(3); setPadMargins(); gPad->SetGrid(); hDrZ_layer->SetTitleSize(letterSize, "XYZ"); hDrZ_layer->Draw("COLZ");                             
    canvas2->cd(4); setPadMargins(); gPad->SetGrid(); hDr_all->SetTitleSize(letterSize, "XYZ"); hDr_all->Draw("COLZ");
    canvas2->SaveAs(addPrefixAfterSlash(outname_png, "extra2_").c_str());              
    canvas2->SaveAs(addPrefixAfterSlash(outname_pdf, "extra2_").c_str()); 

    TCanvas* canvas3 = new TCanvas("canvas3", "Muon analysis extra 3", 1800, 1200);
    canvas3->Divide(2,2);                        
    canvas3->cd(1); setPadMargins(); gPad->SetGrid(); hP_all_mu->SetLineColor(kBlack); hP_all_mu->SetTitleSize(letterSize, "XYZ"); hP_all_mu->Draw();  
    canvas3->cd(2); setPadMargins(); gPad->SetGrid();                                                       
    TH1D* hEff_dr[3]; 
    for (int idr=0; idr<3; ++idr) { 
        hEff_dr[idr] = (TH1D*)hP_pass_dr[idr]->Clone( 
            Form("hEff_dr%.1fcm", DR_CUTS_CM[idr]) 
        ); 
        hEff_dr[idr]->SetDirectory(nullptr);
        hEff_dr[idr]->SetTitle("Matching efficiency vs p_{MC} (Layers > 7); p_{MC} [GeV]; efficiency"); 
        hEff_dr[idr]->Divide(hP_all_mu); 
    } 

    hEff_dr[0]->SetLineColor(kBlue);     hEff_dr[0]->SetMinimum(0.0); hEff_dr[0]->SetMaximum(1.4); 
    hEff_dr[1]->SetLineColor(kRed);   
    hEff_dr[2]->SetLineColor(kGreen+2); 

    hEff_dr[0]->SetTitleSize(letterSize, "XYZ"); hEff_dr[0]->Draw("HIST");            
    hEff_dr[1]->SetTitleSize(letterSize, "XYZ"); hEff_dr[1]->Draw("HIST SAME");       
    hEff_dr[2]->SetTitleSize(letterSize, "XYZ"); hEff_dr[2]->Draw("HIST SAME");       
    { auto leg_dr = new TLegend(0.60,0.70,0.88,0.88);
        for(int idr=0; idr<3; ++idr)
        {leg_dr->AddEntry(hEff_dr[idr],Form("dr < %.1f cm", DR_CUTS_CM[idr]),"l");} 
        leg_dr->Draw(); 
    }  
                                                                   
    canvas3->cd(3); setPadMargins(); gPad->SetGrid();                                                        
    TH1D* hEff_E[3]; 
    for (int ie=0; ie<3; ++ie) { 
        hEff_E[ie] = (TH1D*)hP_pass_ECut[ie]->Clone( 
            Form("hEff_E%.1fcm", E_CUTS_GEV[ie]) 
        ); 
        hEff_E[ie]->SetDirectory(nullptr);
        hEff_E[ie]->SetTitle("Matching efficiency vs p_{MC} (Layers > 7, dr < 13cm); p_{MC} [GeV]; efficiency"); 
        hEff_E[ie]->Divide(hP_all_mu); 
    } 

    hEff_E[0]->SetLineColor(kBlue);     hEff_E[0]->SetMinimum(0.0); hEff_E[0]->SetMaximum(1.4); 
    hEff_E[1]->SetLineColor(kRed);   
    hEff_E[2]->SetLineColor(kGreen+2); 

    hEff_E[0]->SetTitleSize(letterSize, "XYZ"); hEff_E[0]->Draw("HIST");            
    hEff_E[0]->SetTitleSize(letterSize, "XYZ"); hEff_E[1]->Draw("HIST SAME");       
    hEff_E[0]->SetTitleSize(letterSize, "XYZ"); hEff_E[2]->Draw("HIST SAME");       
    { auto leg_E = new TLegend(0.60,0.70,0.88,0.88);
        for(int ie=0; ie<3; ++ie)
        {leg_E->AddEntry(hEff_E[ie],Form("E < %.1f MIP", E_CUTS_GEV[ie]),"l");}  
        leg_E->Draw(); 
    }
    
    canvas3->cd(4); setPadMargins(); gPad->SetGrid();
    TH1D* hEff_Layer[3]; 
    for (int il=0; il<3; ++il) { 
        hEff_Layer[il] = (TH1D*)hP_pass_LayerCut[il]->Clone( 
            Form("hEff_Layer%d", LAYER_CUTS[il]) 
        ); 
        hEff_Layer[il]->SetDirectory(nullptr);
        hEff_Layer[il]->SetTitle("Matching efficiency vs p_{MC} (dr < 13cm); p_{MC} [GeV]; efficiency"); 
        hEff_Layer[il]->Divide(hP_all_mu); 
    } 

    hEff_Layer[0]->SetLineColor(kBlue); hEff_Layer[0]->SetMinimum(0.0); hEff_Layer[0]->SetMaximum(1.4); 
    hEff_Layer[1]->SetLineColor(kRed);   
    hEff_Layer[2]->SetLineColor(kGreen+2); 

    hEff_Layer[0]->SetTitleSize(letterSize, "XYZ"); hEff_Layer[0]->Draw("HIST");            
    hEff_Layer[0]->SetTitleSize(letterSize, "XYZ"); hEff_Layer[1]->Draw("HIST SAME");       
    hEff_Layer[0]->SetTitleSize(letterSize, "XYZ"); hEff_Layer[2]->Draw("HIST SAME");       
    { auto leg_Layer = new TLegend(0.60,0.70,0.88,0.88);
        for(int il=0; il<3; ++il)
        {leg_Layer->AddEntry(hEff_Layer[il],Form("L > %d Layers", LAYER_CUTS[il]),"l");}  
        leg_Layer->Draw(); 
    }  
    canvas3->SaveAs(addPrefixAfterSlash(outname_png, "extra3_").c_str());              
    canvas3->SaveAs(addPrefixAfterSlash(outname_pdf, "extra3_").c_str()); 
    
    TCanvas* canvas4 = new TCanvas("canvas4", "Muon analysis extra 4", 1800, 1200);
    canvas4->Divide(2,2); 

    canvas4->cd(1); setPadMargins(); gPad->SetGrid(); gPad->SetLogy(); hE_z->GetYaxis()->SetMoreLogLabels(); hE_z->GetYaxis()->SetTitle("E [GeV] (log)"); 
                    hE_z->SetTitleSize(letterSize, "XYZ"); hE_z->Draw("COLZ");
    TProfile* pE_z = hE_z->ProfileX("pE_z"); pE_z->SetLineWidth(3); pE_z->SetLineColor(kRed); pE_z->SetMarkerColor(kRed);
                    pE_z->SetDirectory(nullptr); pE_z->Draw("SAME");
    // const double Ecut_GeV = MIP_ENERGY_GEV * t_cm * 100;
    // TLine* cut = new TLine(hE_z->GetXaxis()->GetXmin(), Ecut_GeV, hE_z->GetXaxis()->GetXmax(), Ecut_GeV);
    // cut->SetLineColor(kRed+1); cut->SetLineStyle(2); cut->SetLineWidth(3);
    // cut->Draw("SAME");                              
    //SetRangeUser
    canvas4->cd(2); setPadMargins(); gPad->SetGrid(); gPad->SetLogy(); hEsum_z->GetYaxis()->SetMoreLogLabels(); hEsum_z->GetYaxis()->SetTitle("E_{sum} [GeV] (log)"); 
                    hEsum_z->SetTitleSize(letterSize, "XYZ"); hEsum_z->Draw("COLZ");
    TProfile* pEsum_z = hEsum_z->ProfileX("pEsum_z"); pEsum_z->SetLineWidth(3); pEsum_z->SetLineColor(kRed); pEsum_z->SetMarkerColor(kRed);
                    pEsum_z->SetDirectory(nullptr); pEsum_z->Draw("SAME"); 
    canvas4->cd(3); setPadMargins(); gPad->SetGrid(); hE->Draw();
    canvas4->cd(4); setPadMargins(); gPad->SetGrid(); hEsum->Draw();
    
    canvas4->Update();
    canvas4->SaveAs(addPrefixAfterSlash(outname_png, "extra4_").c_str());              
    canvas4->SaveAs(addPrefixAfterSlash(outname_pdf, "extra4_").c_str()); 

    delete canvas_sim;
    delete canvas_rec;
    delete canvas1;
    delete canvas2;
    delete canvas3;
    delete canvas4;

    return 0;
}