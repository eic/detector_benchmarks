#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
#include "TLorentzVector.h"  
#include "TStyle.h"    

#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"

#include <podio/Frame.h>
#include <podio/CollectionBase.h>
#include "podio/ROOTReader.h"
#include "podio/CollectionIDTable.h"
#include "podio/ObjectID.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticle.h"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHit.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHit.h"

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/Cluster.h"

#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHit.h"

#include <edm4eic/vector_utils_legacy.h>

using namespace std;
using namespace edm4hep;

dd4hep::Detector* det = NULL;

constexpr double ETA_MIN = -4.16, ETA_MAX = -1.16;

inline bool inNHCal(double eta) {return (eta >= ETA_MIN && eta <= ETA_MAX);}

inline string addPrefixAfterSlash(const string& path,
                                   const string& prefix) {
    const auto slash = path.find_last_of('/');
    const auto dot = path.find_last_of('.');
    
    string extension = (dot != string::npos) ? path.substr(dot) : "";
    
    if (slash == string::npos)
        return prefix + extension;
    
    return path.substr(0, slash + 1) + prefix + extension;
}

int sampling_fraction_analysis(const string &filename, string outname_pdf, string outname_png, TString compact_file) 
{

    podio::ROOTReader *reader = new podio::ROOTReader();
    reader->openFile(filename);
    unsigned nEvents = reader->getEntries("events");
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

    constexpr int NBINS = 80; 

    constexpr double E_MIN_GEV = 0.0;
    constexpr double E_MAX_GEV = 12.0;
    constexpr double SAMP_F_MIN = 0.0;
    constexpr double SAMP_F_MAX = 1.0;
    constexpr double SAMP_F_LOW = 0.2;

    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.20);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);

    TH2D *h_sampF_e = new TH2D("h_sampF_e", "nHCal sampling fraction vs. energy (e-); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{hit}]; counts", 
                                            NBINS, E_MIN_GEV, E_MAX_GEV, NBINS, SAMP_F_MIN, SAMP_F_MAX);
   
    TH2D *h_sampF_n = new TH2D("h_sampF_n", "nHCal sampling fraction vs. energy (neutron); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{hit}]; counts", 
                                            NBINS, E_MIN_GEV, E_MAX_GEV, NBINS, SAMP_F_MIN, SAMP_F_MAX);

    TH2D *h_sampF_pi = new TH2D("h_sampF_pi", "nHCal sampling fraction vs. energy (#pi-); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{hit}]; counts", 
                                            NBINS, E_MIN_GEV, E_MAX_GEV, NBINS, SAMP_F_MIN, SAMP_F_MAX);                                           
    
    TH2D *h_sampF_e_Ekin = new TH2D("h_sampF_e_Ekin", "nHCal sampling fraction vs. energy kin (e-); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{kin}]; counts", 
                                            NBINS, E_MIN_GEV, E_MAX_GEV, NBINS, SAMP_F_MIN, SAMP_F_LOW);
   
    TH2D *h_sampF_n_Ekin = new TH2D("h_sampF_n_Ekin", "nHCal sampling fraction vs. energy kin (neutron); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{kin}]; counts", 
                                            NBINS, E_MIN_GEV, E_MAX_GEV, NBINS, SAMP_F_MIN, SAMP_F_LOW);

    TH2D *h_sampF_pi_Ekin = new TH2D("h_sampF_pi_Ekin", "nHCal sampling fraction vs. energy kin (#pi-); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{kin}]; counts", 
                                            NBINS, E_MIN_GEV, E_MAX_GEV, NBINS, SAMP_F_MIN, SAMP_F_LOW);

    for (unsigned ev = 0; ev < nEvents; ev++) 
    {
        double hit_Esum = 0;
        double hit_scint_Esum = 0;
        double singlePart_Ekin = 0;

        auto frameData = reader->readNextEntry(podio::Category::Event);
        if (!frameData) 
        {
            cerr << "Invalid FrameData at event " << ev << endl;
            continue;
        }

        podio::Frame frame(std::move(frameData));

        auto& MCParticles_coll  = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        auto& SimCalorimeterHit_coll = frame.get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");

        if (!SimCalorimeterHit_coll.isValid())    
        {
            cerr << "HcalEndcapNHits collection is not valid!" << endl;
        }
        if (!MCParticles_coll.isValid())    
        {
            cerr << "MCParticles collection is not valid!" << endl;
        }

        edm4hep::MCParticle mcpart =  MCParticles_coll.at(0);
        TLorentzVector v(mcpart.getMomentum().x, mcpart.getMomentum().y, mcpart.getMomentum().z, mcpart.getEnergy());
        if(!inNHCal(v.Eta())) continue;
        singlePart_Ekin = mcpart.getEnergy() - mcpart.getMass();
        int pdg = mcpart.getPDG();

        for (const auto& hit : SimCalorimeterHit_coll) 
        {
                const int hit_slice = decoder->get(hit.getCellID(), slice_index);
                hit_Esum += hit.getEnergy();
                if(hit_slice == 3) hit_scint_Esum += hit.getEnergy();
        }
        
        if (pdg == 11)               // e-
        {
            h_sampF_e->Fill(singlePart_Ekin, hit_scint_Esum/hit_Esum);
            h_sampF_e_Ekin->Fill(singlePart_Ekin, hit_scint_Esum/singlePart_Ekin);
        }
        else if (pdg == -211)        // pi-
        {
            h_sampF_pi->Fill(singlePart_Ekin, hit_scint_Esum/hit_Esum);
            h_sampF_pi_Ekin->Fill(singlePart_Ekin, hit_scint_Esum/singlePart_Ekin);
        }
        else if (pdg == 2112)        // neutron
        {   
            h_sampF_n->Fill(singlePart_Ekin, hit_scint_Esum/hit_Esum);
            h_sampF_n_Ekin->Fill(singlePart_Ekin, hit_scint_Esum/singlePart_Ekin);
        }   
    }

    delete reader; 

    h_sampF_e->Sumw2();
    h_sampF_e_Ekin->Sumw2();
    h_sampF_pi->Sumw2();
    h_sampF_pi_Ekin->Sumw2();
    h_sampF_n->Sumw2();
    h_sampF_n_Ekin->Sumw2();

    TProfile* p_sampF_e = h_sampF_e->ProfileX();                
    p_sampF_e ->SetTitle("nHCal sampling fraction vs. energy (e-); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{hit}]");
    TProfile* p_sampF_e_Ekin = h_sampF_e_Ekin->ProfileX();      
    p_sampF_e_Ekin ->SetTitle("nHCal sampling fraction vs. energy kin (e-); E_{kin} [GeV]; f");
    TProfile* p_sampF_pi = h_sampF_pi->ProfileX();              
    p_sampF_pi ->SetTitle("nHCal sampling fraction vs. energy (#pi-); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{hit}]");
    TProfile* p_sampF_pi_Ekin = h_sampF_pi_Ekin->ProfileX();    
    p_sampF_pi_Ekin ->SetTitle("nHCal sampling fraction vs. energy kin (#pi-); E_{kin} [GeV]; f");
    TProfile* p_sampF_n = h_sampF_n->ProfileX();                
    p_sampF_n ->SetTitle("nHCal sampling fraction vs. energy (neutron); E_{kin} [GeV]; sampling fraction from hits [E_{scint}/E_{hit}]");
    TProfile* p_sampF_n_Ekin = h_sampF_n_Ekin->ProfileX();      
    p_sampF_n_Ekin ->SetTitle("nHCal sampling fraction vs. energy kin (neutron); E_{kin} [GeV]; f");

    TH1D *h_e = p_sampF_e->ProjectionX("h_e");
    TH1D *h_pi = p_sampF_pi->ProjectionX("h_pi");

    TH1D *h_e_over_pi = (TH1D*)h_e->Clone("h_e_over_pi");
    h_e_over_pi->SetTitle("e/h ratio;E_{kin} [GeV];e/h");
    h_e_over_pi->Divide(h_e, h_pi, 1, 1); 


    TCanvas *c_h = new TCanvas("canvas_h", "c_h", 1600, 800);
    c_h->Divide(2,2);
    c_h->cd(1);
    h_sampF_e->Draw("COLZ");

    c_h->cd(2);
    h_sampF_pi->Draw("COLZ");

    c_h->cd(3);
    h_sampF_n->Draw("COLZ");

    c_h->SaveAs(outname_pdf.c_str());
    c_h->SaveAs(outname_png.c_str());

    TCanvas *c_p = new TCanvas("canvas_p", "c_p", 1600, 800);
    c_p->Divide(2,2);
    c_p->cd(1);
    p_sampF_e->SetLineWidth(3); p_sampF_e->SetLineColor(kRed); p_sampF_e->SetMarkerColor(kRed);
    p_sampF_e->Draw();
    c_p->cd(2);
    p_sampF_pi->SetLineWidth(3); p_sampF_pi->SetLineColor(kRed); p_sampF_pi->SetMarkerColor(kRed);
    p_sampF_pi->Draw();
    c_p->cd(3);
    p_sampF_n->SetLineWidth(3); p_sampF_n->SetLineColor(kRed); p_sampF_n->SetMarkerColor(kRed);
    p_sampF_n->Draw();
    c_p->cd(4);
    h_e_over_pi->Draw();

    c_p->SaveAs(addPrefixAfterSlash(outname_png, "prof_sampf_vs_Ehit").c_str());
    c_p->SaveAs(addPrefixAfterSlash(outname_pdf, "prof_sampf_vs_Ehit").c_str());

    TH1D *h_e_Ekin = p_sampF_e_Ekin->ProjectionX("h_e_Ekin");
    TH1D *h_pi_Ekin = p_sampF_pi_Ekin->ProjectionX("h_pi_Ekin");

    TH1D *h_e_over_pi_Ekin = (TH1D*)h_e_Ekin->Clone("h_e_over_pi_Ekin");
    h_e_over_pi_Ekin->SetTitle("e/#pi ratio;E_{kin} [GeV];e/#pi");
    h_e_over_pi_Ekin->Divide(h_e_Ekin, h_pi_Ekin, 1, 1); 

    TCanvas *c_h_Ekin = new TCanvas("canvas_h_Ekin", "c_h_Ekin", 1600, 800);
    c_h_Ekin->Divide(2,2);
    c_h_Ekin->cd(1);
    h_sampF_e_Ekin->Draw("COLZ");

    c_h_Ekin->cd(2);
    h_sampF_pi_Ekin->Draw("COLZ");

    c_h_Ekin->cd(3);
    h_sampF_n_Ekin->Draw("COLZ");

    c_h_Ekin->SaveAs(addPrefixAfterSlash(outname_png, "hist_sampf_vs_Ekin").c_str());
    c_h_Ekin->SaveAs(addPrefixAfterSlash(outname_pdf, "hist_sampf_vs_Ekin").c_str());

    TCanvas *c_p_Ekin = new TCanvas("canvas_p_Ekin", "c_p_Ekin", 1600, 800);
    c_p_Ekin->Divide(2,2);
    c_p_Ekin->cd(1);
    p_sampF_e_Ekin->SetLineWidth(3); p_sampF_e_Ekin->SetLineColor(kRed); p_sampF_e_Ekin->SetMarkerColor(kRed);
    p_sampF_e_Ekin->Draw();
    c_p_Ekin->cd(2);
    p_sampF_pi_Ekin->SetLineWidth(3); p_sampF_pi_Ekin->SetLineColor(kRed); p_sampF_pi_Ekin->SetMarkerColor(kRed);
    p_sampF_pi_Ekin->Draw();
    c_p_Ekin->cd(3);
    p_sampF_n_Ekin->SetLineWidth(3); p_sampF_n_Ekin->SetLineColor(kRed); p_sampF_n_Ekin->SetMarkerColor(kRed);
    p_sampF_n_Ekin->Draw();
    c_p_Ekin->cd(4);
    h_e_over_pi_Ekin->Draw();

    c_p_Ekin->SaveAs(addPrefixAfterSlash(outname_png, "prof_sampf_vs_Ekin").c_str());
    c_p_Ekin->SaveAs(addPrefixAfterSlash(outname_pdf, "prof_sampf_vs_Ekin").c_str());    

    return 0;
}