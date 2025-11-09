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


#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"

#include <podio/Frame.h>
#include <podio/CollectionBase.h>
#include "podio/ROOTReader.h"
//#include <podio/ROOTFrameReader.h>
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


#include <edm4eic/vector_utils_legacy.h>
#include <edm4hep/Vector3f.h>
#include <root/TStyle.h>

using namespace std;
using namespace ROOT;
using namespace TMath;
using namespace edm4hep;

inline string addPrefixAfterSlash(const string& path, const string& prefix) {
    const auto slash = path.find_last_of('/');
    if (slash == string::npos)
        return prefix + path;
    return path.substr(0, slash + 1) + prefix + path.substr(slash + 1);
}

string changeExtension(const string& path, const string& new_ext)
{
    size_t pos = path.find_last_of('.');
    if (pos != string::npos)
        return path.substr(0, pos) + new_ext;
    return path + new_ext;
}

int basic_distribution_analysis(const string &filename, string outname_png, string outname_root) 
{
    constexpr double NBINS = 50;
    constexpr double MIN_X_MM = -3000, MAX_X_MM = 3000;
    constexpr double MIN_Y_MM = -3000, MAX_Y_MM = 3000;
    constexpr double MIN_Z_MM = -4500, MAX_Z_MM = -3600;
    constexpr double MIN_R_MM = 0, MAX_R_MM = 3000;
    constexpr double MIN_TOTAL_ENERGY_GEV = 0, MAX_TOTAL_ENERGY_GEV = 0.3;
    constexpr double MIN_HITS = 0, MAX_HITS = 300;
    constexpr double MIN_LAYER_HITS = 0, MAX_LAYER_HITS = 50;
    constexpr double MIN_LAYER_ENERGY_GEV = 0, MAX_LAYER_ENERGY_GEV = 0.2;
    constexpr double MIN_KINETIC_GEV = 0, MAX_KINETIC_GEV = 6;

    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);

    string outname_pdf = changeExtension(outname_png, ".pdf");

    podio::ROOTReader *reader = new podio::ROOTReader();
    reader->openFile(filename);
    unsigned nEvents = reader->getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    TH1D *h_energyTotal = new TH1D("h_energyTotal", "Total Energy per Event; Energy per Event [GeV]", 
                                    NBINS, MIN_TOTAL_ENERGY_GEV, MAX_TOTAL_ENERGY_GEV);
    
    TH2D *h_layerEnergy = new TH2D("h_layerEnergy", "Energy in Layers; Z [mm]; Energy [GeV]", 
                                    NBINS, MIN_Z_MM, MAX_Z_MM, NBINS, MIN_LAYER_ENERGY_GEV, MAX_LAYER_ENERGY_GEV);
    
    TProfile *p_layerEnergy = new TProfile("p_layerEnergy", "Energy in Layers; Z [mm]; Mean Energy [GeV]",
                                            NBINS, MIN_Z_MM, MAX_Z_MM);
    
    TH1D *h_hitCount = new TH1D("h_hitCount", "Hits per Event; Hits per Event; N", 
                                NBINS, MIN_HITS, MAX_HITS);
    
    TH2D *h_layerHits = new TH2D("h_layerHits", "Hits in Layers; Z [mm]; Hits", 
                                  NBINS, MIN_Z_MM, MAX_Z_MM, NBINS, MIN_LAYER_HITS, MAX_LAYER_HITS);
    
    TProfile *p_layerHits = new TProfile("p_layerHits", "Hits in Layers; Z [mm]; Mean Hits",
                                          NBINS, MIN_Z_MM, MAX_Z_MM);
    
    TH2D *h_XYPos = new TH2D("h_XYPos", "Hits position;X [mm];Y [mm]", 
                              NBINS, MIN_X_MM, MAX_X_MM, NBINS, MIN_Y_MM, MAX_Y_MM);
    
    TH2D *h_ZRPos = new TH2D("h_ZRPos", "Hits position;Z [mm];R [mm]", 
                              NBINS, MIN_Z_MM, MAX_Z_MM, NBINS, MIN_R_MM, MAX_R_MM);
    
    TH2D *h_XYEnergy = new TH2D("h_XYEnergy", "Hits energy;X [mm];Y [mm]", 
                                 NBINS, MIN_X_MM, MAX_X_MM, NBINS, MIN_Y_MM, MAX_Y_MM);

    TH2D *h_energyRes = new TH2D("h_energyRes", "Kinetic energy vs sum hit energy;Ekin [GeV]; Ehit sum [GeV]", 
                                  NBINS, MIN_KINETIC_GEV, MAX_KINETIC_GEV,
                                  NBINS, MIN_TOTAL_ENERGY_GEV, MAX_TOTAL_ENERGY_GEV);
    
    TProfile *p_energyRes = new TProfile("p_energyRes", "Kinetic energy vs sum hit energy",
                                          NBINS, MIN_KINETIC_GEV, MAX_KINETIC_GEV, "s");
    
    for (unsigned ev = 0; ev < nEvents; ev++) 
    {
        auto frameData = reader->readNextEntry(podio::Category::Event);
        if (!frameData) 
        {
            cerr << "Invalid FrameData at event " << ev << endl;
            continue;
        }

        podio::Frame frame(std::move(frameData));

        auto& hits = frame.get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
        auto& mcCol = frame.get<edm4hep::MCParticleCollection>("MCParticles");
        if (!hits.isValid() && !mcCol.isValid())    
        {
            cerr << "HcalEndcapNHits or MCParticles collection is not valid!" << endl;
        }

        double Ekin = mcCol[0].getEnergy() - mcCol[0].getMass();

        map<double, pair<int, double>> layerData;
        double totalEnergy = 0;
        int totalHits = 0;

        for (const auto& hit : hits) 
        {
            totalHits++;
            totalEnergy += hit.getEnergy();

            auto pos = hit.getPosition();
            double r = sqrt(pos.x * pos.x + pos.y * pos.y);
            
            h_XYPos->Fill(pos.x, pos.y);
            h_ZRPos->Fill(pos.z, r);
            h_XYEnergy->Fill(pos.x, pos.y, hit.getEnergy());

            double zBin = round(pos.z);
            layerData[zBin].first++;         
            layerData[zBin].second += hit.getEnergy();
        }

        h_energyTotal->Fill(totalEnergy);
        h_hitCount->Fill(totalHits);
        h_energyRes->Fill(Ekin, totalEnergy);
        p_energyRes->Fill(Ekin, totalEnergy);

        for (const auto& [zValue, stats] : layerData)
        {
            h_layerHits->Fill(zValue, stats.first);
            p_layerHits->Fill(zValue, stats.first);
            h_layerEnergy->Fill(zValue, stats.second);
            p_layerEnergy->Fill(zValue, stats.second);
        }
    }
    
    TFile *outFile = new TFile(outname_root.c_str(), "RECREATE");
    h_energyRes->Write();
    p_energyRes->Write();
    outFile->Close();
    delete outFile;

    TCanvas *c_evLayers = new TCanvas("c_evLayers", "c_evLayers", 1600, 800);
    c_evLayers->Divide(2,2);
    c_evLayers->cd(1);
    h_energyTotal->Draw();
    c_evLayers->cd(2);
    h_layerEnergy->Draw("COLZ");
    p_layerEnergy->SetLineWidth(3); p_layerEnergy->SetLineColor(kRed); p_layerEnergy->SetMarkerColor(kRed);
    p_layerEnergy->Draw("SAME");
    c_evLayers->cd(3);
    h_hitCount->Draw();
    c_evLayers->cd(4);
    h_layerHits->Draw("COLZ");
    p_layerHits->SetLineWidth(3); p_layerHits->SetLineColor(kRed); p_layerHits->SetMarkerColor(kRed);
    p_layerHits->Draw("SAME");

    c_evLayers->SaveAs(outname_png.c_str());
    c_evLayers->SaveAs(outname_pdf.c_str());

    TCanvas *c_hit_posE = new TCanvas("c_hit_posE", "c_hit_posE", 1600, 800);
    c_hit_posE->Divide(2,2);
    c_hit_posE->cd(1);
    h_XYPos->Draw("COLZ");
    c_hit_posE->cd(2);
    h_ZRPos->Draw("COLZ");
    c_hit_posE->cd(3);
    h_XYEnergy->Draw("COLZ");
    c_hit_posE->cd(4);
    
    c_hit_posE->SaveAs(addPrefixAfterSlash(outname_png, "hit_posE_").c_str());
    c_hit_posE->SaveAs(addPrefixAfterSlash(outname_pdf, "hit_posE_").c_str());

    delete reader;
    delete c_evLayers;
    delete c_hit_posE;

    return 0;
}