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

#include "TStyle.h"

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

// #include "FileList.h"
// #include "EICutil.h"
// #include "BasicUtil.h"

// #pragma link C++ class vector<edm4hep::MCParticleData>+;
// #pragma link C++ class vector<podio::ObjectID>+;
// #pragma link C++ class vector<edm4hep::SimCalorimeterHitData>+;

using namespace std;
using namespace ROOT;
using namespace TMath;
using namespace edm4hep;

int basic_distribution_analysis(const string &filename, string outname_pdf, string outname_png) 
{
    podio::ROOTReader *reader = new podio::ROOTReader();
    reader->openFile(filename);
    unsigned nEvents = reader->getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    vector<double> xPosition;
    vector<double> yPosition;
    vector<double> zPosition;
    vector<double> rPosition;
    vector<double> energy;
    vector<double> totalEventEnergy;
    vector<double> totalEventHits;
    vector<double> zLayerValues;
    vector<int> nHitsPerLayer;
    vector<double> nEnergyPerLayer;
    
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
        if (!hits.isValid())    
        {
            cerr << "HcalEndcapNHits collection is not valid!" << endl;
        }

        map<double, pair<int, double>> layerData;
        double totalEnergy = 0;
        int totalHits = 0;

        for (const auto& hit : hits) 
        {
            totalHits++;
            energy.push_back(hit.getEnergy());
            totalEnergy += hit.getEnergy();

            auto pos = hit.getPosition();
            xPosition.push_back(pos.x);
            yPosition.push_back(pos.y);
            zPosition.push_back(pos.z);
            rPosition.push_back(sqrt(pos.x * pos.x + pos.y * pos.y));

            double zBin = round(pos.z);
            layerData[zBin].first++;         
            layerData[zBin].second += hit.getEnergy();
        }

        totalEventEnergy.push_back(totalEnergy);
        totalEventHits.push_back(totalHits);

        for (const auto& [zValue, stats] : layerData)
        {
            zLayerValues.push_back(zValue);
            nHitsPerLayer.push_back(stats.first);
            nEnergyPerLayer.push_back(stats.second);
        }

    }

    auto minmax_xPosition        = minmax_element(xPosition.begin(), xPosition.end());
    auto minmax_yPosition        = minmax_element(yPosition.begin(), yPosition.end());
    auto minmax_zPosition        = minmax_element(zPosition.begin(), zPosition.end());
    auto minmax_rPosition        = minmax_element(rPosition.begin(), rPosition.end());
    auto minmax_totalEventEnergy = minmax_element(totalEventEnergy.begin(), totalEventEnergy.end());
    auto minmax_totalEventHits   = minmax_element(totalEventHits.begin(), totalEventHits.end());
    auto minmax_energy           = minmax_element(energy.begin(), energy.end());
    auto minmax_zLayerValues     = minmax_element(zLayerValues.begin(), zLayerValues.end());
    auto minmax_nHitsPerLayer    = minmax_element(nHitsPerLayer.begin(), nHitsPerLayer.end());
    auto minmax_nEnergyPerLayer  = minmax_element(nEnergyPerLayer.begin(), nEnergyPerLayer.end());

    int nbins = nEvents/1000;

    TH1D *h_energyTotal = new TH1D("h_energyTotal", "Total Energy per Event; Energy per Event", 
                                    nbins/2, *minmax_totalEventEnergy.first, *minmax_totalEventEnergy.second);
    TH2D *h_layerEnergy = new TH2D("h_layerEnergy", "Energy in Layers (Z); Z; Energy", 
                                    nbins/5, *minmax_zLayerValues.first, *minmax_zLayerValues.second, 
                                    nbins/4, *minmax_nEnergyPerLayer.first, *minmax_nEnergyPerLayer.second);
    TProfile *p_layerEnergy = new TProfile("p_layerEnergy", "Energy in Layers (Z); Z; Mean Energy",
                                    nbins/5, *minmax_zLayerValues.first, *minmax_zLayerValues.second);
    TH1D *h_hitCount    = new TH1D("h_hitCount", "Number of Hits per Event; Hits per Event", 
                                    nbins/2, *minmax_totalEventHits.first, *minmax_totalEventHits.second);
    TH2D *h_layerHits   = new TH2D("h_layerHits", "Number of Hits in Layers (Z); Layer(Z); Hits", 
                                    nbins/5, *minmax_zLayerValues.first, *minmax_zLayerValues.second, 
                                    *minmax_nHitsPerLayer.second, *minmax_nHitsPerLayer.first, *minmax_nHitsPerLayer.second);
    TProfile *p_layerHits = new TProfile("p_layerHits", "Number of Hits in Layers (Z); Z; Mean Hits",
                                    nbins/5, *minmax_zLayerValues.first, *minmax_zLayerValues.second);
    TH2D *h_XYPos       = new TH2D("h_XYPos", "Hits position X,Y;X [mm];Y [mm]", 
                                    nbins, *minmax_xPosition.first, *minmax_xPosition.second, 
                                    nbins, *minmax_yPosition.first, *minmax_yPosition.second);
    TH2D *h_ZRPos       = new TH2D("h_ZRPos", "Hits position Z,R;Z [mm];R [mm]", 
                                    nbins/5, *minmax_zPosition.first, *minmax_zPosition.second, 
                                    nbins, *minmax_rPosition.first, *minmax_rPosition.second);
    TH2D *h_XYEnergy    = new TH2D("h_XYEnergy", "Hits energy X,Y;X [mm];Y [mm]", 
                                    nbins, *minmax_xPosition.first, *minmax_xPosition.second, 
                                    nbins, *minmax_yPosition.first, *minmax_yPosition.second);

    for(int i = 0; i < xPosition.size(); i++)
    {
        
        h_XYPos->Fill(xPosition[i], yPosition[i]);
        h_ZRPos->Fill(zPosition[i], rPosition[i]);
        h_XYEnergy->Fill(xPosition[i], yPosition[i], energy[i]);
    }

    for (int i = 0; i < zLayerValues.size(); i++)
    {
        h_layerHits->Fill(zLayerValues[i], nHitsPerLayer[i]);
        p_layerHits->Fill(zLayerValues[i], nHitsPerLayer[i]);

        h_layerEnergy->Fill(zLayerValues[i], nEnergyPerLayer[i]);
        p_layerEnergy->Fill(zLayerValues[i], nEnergyPerLayer[i]);
    }

    for(int i = 0; i < nEvents; i++)
    {
        h_energyTotal->Fill(totalEventEnergy[i]);
        h_hitCount->Fill(totalEventHits[i]);
    }

    gStyle->SetPadLeftMargin(0.13);

    TCanvas *canvas = new TCanvas("canvas", "canvas", 1600, 800);
    canvas->Divide(4,2);
    canvas->cd(1);
    h_energyTotal->Draw();
    canvas->cd(2);
    h_layerEnergy->Draw("COLZ");
    p_layerEnergy->Draw("SAME");
    canvas->cd(3);
    h_hitCount->Draw();
    canvas->cd(4);
    h_layerHits->Draw("COLZ");
    p_layerHits->Draw("SAME");
    canvas->cd(5);
    h_XYPos->Draw("COLZ");
    canvas->cd(6);
    h_ZRPos->Draw("COLZ");
    canvas->cd(7);
    h_XYEnergy->Draw("COLZ");

    canvas->SaveAs(outname_pdf.c_str());
    canvas->SaveAs(outname_png.c_str());

    return 0;
}
