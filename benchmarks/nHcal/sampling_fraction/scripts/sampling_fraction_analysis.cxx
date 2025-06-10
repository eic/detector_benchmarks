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

dd4hep::Detector* det = NULL;

int sampling_fraction_analysis() 
{
    const char *inputFile = "sim_epic_backward_hcal_only_E5.0GeV_combined_10files.edm4hep.root";   
    const char *outputFile = "output.pdf";

    podio::ROOTReader *reader = new podio::ROOTReader();
    reader->openFile(inputFile);
    unsigned nEvents = reader->getEntries("events");
    cout << "Number of events: " << nEvents << endl;

    TString compact_file = "/opt/software/linux-debian12-x86_64_v2/gcc-12.2.0/epic-git.4fb92ec2884291bffddbf330c117cefc9958e9af_main-3kkvnmpwmh47wk75lozg3bum5yo35v3f/share/epic/epic_backward_hcal_only.xml";
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

    TH2D *h_sampling_fraction = new TH2D("h_sampling_fraction", "nHCal sampling fraction vs. energy; E [GeV]; sampling fraction; counts", 200, 0.0, 20.0, 500, 0.0, 0.05);
    TProfile *p_sampling_fraction = new TProfile("p_sampling_fraction", "nHCal sampling fraction vs. energy; E [GeV]; sampling fraction", 200, 0.0, 20.0, 0.0, 0.05);
    

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
        singlePart_Ekin = mcpart.getEnergy()-mcpart.getMass();
        
        for (const auto& hit : SimCalorimeterHit_coll) 
        {
            // if(mcpart.getGeneratorStatus()==1)
            // {
                const int hit_slice = decoder->get(hit.getCellID(), slice_index);
                hit_Esum += hit.getEnergy();
                if(hit_slice == 3) hit_scint_Esum += hit.getEnergy();
            //}
        }        

        h_sampling_fraction->Fill(singlePart_Ekin, hit_scint_Esum/hit_Esum);
		p_sampling_fraction->Fill(singlePart_Ekin, hit_scint_Esum/hit_Esum);
    }
    
    delete reader; 

    TCanvas *canvas = new TCanvas("canvas", "canvas", 1600, 800);
    canvas->Divide(1,1);
    h_sampling_fraction->Draw("COLZ");
    p_sampling_fraction->Draw("SAME");

    canvas->Print(outputFile);

    return 0;
}
