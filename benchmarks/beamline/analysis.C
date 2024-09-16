// Script to plot the x and y positions and momentum of beamline particles as they pass through the magnets

#include <iostream>
#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/SimCalorimeterHit.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RVec.hxx"
#include "functors.h"

using RVecS       = ROOT::VecOps::RVec<string>;
using RNode       = ROOT::RDF::RNode;

void analysis(  std::string inFile      = "/scratch/EIC/G4out/beamline/beamlineTest.edm4hep.root",
                std::string compactName = "/home/simong/EIC/epic/install/share/epic/epic_ip6_extended.xml"){

    //Set implicit multi-threading
    ROOT::EnableImplicitMT();

    //Load the detector config
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(compactName);
 
    ROOT::RDataFrame d0("events",inFile, {"BackwardsBeamlineHits"});
    RNode d1 = d0;
    RVecS colNames = d1.GetColumnNames();
    
    //Get the collection 
    std::string readoutName = "BackwardsBeamlineHits";

    if(Any(colNames==readoutName)){

        auto ids = detector.readout(readoutName).idSpec().fields();

        for(auto &[key,id]: ids){
            TString colName = key+"ID";
            d1              = d1.Define(colName,getSubID(key,detector),{readoutName});
        }
    }
    else{
        std::cout << "Collection " << readoutName << " not found in file" << std::endl;
        return;
    }
    
    d1.Snapshot("events","output.root");

}
