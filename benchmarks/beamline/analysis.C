// Script to plot the x and y positions and momentum of beamline particles as they pass through the magnets

#include <iostream>
#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/Vector3f.h"
#include "edm4hep/Vector3d.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RVec.hxx"
#include "functors.h"
#include "TCanvas.h"
#include "TStyle.h"

using RVecS       = ROOT::VecOps::RVec<string>;
using RNode       = ROOT::RDF::RNode;

void analysis(  TString inFile      = "/scratch/EIC/G4out/beamline/beamlineTest.edm4hep.root",
                TString outFile     = "output.root",
                std::string compactName = "/home/simong/EIC/epic/install/share/epic/epic_ip6_extended.xml"){

    //Set ROOT style    
    gStyle->SetPadLeftMargin(0.25);  // Set left margin
    gStyle->SetPadRightMargin(0.15); // Set right margin
    gStyle->SetPadTopMargin(0.05);   // Set top margin
    gStyle->SetPadBottomMargin(0.15);// Set bottom margin
    gStyle->SetTitleYOffset(2);    // Adjust y-axis title offset
    gStyle->SetOptStat(1110);

    //Set implicit multi-threading
    ROOT::EnableImplicitMT();
       
    //Load the detector config
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(compactName);
 
    ROOT::RDataFrame d0("events",inFile, {"BackwardsBeamlineHits"});
    RNode d1 = d0;
    RVecS colNames = d1.GetColumnNames();

    //Set number of entries to process
    // d1 = d1.Range(0,1000);
    
    //Get the collection 
    std::string readoutName = "BackwardsBeamlineHits";

    std::cout << "Running lazy RDataframe execution" << std::endl;

    if(Any(colNames==readoutName)){

        d1 = d1.Define("pipeID",getSubID("pipe",detector),{readoutName})
               .Define("endID",getSubID("end",detector),{readoutName});

        //global x,y,z position and momentum
        d1 = d1.Define("xpos_global","BackwardsBeamlineHits.position.x")
                .Define("ypos_global","BackwardsBeamlineHits.position.y")
                .Define("zpos_global","BackwardsBeamlineHits.position.z")
                .Define("px_global","BackwardsBeamlineHits.momentum.x")
                .Define("py_global","BackwardsBeamlineHits.momentum.y")
                .Define("pz_global","BackwardsBeamlineHits.momentum.z");

        d1 = d1.Define("hitPosMom",globalToLocal(detector),{readoutName})
                .Define("xpos","hitPosMom[0]")
                .Define("ypos","hitPosMom[1]")
                .Define("zpos","hitPosMom[2]")
                .Define("xmom","hitPosMom[3]")
                .Define("ymom","hitPosMom[4]")
                .Define("zmom","hitPosMom[5]");        

    }
    else{
        std::cout << "Collection " << readoutName << " not found in file" << std::endl;
        return;
    }    

    std::cout << "Executing Analysis and creating histograms" << std::endl;

    //Create array of histogram results
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsxpx;
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsypy;


    //Create histograms
    for(int i=0; i<=7; i++){

        std::string name = "pipeID";
        name += std::to_string(i);
        auto filterDF = d1.Define("xposf","xpos["+std::to_string(i)+"]")
                          .Define("yposf","ypos["+std::to_string(i)+"]")
                          .Define("xmomf","xmom["+std::to_string(i)+"]")
                          .Define("ymomf","ymom["+std::to_string(i)+"]");

        //Calculate Min and Max values
        auto xmin = filterDF.Min("xposf").GetValue();
        auto xmax = filterDF.Max("xposf").GetValue();
        auto ymin = filterDF.Min("yposf").GetValue();
        auto ymax = filterDF.Max("yposf").GetValue();
        auto pxmin = filterDF.Min("xmomf").GetValue();
        auto pxmax = filterDF.Max("xmomf").GetValue();
        auto pymin = filterDF.Min("ymomf").GetValue();
        auto pymax = filterDF.Max("ymomf").GetValue();

        TString xname = name+";x offset [mm]; x trajectory component";
        TString yname = name+";y offset [mm]; y trajectory component";
        hHistsxpx[name] = filterDF.Histo2D({xname,xname,100,xmin,xmax,100,pxmin,pxmax},"xposf","xmomf");
        hHistsypy[name] = filterDF.Histo2D({yname,yname,100,ymin,ymax,100,pymin,pymax},"yposf","ymomf");

    }

    TCanvas *cX = new TCanvas("x_px_canvas","x_px_canvas",3000,1600);
    cX->Divide(4,2);

    int i=1;
    //Write histograms to file
    for(auto& h : hHistsxpx){
        cX->cd(i++);
        h.second->Draw("colz");
    }


    TCanvas *cY = new TCanvas("y_py_canvas","y_py_canvas",3000,1600);
    cY->Divide(4,2);

    i=1;
    for(auto& h : hHistsypy){
        cY->cd(i++);
        h.second->Draw("colz");
    }

    TFile *f = new TFile(outFile,"RECREATE");
    cX->Write();
    cY->Write();

    f->Close();

    std::cout << "Saving events to file" << std::endl;

    // ROOT::RDF::RSnapshotOptions opts;
    // opts.fMode = "UPDATE";
    // d1.Snapshot("events",outFile,{"pipeID","endID","xpos","ypos","zpos","xmom","ymom","zmom","xpos_global","ypos_global","zpos_global","px_global","py_global","pz_global"},opts);
}
