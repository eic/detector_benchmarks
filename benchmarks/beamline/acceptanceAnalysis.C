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
#include "shared_functions.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEllipse.h"

using RVecS       = ROOT::VecOps::RVec<string>;
using RNode       = ROOT::RDF::RNode;

int acceptanceAnalysis( TString inFile             = "/home/simong/EIC/detector_benchmarks_anl/sim_output/beamline/acceptanceTestXS2.edm4hep.root",
                        TString outFile            = "output.root",
                        std::string compactName    = "/home/simong/EIC/epic/install/share/epic/epic_ip6_extended.xml",
                        TString beampipeCanvasName = "acceptance_in_beampipe.png",
                        TString EThetaCanvasName   = "acceptance_energy_theta.png",
                        TString EThetaAccCanvasName= "acceptance_energy_theta_acceptance.png",
                        TString entryFractionCanvasName = "acceptance_entries.png") {

    //Set ROOT style    
    gStyle->SetPadLeftMargin(0.1);  // Set left margin
    gStyle->SetPadRightMargin(0.0); // Set right margin
    gStyle->SetPadTopMargin(0.0);   // Set top margin
    gStyle->SetPadBottomMargin(0.1);// Set bottom margin
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleX(0.12);          // Place the title on the top right
    gStyle->SetTitleY(0.985);         // Place the title on the top right
    gStyle->SetTitleSize(0.08, "t");
    gStyle->SetTitleXOffset(1.0);    // Adjust y-axis title offset
    gStyle->SetTitleYOffset(1.0);    // Adjust y-axis title offset
    gStyle->SetOptStat(0);

    int pass = 0;

    //Set implicit multi-threading
    // ROOT::EnableImplicitMT();
       
    //Load the detector config
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(compactName);
 
    ROOT::RDataFrame d0("events",inFile, {"BackwardsBeamlineHits"});
    RNode d1 = d0;
    RVecS colNames = d1.GetColumnNames();

    //Get total number of entries
    int nEntries = d1.Count().GetValue();
    //Set number of entries to process
    // d1 = d1.Range(0,1000);
    
    //Get the collection 
    std::string mcParticlesName = "MCParticles";
    std::string readoutName = "BackwardsBeamlineHits";  

    std::cout << "Running lazy RDataframe execution" << std::endl;

    if(Any(colNames==mcParticlesName)){

        // Get theta and energy of the particles
        d1 = d1.Define("theta",[](const vector<edm4hep::MCParticleData>& mcParticles) {
            // Calculate theta from momentum components as angle from negative z axis
            double px = mcParticles[0].momentum.x;
            double py = mcParticles[0].momentum.y;
            double pz = mcParticles[0].momentum.z;
            double p = std::sqrt(px*px + py*py + pz*pz);
            double theta = M_PI-std::acos(pz / p); // Angle from the z-axis
            
            return theta;
        }, {mcParticlesName})
        .Define("energy",[](const vector<edm4hep::MCParticleData>& mcParticles) {
           
            //Calculate energy from mass and momentum
            double mass = mcParticles[0].mass;
            double px = mcParticles[0].momentum.x;
            double py = mcParticles[0].momentum.y;
            double pz = mcParticles[0].momentum.z;
            double energy = std::sqrt(px*px + py*py + pz*pz + mass*mass);
            
            return energy;
        }, {mcParticlesName});

    }
    else{
        std::cout << "Collection " << mcParticlesName << " not found in file" << std::endl;
        return 1;
    }

    int eBins = 100;
    int thetaBins = 100;

    auto totalETheta = d1.Histo2D({"Energy vs Theta","Energy vs Theta",eBins,6,18,thetaBins,0,0.011},"energy","theta");

    if(Any(colNames==readoutName)){

        d1 = d1.Define("pipeID",getSubID("pipe",detector),{readoutName})
               .Define("endID",getSubID("end",detector),{readoutName})
               .Define("pipeParameters",getVolumeParametersFromCellID(detector),{readoutName})
               .Define("pipeRadius",[](const ROOT::VecOps::RVec<volParams>& params) {
                ROOT::VecOps::RVec<double> radii;
                for (const auto& param : params) {
                    radii.push_back(param.radius);
                }
                return radii;
                }, {"pipeParameters"});
                

        //global x,y,z position and momentum
        d1 = d1 .Define("NHits","BackwardsBeamlineHits.size()");
        
        d1 = d1.Define("hitPosMom",globalToLocal(detector),{readoutName})
                .Define("xpos","hitPosMom[0]")
                .Define("ypos","hitPosMom[1]")
                .Define("zpos","hitPosMom[2]");        

    }
    else{
        std::cout << "Collection " << readoutName << " not found in file" << std::endl;
        return 1;
    }    

    // Calculate the maximum pipe radius for plotting
    auto maxPipeRadius = 2*d1.Max("pipeRadius").GetValue();

    std::cout << "Executing Analysis and creating histograms" << std::endl;

    //Create array of histogram results
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsxy;
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsETheta;
    

    std::map<TString,ROOT::RDF::RResultPtr<double>> pipeRadii;
    std::map<TString,double> filterEntries;
    std::map<TString,double> filterAcceptanceIntegral;
    
    //Create histograms
    for(int i=0; i<=7; i++){

        std::string name = "pipeID";
        std::string str_i = std::to_string(i);
        name += str_i;
        auto filterDF = d1.Define("xposf","xpos[pipeID=="+str_i+"]")
                          .Define("yposf","ypos[pipeID=="+str_i+"]")
                          .Define("pipeRadiusf","pipeRadius[pipeID=="+str_i+"]");
                   

        TString beamspotName = "Beamspot ID"+str_i+";x offset [cm]; y offset [cm]";
        TString xyname = name+";x Offset [cm]; y Offset [cm]";
        TString xname = name+";x Offset [cm]; x trajectory component";
        TString yname = name+";y Offset [cm]; y trajectory component";
        hHistsxy[name] = filterDF.Histo2D({beamspotName,beamspotName,400,-maxPipeRadius,maxPipeRadius,400,-maxPipeRadius,maxPipeRadius},"xposf","yposf");

        auto extraFilterDF = filterDF.Filter("All(pipeID<"+std::to_string(i+1)+")&&Any(pipeID=="+std::to_string(i)+")" );
        TString EThetaName = "Energy vs Theta ID"+str_i+";Energy [GeV]; Theta [rad]";
        TString nameETheta = name+";Energy [GeV]; Theta [rad]";
        hHistsETheta[name] = extraFilterDF.Histo2D({EThetaName,EThetaName,eBins,2,18,thetaBins,0,0.011},"energy","theta");

        //Parameters of the pipe
        pipeRadii[name]    = filterDF.Max("pipeRadiusf");        
        // std::cout << "Pipe ID: " << name << " Radius: " << pipeRadii[name] << " " << filterDF.Min("pipeRadiusf").GetValue() << std::endl;
    
    }

    // Create histograms of the beamspot
    TCanvas *cXY = new TCanvas("beamspot_canvas","beamspot_canvas",3000,1600);
    cXY->Divide(4,2);
    int i=1;
    for(auto [name,h] : hHistsxy){
        // Get the pipe radius for this histogram
        auto pipeRadius = pipeRadii[name].GetValue();
        cXY->cd(i++);

        h->Draw("col");
        //Draw cicle
        TEllipse *circle = new TEllipse(0,0,pipeRadius);
        circle->SetLineColor(kRed);
        circle->SetLineWidth(2);
        circle->SetFillStyle(0);
        circle->Draw("same");

    }

    // Cnavas for energy vs theta
    TCanvas *cETheta = new TCanvas("energy_theta_canvas","energy_theta_canvas",3000,1600);
    cETheta->Divide(4,2);
    i=1;
    for(auto [name,h] : hHistsETheta){        
        cETheta->cd(i++);
        h->Draw("colz");
        filterEntries[name] = h->GetEntries()/ nEntries;
    }
  
    // Canvas for energy vs theta acceptance
    TCanvas *cEThetaAcc = new TCanvas("energy_theta_acceptance_canvas","energy_theta_acceptance_canvas",3000,1600);
    cEThetaAcc->Divide(4,2);
    i=1;
    for(auto [name,h] : hHistsETheta){
        cEThetaAcc->cd(i++);
        // Clone the histogram before dividing to avoid modifying the original
        TH2D* hAcc = (TH2D*)h->Clone(Form("%s_acceptance", name.Data()));
        hAcc->Divide(totalETheta.GetPtr());
        hAcc->SetMaximum(1);
        hAcc->Draw("colz");
        filterAcceptanceIntegral[name] = hAcc->Integral()/ (eBins*thetaBins); // Normalize by the number of bins
    }

    TH1F* hPipeEntries = CreateFittedHistogram("hNumberOfHits",
        "Number of Hits per Pipe ID",
        filterEntries,
        {},
        "Pipe ID");
    TH1F* hPipeAcceptance = CreateFittedHistogram("hPipeAcceptance",
        "Acceptance per Pipe ID",
        filterAcceptanceIntegral,
        {},
        "Pipe ID");

    TCanvas *cPipeAcceptance = new TCanvas("cPipeAcceptance", "Pipe Acceptance", 1200, 400);    
    cPipeAcceptance->Divide(2, 1);
    cPipeAcceptance->cd(1);
    hPipeEntries->Draw("");
    cPipeAcceptance->cd(2);
    hPipeAcceptance->Draw("");
    cPipeAcceptance->SetGrid();
    cPipeAcceptance->Update();

    // Save 2D canvases as pngs
    cXY->SaveAs(beampipeCanvasName);
    cETheta->SaveAs(EThetaCanvasName);
    cEThetaAcc->SaveAs(EThetaAccCanvasName);
    cPipeAcceptance->SaveAs(entryFractionCanvasName);

    TFile *f = new TFile(outFile,"RECREATE");
    cXY->Write();
    cETheta->Write();
    cEThetaAcc->Write();
    cPipeAcceptance->Write();

    f->Close();

    std::cout << "Analysis complete. Results saved to " << outFile << std::endl;
    return pass;

    // std::cout << "Saving events to file" << std::endl;

    //  ROOT::RDF::RSnapshotOptions opts;
    // opts.fMode = "UPDATE";
    // d1.Snapshot("events",outFile,{"pipeID","endID","pipeRadius","xpos","ypos","zpos","xmom","ymom","zmom","xpos_global","ypos_global","zpos_global","px_global","py_global","pz_global"},opts);
}
