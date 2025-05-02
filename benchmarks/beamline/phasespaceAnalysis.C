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

void phasespaceAnalysis(  TString inFile      = "/scratch/EIC/G4out/beamline/beamlineTest.edm4hep.root",
                TString outFile     = "output.root",
                std::string compactName = "/home/simong/EIC/epic/install/share/epic/epic_ip6_extended.xml"){

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
               .Define("endID",getSubID("end",detector),{readoutName})
               .Define("pipeParameters",getVolumeParametersFromCellID(detector),{readoutName})
               .Define("pipeRadius",[](const ROOT::VecOps::RVec<volParams>& params) {
                ROOT::VecOps::RVec<double> radii;
                for (const auto& param : params) {
                    radii.push_back(param.radius);
                }
                return radii;
                }, {"pipeParameters"})
                .Define("xdet",[](const ROOT::VecOps::RVec<volParams>& params) {
                ROOT::VecOps::RVec<double> xPos;
                for (const auto& param : params) {
                    xPos.push_back(param.xPos);
                }
                return xPos;
                }, {"pipeParameters"})
                .Define("zdet",[](const ROOT::VecOps::RVec<volParams>& params) {
                ROOT::VecOps::RVec<double> zPos;
                for (const auto& param : params) {
                    zPos.push_back(param.zPos);
                }
                return zPos;
                }, {"pipeParameters"})
                .Define("rotation",[](const ROOT::VecOps::RVec<volParams>& params) {
                ROOT::VecOps::RVec<double> rotation;
                for (const auto& param : params) {
                    rotation.push_back(param.rotation);
                }
                return rotation;
                }, {"pipeParameters"});
                

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
                .Define("xmomMag","hitPosMom[3]")
                .Define("ymomMag","hitPosMom[4]")
                .Define("zmomMag","hitPosMom[5]")
                .Define("momMag","sqrt(xmomMag*xmomMag+ymomMag*ymomMag+zmomMag*zmomMag)")
                .Define("xmom","xmomMag/momMag")
                .Define("ymom","ymomMag/momMag")
                .Define("zmom","zmomMag/momMag");        

    }
    else{
        std::cout << "Collection " << readoutName << " not found in file" << std::endl;
        return;
    }    

    // Calculate the maximum pipe radius for plotting
    auto maxPipeRadius = 1.2*d1.Max("pipeRadius").GetValue();

    std::cout << "Executing Analysis and creating histograms" << std::endl;

    //Create array of histogram results
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsxy;
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsxyZoom;
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsxpx;
    std::map<TString,ROOT::RDF::RResultPtr<TH2D>> hHistsypy;

    std::map<TString,double> xMeans;
    std::map<TString,double> yMeans;
    std::map<TString,double> xStdDevs;
    std::map<TString,double> yStdDevs;
    std::map<TString,double> pxMeans;
    std::map<TString,double> pyMeans;
    std::map<TString,double> pxStdDevs;
    std::map<TString,double> pyStdDevs;

    //Fit paremeter and error maps
    std::map<TString,double> xMeanFit;
    std::map<TString,double> yMeanFit;
    std::map<TString,double> xMeanFitErr;
    std::map<TString,double> yMeanFitErr;
    std::map<TString,double> xStdDevFit;
    std::map<TString,double> yStdDevFit;
    std::map<TString,double> xStdDevFitErr;
    std::map<TString,double> yStdDevFitErr;
    std::map<TString,double> pxMeanFit;
    std::map<TString,double> pyMeanFit;
    std::map<TString,double> pxMeanFitErr;
    std::map<TString,double> pyMeanFitErr;
    std::map<TString,double> pxStdDevFit;
    std::map<TString,double> pyStdDevFit;
    std::map<TString,double> pxStdDevFitErr;
    std::map<TString,double> pyStdDevFitErr;

    std::map<TString,double> pipeRadii;
    std::map<TString,double> pipeXPos;
    std::map<TString,double> pipeZPos;
    std::map<TString,double> pipeRotation;

    auto xmin = d1.Min("xpos").GetValue();
    auto xmax = d1.Max("xpos").GetValue();
    auto ymin = d1.Min("ypos").GetValue();
    auto ymax = d1.Max("ypos").GetValue();
    auto pxmin = d1.Min("xmom").GetValue();
    auto pxmax = d1.Max("xmom").GetValue();
    auto pymin = d1.Min("ymom").GetValue();
    auto pymax = d1.Max("ymom").GetValue();
    
    //Create histograms
    for(int i=0; i<=7; i++){

        std::string name = "pipeID";
        name += std::to_string(i);
        auto filterDF = d1.Define("xposf","xpos["+std::to_string(i)+"]")
                          .Define("yposf","ypos["+std::to_string(i)+"]")
                          .Define("xmomf","xmom["+std::to_string(i)+"]")
                          .Define("ymomf","ymom["+std::to_string(i)+"]")
                          .Define("pipeRadiusf","pipeRadius["+std::to_string(i)+"]")
                          .Define("xdetf","xdet["+std::to_string(i)+"]")
                          .Define("zdetf","zdet["+std::to_string(i)+"]")
                          .Define("rotationf","rotation["+std::to_string(i)+"]");
                          

        //Calculate Min and Max values
        auto xminf = filterDF.Min("xposf").GetValue();
        auto xmaxf = filterDF.Max("xposf").GetValue();
        auto yminf = filterDF.Min("yposf").GetValue();
        auto ymaxf = filterDF.Max("yposf").GetValue();
        auto pxminf = filterDF.Min("xmomf").GetValue();
        auto pxmaxf = filterDF.Max("xmomf").GetValue();
        auto pyminf = filterDF.Min("ymomf").GetValue();
        auto pymaxf = filterDF.Max("ymomf").GetValue();
        // Calculate means and standard deviations
        xMeans[name] = filterDF.Mean("xposf").GetValue();
        yMeans[name] = filterDF.Mean("yposf").GetValue();
        xStdDevs[name] = filterDF.StdDev("xposf").GetValue();
        yStdDevs[name] = filterDF.StdDev("yposf").GetValue();
        pxMeans[name] = filterDF.Mean("xmomf").GetValue();
        pyMeans[name] = filterDF.Mean("ymomf").GetValue();
        pxStdDevs[name] = filterDF.StdDev("xmomf").GetValue();
        pyStdDevs[name] = filterDF.StdDev("ymomf").GetValue();

        // Calculate axes for zoomed beamspot histogram so that it is quare around the mean x and y
        double halfrange = std::max({xMeans[name]-xminf, xmaxf-xMeans[name], yMeans[name]-yminf, ymaxf-yMeans[name]});
        double xMinZoom = xMeans[name] - halfrange;
        double xMaxZoom = xMeans[name] + halfrange;
        double yMinZoom = yMeans[name] - halfrange;
        double yMaxZoom = yMeans[name] + halfrange;

        TString beamspotName = "Beamspot ID"+std::to_string(i)+";x offset [cm]; y offset [cm]";
        TString xyname = name+";x Offset [cm]; y Offset [cm]";
        TString xname = name+";x Offset [cm]; x trajectory component";
        TString yname = name+";y Offset [cm]; y trajectory component";
        hHistsxy[name] = filterDF.Histo2D({beamspotName,beamspotName,400,-maxPipeRadius,maxPipeRadius,400,-maxPipeRadius,maxPipeRadius},"xposf","yposf");
        hHistsxyZoom[name] = filterDF.Histo2D({xyname,xyname,100,xMinZoom,xMaxZoom,100,yMinZoom,yMaxZoom},"xposf","yposf");
        hHistsxpx[name]    = filterDF.Histo2D({xname,xname,400,xmin,xmax,400,pxmin,pxmax},"xposf","xmomf");
        hHistsypy[name]    = filterDF.Histo2D({yname,yname,400,ymin,ymax,400,pymin,pymax},"yposf","ymomf");

        //Parameters of the pipe
        pipeRadii[name]    = filterDF.Max("pipeRadiusf").GetValue();
        pipeXPos[name]     = filterDF.Max("xdetf").GetValue();
        pipeZPos[name]     = filterDF.Max("zdetf").GetValue();
        pipeRotation[name] = filterDF.Max("rotationf").GetValue();

        //Fit gaussian to the x, y, px and py distributions
        auto xhist = hHistsxy[name]->ProjectionX();
        auto yhist = hHistsxy[name]->ProjectionY();
        auto pxhist = hHistsxpx[name]->ProjectionY();
        auto pyhist = hHistsypy[name]->ProjectionY();
        xhist->Fit("gaus","Q");
        yhist->Fit("gaus","Q");
        pxhist->Fit("gaus","Q");
        pyhist->Fit("gaus","Q");
        //Get the fit parameters and errors
        xMeanFit[name] = xhist->GetFunction("gaus")->GetParameter(1);
        yMeanFit[name] = yhist->GetFunction("gaus")->GetParameter(1);
        xMeanFitErr[name] = xhist->GetFunction("gaus")->GetParError(1);
        yMeanFitErr[name] = yhist->GetFunction("gaus")->GetParError(1);
        xStdDevFit[name] = xhist->GetFunction("gaus")->GetParameter(2);
        yStdDevFit[name] = yhist->GetFunction("gaus")->GetParameter(2);
        xStdDevFitErr[name] = xhist->GetFunction("gaus")->GetParError(2);
        yStdDevFitErr[name] = yhist->GetFunction("gaus")->GetParError(2);
        pxMeanFit[name] = pxhist->GetFunction("gaus")->GetParameter(1);
        pyMeanFit[name] = pyhist->GetFunction("gaus")->GetParameter(1);
        pxMeanFitErr[name] = pxhist->GetFunction("gaus")->GetParError(1);
        pyMeanFitErr[name] = pyhist->GetFunction("gaus")->GetParError(1);
        pxStdDevFit[name] = pxhist->GetFunction("gaus")->GetParameter(2);
        pyStdDevFit[name] = pyhist->GetFunction("gaus")->GetParameter(2);
        pxStdDevFitErr[name] = pxhist->GetFunction("gaus")->GetParError(2);
        pyStdDevFitErr[name] = pyhist->GetFunction("gaus")->GetParError(2);
     
    }

    // Create histograms of the beamspot
    TCanvas *cXY = new TCanvas("beamspot_canvas","beamspot_canvas",3000,1600);
    cXY->Divide(4,2);
    int i=1;
    for(auto [name,h] : hHistsxy){
        // Get the pipe radius for this histogram
        auto pipeRadius = pipeRadii[name];
        cXY->cd(i++);

        h->Draw("col");
        //Draw cicle
        TEllipse *circle = new TEllipse(0,0,pipeRadius);
        circle->SetLineColor(kRed);
        circle->SetLineWidth(2);
        circle->SetFillStyle(0);
        circle->Draw("same");

        // Add zoomed version in the top-right corner
        TPad *pad = new TPad("zoomPad", "Zoomed View", 0.65, 0.65, 1.0, 1.0);
        pad->SetFillStyle(0); // Transparent background
        pad->Draw();
        pad->cd();
        // Draw the zoomed histogram without its title or axis titles
        hHistsxyZoom[name]->SetTitle("");
        hHistsxyZoom[name]->GetXaxis()->SetTitle("");
        hHistsxyZoom[name]->GetYaxis()->SetTitle("");
        hHistsxyZoom[name]->Draw("col");
        cXY->cd(); // Return to the main canvas
    }

    // x-px canvas
    TCanvas *cX = new TCanvas("x_px_canvas","x_px_canvas",3000,1600);
    cX->Divide(4,2);

    i=1;
    //Write histograms to file
    for(auto& h : hHistsxpx){
        cX->cd(i++);
        h.second->Draw("col");
    }

    // y-py canvas
    TCanvas *cY = new TCanvas("y_py_canvas","y_py_canvas",3000,1600);
    cY->Divide(4,2);

    i=1;
    for(auto& h : hHistsypy){
        cY->cd(i++);
        h.second->Draw("col");
    }

    // Save 2D canvases as pngs
    cXY->SaveAs("beamspot.png");
    cX->SaveAs("x_px.png");
    cY->SaveAs("y_py.png");

    // ---------------------------------------------------------------------------
    // Create histograms showing the fitted means and standard deviations of the positions and momenta
    // ---------------------------------------------------------------------------

    // Create histograms for fitted X means and standard deviations
    TH1F* hFittedXMeans = CreateFittedHistogram("hFittedXMeans", 
        "Mean X Offset [cm]", 
        xMeanFit, 
        xMeanFitErr, 
        "Pipe ID");

    TH1F* hFittedXStdDevs = CreateFittedHistogram("hFittedXStdDevs", 
        "Std Deviation X Offset [cm]", 
        xStdDevFit, 
        xStdDevFitErr, 
        "Pipe ID");

    // Create histograms for fitted Y means and standard deviations
    TH1F* hFittedYMeans = CreateFittedHistogram("hFittedYMeans", 
        "Mean Y Offset [cm]",
        yMeanFit, 
        yMeanFitErr, 
        "Pipe ID");

    TH1F* hFittedYStdDevs = CreateFittedHistogram("hFittedYStdDevs", 
        "Std Deviation Y Offset [cm]",
        yStdDevFit, 
        yStdDevFitErr, 
        "Pipe ID");

    TH1F* hFittedPxMeans = CreateFittedHistogram("hFittedPxMeans", 
        "Mean Px", 
        pxMeanFit, 
        pxMeanFitErr, 
        "Pipe ID");
    TH1F* hFittedPyMeans = CreateFittedHistogram("hFittedPyMeans",
        "Mean Py", 
        pyMeanFit, 
        pyMeanFitErr, 
        "Pipe ID");
    TH1F* hFittedPxStdDevs = CreateFittedHistogram("hFittedPxStdDevs",
        "Std Deviation Px", 
        pxStdDevFit, 
        pxStdDevFitErr, 
        "Pipe ID");
    TH1F* hFittedPyStdDevs = CreateFittedHistogram("hFittedPyStdDevs",
        "Std Deviation Py",
        pyStdDevFit,
        pyStdDevFitErr,
        "Pipe ID");


    // Create a canvas for the fitted beamspot means and standard deviations
    TCanvas *cFittedMeans = new TCanvas("cFittedMeans", "Fitted Beamspot Means and Std Deviation", 1200, 800);
    cFittedMeans->Divide(2, 2);
    cFittedMeans->cd(1);
    hFittedXMeans->Draw("E1"); // "E1" draws error bars
    cFittedMeans->cd(2);
    hFittedXStdDevs->Draw("E1"); // "E1" draws error bars
    cFittedMeans->cd(3);
    hFittedYMeans->Draw("E1"); // "E1" draws error bars
    cFittedMeans->cd(4);
    hFittedYStdDevs->Draw("E1"); // "E1" draws error bars
    cFittedMeans->SetGrid();
    cFittedMeans->Update();
    // Save the canvas as a PNG file
    cFittedMeans->SaveAs("fitted_means_stddevs.png");

    // Create a canvas for the fitted momentum means and standard deviations
    TCanvas *cFittedMomentumMeans = new TCanvas("cFittedMomentumMeans", "Fitted Momentum Means and Std Deviation", 1200, 800);
    cFittedMomentumMeans->Divide(2, 2);
    cFittedMomentumMeans->cd(1);
    hFittedPxMeans->Draw("E1"); // "E1" draws error bars
    cFittedMomentumMeans->cd(2);
    hFittedPxStdDevs->Draw("E1"); // "E1" draws error bars
    cFittedMomentumMeans->cd(3);
    hFittedPyMeans->Draw("E1"); // "E1" draws error bars
    cFittedMomentumMeans->cd(4);
    hFittedPyStdDevs->Draw("E1"); // "E1" draws error bars
    cFittedMomentumMeans->SetGrid();
    cFittedMomentumMeans->Update();
    // Save the canvas as a PNG file
    cFittedMomentumMeans->SaveAs("fitted_momentum_means_stddevs.png");


    // -----------------------------------------------------------------------------
    // Create histograms of the beampipe parameters
    // -----------------------------------------------------------------------------
    TH1F* hPipeRadii = CreateFittedHistogram("hPipeRadii",
        "Radius [cm]",
        pipeRadii,
        {},
        "Pipe ID");
    TH1F* hPipeXPos = CreateFittedHistogram("hPipeXPos",
        "X Position [cm]",
        pipeXPos,
        {},
        "Pipe ID");
    TH1F* hPipeZPos = CreateFittedHistogram("hPipeZPos",
        "Z Position [cm]",
        pipeZPos,
        {},
        "Pipe ID");
    TH1F* hPipeRotations = CreateFittedHistogram("hPipeRotations",
        "Rotation [rad]",
        pipeRotation,
        {},
        "Pipe ID");

    // Create a canvas for the pipe parameters
    TCanvas *cPipeParams = new TCanvas("cPipeParams", "Pipe Parameters", 1200, 400);    
    cPipeParams->Divide(4, 1);
    cPipeParams->cd(1);
    hPipeRadii->Draw("");
    cPipeParams->cd(2);
    hPipeXPos->Draw("");
    cPipeParams->cd(3); 
    hPipeZPos->Draw("");
    cPipeParams->cd(4);
    hPipeRotations->Draw("");
    cPipeParams->SetGrid();
    cPipeParams->Update();
    // Save the canvas as a PNG file
    cPipeParams->SaveAs("pipe_parameters.png");


    TFile *f = new TFile(outFile,"RECREATE");
    cXY->Write();
    cX->Write();
    cY->Write();
    cFittedMeans->Write();
    cFittedMomentumMeans->Write();
    cPipeParams->Write();

    f->Close();

    std::cout << "Saving events to file" << std::endl;

     ROOT::RDF::RSnapshotOptions opts;
    //  opts.fMode = "UPDATE";
    //  d1.Snapshot("events",outFile,{"pipeID","endID","pipeRadius","xpos","ypos","zpos","xmom","ymom","zmom","xpos_global","ypos_global","zpos_global","px_global","py_global","pz_global"},opts);
}
