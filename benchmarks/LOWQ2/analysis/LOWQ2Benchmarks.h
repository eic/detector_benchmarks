#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "LOWQ2hits.h"
#include "LOWQ2acceptance.h"
#include "LOWQ2clusters.h"
#include "LOWQ2reconstruction.h"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using H3ResultPtr = ROOT::RDF::RResultPtr<TH3D>;

using RVecS       = ROOT::VecOps::RVec<string>;
using RVecI       = ROOT::VecOps::RVec<int>;

std::map<TString,std::tuple<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>,std::map<TString,H3ResultPtr>>> histMap;



//---------------------------------------------------------------------------------------------
// Create dataframe from input file(s)
//---------------------------------------------------------------------------------------------
RNode initialise( string fileName ){

  ROOT::RDataFrame d0("events",fileName);
  return d0;

}

//---------------------------------------------------------------------------------------------
// Create histograms showing the variation of flux in various system components
//---------------------------------------------------------------------------------------------
void WriteFluxPlots(H2ResultPtr hist, std::string parts){
  // Make projection of 2D pixel binned histogram
  auto nBins      = hist->GetNcells();
  TString pixelFluxName = TString(parts)+"Flux";
  TString logPixelFluxName = TString(parts)+"LogFlux";

  //Get maximum bin content to set range
  double maxBinContent = hist->GetMaximum();
  double logMax = log10(maxBinContent);

  //Create pixel flux histogram
  TH1D* pixelFlux = new TH1D(pixelFluxName,pixelFluxName,100,0,maxBinContent);

  //Create log pixel flux histogram
  TH1D* logPixelFlux = new TH1D(logPixelFluxName,logPixelFluxName,100,0,logMax);
  for(int i=0; i<nBins; i++){
    pixelFlux->Fill(hist->GetBinContent(i));
    logPixelFlux->Fill(log10(hist->GetBinContent(i)));
  }
  pixelFlux->GetXaxis()->SetTitle("N_{hits}");
  logPixelFlux->GetXaxis()->SetTitle("log(N_{hits})");
  pixelFlux->GetYaxis()->SetTitle("N_{pixels}");
  logPixelFlux->GetYaxis()->SetTitle("N_{pixels}");
      
  pixelFlux->Write();
  logPixelFlux->Write();

}

//---------------------------------------------------------------------------------------------
//Create histograms showing the variation of resolution accross the kinematic phase space
//---------------------------------------------------------------------------------------------
void WriteResolutionPlots(H2ResultPtr hist, std::string parts){

  hist->FitSlicesX(nullptr,0,-1,5);
  TH1D* h1 = (TH1D*)gDirectory->Get(TString(hist->GetName())+"_2");
  h1->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  h1->GetYaxis()->SetTitle(TString(hist->GetYaxis()->GetTitle())+"Resolution");

  h1->Write();

}


//---------------------------------------------------------------------------------------------
// Format and write plots
//---------------------------------------------------------------------------------------------
void writePlots( TString outName ){

  TFile* rootfile = new TFile(outName,"RECREATE");
  auto benchmarkDir = rootfile->mkdir("LOWQ2");
  
  //Loop over 1D histograms saving
  for(auto &[bmkey,hists]: histMap){

    //Split mapped distributions name into parts
    std::stringstream ss(bmkey.Data());
    std::string part;
    std::vector<std::string> parts;
    while (std::getline(ss, part, '/')) {
      parts.push_back(part);
    }

    TDirectory* dir = benchmarkDir;
    for (size_t i = 0; i < parts.size(); ++i) {
      if (i == parts.size() - 1) {
        // This is the last part, write the histogram
        dir->cd();
      } else {
        // This is not the last part, create or get the directory
        if (!dir->GetDirectory(parts[i].c_str())) {
          dir = dir->mkdir(parts[i].c_str());
        } else {
          dir = dir->GetDirectory(parts[i].c_str());
        }
        dir->cd();
      }
    }
        
    TDirectory* currentDir = gDirectory;

    for(auto &[key,hist]: get<0>(hists)){

      //Split histogram name into parts
      std::stringstream ss(key.Data());
      std::string part;
      std::vector<std::string> parts;
      while (std::getline(ss, part, '/')) {
        parts.push_back(part);
      }

      TDirectory* dir = currentDir;
      for (size_t i = 0; i < parts.size(); ++i) {
        if (i == parts.size() - 1) {
          // This is the last part, write the histogram
          hist->SetMinimum(0);
          hist->Write(parts[i].c_str());
        } else {
          // This is not the last part, create or get the directory
          if (!dir->GetDirectory(parts[i].c_str())) {
            dir = dir->mkdir(parts[i].c_str());
          } else {
            dir = dir->GetDirectory(parts[i].c_str());
          }
          dir->cd();
        }
      }
      currentDir->cd();
      // hist->Write(key);
    }
    
    //Loop over 2D histograms saving and calculating fluxes as needed
    for(auto &[key,hist]: get<1>(hists)){
      
      TDirectory* currentDir = gDirectory;
      std::stringstream ss(key.Data());
      std::string part;
      TDirectory* dir = currentDir;
      std::vector<std::string> parts;
      while (std::getline(ss, part, '/')) {
        parts.push_back(part);
      }
      for (size_t i = 0; i < parts.size(); ++i) {
        if (i == parts.size() - 1) {          
          // If histogram name contains rate or hits then calculate flux
          if(parts[i].find("Rate") != std::string::npos || parts[i].find("Hits") != std::string::npos){
            WriteFluxPlots(hist,parts[i].c_str());            
          }

          // If histogram name contains Res then create a profile
          if(parts[i].find("Res") != std::string::npos){
            WriteResolutionPlots(hist,parts[i].c_str());
          }

          hist->Sumw2(false);
          hist->Write();
        } else {
          // This is not the last part, create or get the directory
          if (!dir->GetDirectory(parts[i].c_str())) {
            dir = dir->mkdir(parts[i].c_str());
          } else {
            dir = dir->GetDirectory(parts[i].c_str());
          }
          dir->cd();
        }
      }
      currentDir->cd();      
    }

    //Loop over 3D histograms saving
    for(auto &[key,hist]: get<2>(hists)){
      TDirectory* currentDir = gDirectory;
      std::stringstream ss(key.Data());
      std::string part;
      TDirectory* dir = currentDir;
      std::vector<std::string> parts;
      while (std::getline(ss, part, '/')) {
        parts.push_back(part);
      }
      for (size_t i = 0; i < parts.size(); ++i) {
        if (i == parts.size() - 1) {
          // This is the last part, write the histogram
          hist->Write(parts[i].c_str());

          //Fit histogram z slices and save 2d histogram
          hist->FitSlicesZ(nullptr,0,-1,5);
          TH2D* h2 = (TH2D*)gDirectory->Get(TString(hist->GetName())+"_2");
          h2->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
          h2->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
          h2->GetZaxis()->SetTitle(TString(hist->GetZaxis()->GetTitle())+"Resolution");
          h2->SetTitle(hist->GetTitle());
          h2->Write();          

        } else {
          // This is not the last part, create or get the directory
          if (!dir->GetDirectory(parts[i].c_str())) {
            dir = dir->mkdir(parts[i].c_str());
          } else {
            dir = dir->GetDirectory(parts[i].c_str());
          }
          dir->cd();
        }
      }
      currentDir->cd();
      // hist->Write(key);
    }

    benchmarkDir->cd();

  }

  rootfile->Close();

}

//---------------------------------------------------------------------------------------------
// Create the benchmark plots
//---------------------------------------------------------------------------------------------
void LOWQ2Benchmarks( string inName = "/scratch/EIC/G4out/qr_18x275_new.edm4hep*.root",
		      TString outName = "LOWQ2QRRates.root", dd4hep::Detector& detector=dd4hep::Detector::getInstance(), double eventRate=0.0 ){

  auto node = initialise( inName );
  
  int events = *node.Count();
  double eventWeight = eventRate / events;

  RVecS colNames = node.GetColumnNames();

  std::string readoutName = "TaggerTrackerHits";

  if(Any(colNames==readoutName)){
    //-----------------------------------------
    // Hit detector IDs
    //-----------------------------------------
    auto ids = detector.readout(readoutName).idSpec().fields();

    node = node.Define("hitParticleIndex","_TaggerTrackerHits_MCParticle.index")
               .Define("generatorStat","MCParticles.generatorStatus")
               .Define("primary",[](RVecI index, RVecI status){
                        RVecI result(index.size());
                        for (size_t i = 0; i < index.size(); ++i) {
                            result[i] = status[index[i]] == 1;
                        }
                        return result;
                },{"hitParticleIndex","generatorStat"});

    std::string filterName = "TaggerTrackerHitsFilter";
    auto allnode       = node.Define(filterName,"TaggerTrackerHits");
    auto primarynode   = node.Define(filterName,"TaggerTrackerHits[primary==1]");
    auto secondarynode = node.Define(filterName,"TaggerTrackerHits[primary==0]");

    for(auto &[key,id]: ids){
      TString colName = key+"ID";
      node          = node.Define(colName,getSubID(key,detector),{readoutName});
      primarynode   = primarynode.Define(colName,getSubID(key,detector),{filterName});
      secondarynode = secondarynode.Define(colName,getSubID(key,detector),{filterName});
    }

    //Create Rate Plots
    histMap["Rates/AllHits"]        = createHitPlots(node,eventWeight);
    histMap["Rates/PrimaryHits"]    = createHitPlots(primarynode,eventWeight);
    histMap["Rates/SecondaryHits"]  = createHitPlots(secondarynode,eventWeight); 

  }

  if((Any(colNames==readoutName) || Any(colNames=="InclusiveKinematicsElectron"))  && Any(colNames=="MCParticles")){  
    histMap["AcceptanceDistributions"] = createAcceptancePlots(node);
  }

  if(Any(colNames=="TaggerTrackerClusterPositions")){  
    histMap["ClusterDistributions"] =  createClusterPlots(node);
  }

  if(Any(colNames=="LowQ2TrackParameters") && Any(colNames=="MCParticles")){  
    histMap["ReconstructedDistributions"] = createReconstructionPlots(node);
  }

  writePlots( outName );

}