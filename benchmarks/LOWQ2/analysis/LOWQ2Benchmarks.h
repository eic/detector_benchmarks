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

using RVecS       = ROOT::VecOps::RVec<string>;
using RVecI       = ROOT::VecOps::RVec<int>;

std::map<TString,std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>>> histMap;




// Create dataframe from input file(s)
RNode initialise( string fileName ){

  ROOT::RDataFrame d0("events",fileName);
  return d0;

}

// Format and write plots
void writePlots( TString outName ){

  TFile* rootfile = new TFile(outName,"RECREATE");
  auto benchmarkdir = rootfile->mkdir("LOWQ2");
  
  for(auto &[bmkey,hists]: histMap){

    auto histsdir   = benchmarkdir->mkdir(bmkey)->cd();

    for(auto &[key,hist]: hists.first){
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
    
    for(auto &[key,hist]: hists.second){
      
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
          // hist->Write(parts[i].c_str());
              
          // Make projection of 2D pixel binned histogram
          auto nBins      = hist->GetNcells();
          TString pixelFluxName = TString(parts[i].c_str())+"Flux";
          //Get maximum bin content to set range
          double maxBinContent = hist->GetMaximum();
          double logMax = log10(maxBinContent);

          TH1D* pixelFlux = new TH1D(pixelFluxName,pixelFluxName,100,0,logMax);
          for(int i=0; i<nBins; i++){
            pixelFlux->Fill(log10(hist->GetBinContent(i)));
          }
          pixelFlux->GetXaxis()->SetTitle("log(N_{hits})");
              
          hist->Write();
          pixelFlux->Write();
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
  }

  rootfile->Close();

}

void LOWQ2Benchmarks( string inName = "/scratch/EIC/G4out/qr_18x275_new.edm4hep*.root",
		      TString outName = "LOWQ2QRRates.root", dd4hep::Detector& detector=dd4hep::Detector::getInstance(), double eventRate=0.0 ){

  auto node = initialise( inName );

  RVecS colNames = node.GetColumnNames();

  std::string readoutName = "TaggerTrackerHits";

  if(Any(colNames==readoutName)){
    //-----------------------------------------
    // Hit detector IDs
    //-----------------------------------------
    auto ids = detector.readout(readoutName).idSpec().fields();
    for(auto &[key,id]: ids){
      TString colName = key+"ID";
      node = node.Define(colName,getSubID(key,detector),{readoutName});
    }
  }

  //Create Plots
  if(Any(colNames==readoutName)){
    histMap["SimDistributions"]  = createHitPlots(node,eventRate);
    
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