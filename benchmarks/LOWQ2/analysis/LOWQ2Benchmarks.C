#include "TString.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "LOWQ2hits.h"
#include "LOWQ2clusters.h"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;

using RVecS       = ROOT::VecOps::RVec<string>;
using RVecI       = ROOT::VecOps::RVec<int>;

// std::map<TString,H1ResultPtr> hHists1D;
// std::map<TString,H2ResultPtr> hHists2D;

std::map<TString,std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>>> histMap;


// Create dataframe from input file(s)
RNode initialise( string fileName = "/scratch/EIC/ReconOut/tempEvents10x100-compare.root" ){

  ROOT::RDataFrame d0("events",fileName);
  return d0;

}


// Format and write plots
void writePlots( TString outName = "LOWQ2Benchmarks.root"){

  TFile* rootfile = new TFile(outName,"RECREATE");
  auto benchmarkdir = rootfile->mkdir("LOWQ2");
  
  for(auto &[bmkey,hists]: histMap){

    auto histsdir   = benchmarkdir->mkdir(bmkey)->cd();

    for(auto &[key,hist]: hists.first){
      hist->Write();
    }
    
    for(auto &[key,hist]: hists.second){
      
      // Make projection of 2D pixel binned histogram
      auto nBins      = hist->GetNcells();
      TString pixelFluxName = TString(hist->GetTitle())+"Flux";
      TH1D* pixelFlux = new TH1D(pixelFluxName,pixelFluxName,20,0,20);
      for(int i=0; i<nBins; i++){
	pixelFlux->Fill(hist->GetBinContent(i));
      }
      
      hist->Write();
      pixelFlux->Write();
      
    }
  }

  rootfile->Close();

}

// Main method
void LOWQ2Benchmarks(){

  auto node = initialise();

  RVecS colNames = node.GetColumnNames();

  if(Any(colNames=="TaggerTrackerHits")){
    histMap["SimDistributions"] = createHitPlots(node);
  }

  if(Any(colNames=="TaggerTrackerClusterPositions")){  
    histMap["ClusterDistributions"] =  createClusterPlots(node);
  }

  writePlots();

}
