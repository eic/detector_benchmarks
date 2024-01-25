#pragma once

#include "functors.h"
#include "DD4hep/Detector.h" // Add the missing include

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using RVecI       = ROOT::VecOps::RVec<int>;

// Lazy execution methods
std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>> createHitPlots( RNode d1, string compactName = "/home/simong/EIC/epic/epic_18x275.xml" ){

  int xChip = 448;
  int yChip = 512;
  std::map<int,int> xPixelMin = {{1,-xChip*3},{2,-xChip*3}};
  std::map<int,int> xPixelMax = {{1, xChip*3},{2, xChip*3}};

  std::map<int,int> yPixelMin = {{1,-yChip*3},{2,-yChip*3}};
  std::map<int,int> yPixelMax = {{1, yChip*3},{2, yChip*3}};


  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
    

  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compactName);
  //-----------------------------------------
  // Hit detector IDs
  //-----------------------------------------
  auto ids = detector.readout("TaggerTrackerHits").idSpec().fields();
  for(auto &[key,id]: ids){
    TString colName = key+"ID";
    d1 = d1.Define(colName,getSubID(key,detector),{"TaggerTrackerHits"});
  }

  d1 = d1.Define("boardID",[xChip](RVecI xID){ return (xID + 3*xChip) / (2 * xChip); },{"xID"})
        .Define("xChipID", [xChip](RVecI xID){ return (xID + 3*xChip) / (xChip); },{"xID"})
        .Define("xColumnID", [xChip](RVecI xID){ return (xID + 3*xChip) / 2; },{"xID"})
        .Define("yChipID", [yChip](RVecI yID){ return (yID + 3*yChip) / (yChip); },{"yID"});

  hHists1D["hmoduleID"] = d1.Histo1D({"hmoduleID",  "hmoduleID", 3,     0,    3   }, "moduleID");

  // Module Loop
  for(int i=1; i<=2; i++){

    double xMin   = xPixelMin[i];
    double xMax   = xPixelMax[i];
    double yMin   = yPixelMin[i];
    double yMax   = yPixelMax[i];
    int    xRange = xMax-xMin;
    int    yRange = yMax-yMin;

    TString modNum         = std::to_string(i);
    TString moduleTag      = "module"+modNum;
    TString layerName      = "layerID"+moduleTag;
    TString layerHistName  = "h"+layerName;
    

    auto d2 = d1.Define("ModuleFilter",[i](RVecI modID){return modID==i;},{"moduleID"})
      .Define(layerName,"layerID[ModuleFilter]");
        
    hHists1D[layerHistName]  = d2.Histo1D({layerHistName,   layerHistName,  4,     0,    4   }, layerName );
    
    // Layer Loop
    for(int j=0; j<=3; j++){
      
      TString layerNum      = std::to_string(j);
      TString layerTag      = "layer"+layerNum;
      
      TString xName         = "xID"+moduleTag+layerTag;
      TString xHistName     = "h"+xName;

      TString yName         = "yID"+moduleTag+layerTag;
      TString yHistName     = "h"+yName;

      TString xChipName     = "xChipID"+moduleTag+layerTag;
      TString yChipName     = "yChipID"+moduleTag+layerTag;

      TString xyHistName    = "h"+xName+yName;

      TString xColumnName    = "xColumnID"+moduleTag+layerTag;

      TString boardName     = "boardID"+moduleTag+layerTag;

      std::vector<string>  layerSizeInput = {xName.Data()};
      TString layerSizeName = "HitsPerEvent"+moduleTag+layerTag;
      TString sizeHistName  = "h"+layerSizeName;
      
      auto d3 = d2.Define("LayerFilter",[j](RVecI layID){return layID==j;},{"layerID"})
	                .Define(xName,"xID[LayerFilter&&ModuleFilter]")
	                .Define(yName,"yID[LayerFilter&&ModuleFilter]")
                  .Define(boardName,"boardID[LayerFilter&&ModuleFilter]")
                  .Define(xChipName,"xChipID[LayerFilter&&ModuleFilter]")
                  .Define(yChipName,"yChipID[LayerFilter&&ModuleFilter]")
                  .Define(xColumnName,"xColumnID[LayerFilter&&ModuleFilter]")
                  .Define(layerSizeName,[](RVecI lay){return lay.size();},layerSizeInput);

      
    
      hHists1D[xHistName]    = d3.Histo1D({xHistName,   xHistName,  xRange, xMin, xMax   }, xName );
      hHists1D[yHistName]    = d3.Histo1D({yHistName,   yHistName,  yRange, yMin, yMax   }, yName );
            
      hHists2D[xyHistName]   = d3.Histo2D({xyHistName,  xyHistName, xRange, xMin, xMax, yRange, yMin, yMax}, xName, yName );
      
      hHists1D[sizeHistName] = d3.Histo1D({sizeHistName,   sizeHistName,  100, 0, 100   }, layerSizeName );

      // Plot number of hits per boardID for each layer
      TString boardIDHistName = "hBoardID" +moduleTag + layerTag;
      hHists1D[boardIDHistName] = d3.Histo1D({boardIDHistName, boardIDHistName, 3, 0, 3}, boardName);
      
      // Plot number of hits per chipID for each layer
      TString chipIDHistName = "hChipID" +moduleTag + layerTag;
      hHists2D[chipIDHistName] = d3.Histo2D({chipIDHistName, chipIDHistName, 6, 0, 6, 6, 0, 6}, xChipName, yChipName);

      // Plot number of hits per chipID for each layer
      TString xColumnIDHistName = "hxColumnID" +moduleTag + layerTag;
      hHists2D[xColumnIDHistName] = d3.Histo2D({xColumnIDHistName, xColumnIDHistName, 3*xChip, 0, 3.0*xChip, 6, 0, 6}, xColumnName, yChipName);
      

    }
  }

  return {hHists1D,hHists2D};
 
}
