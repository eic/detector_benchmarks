#pragma once

#include "functors.h"
#include "DD4hep/Detector.h"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using RVecI       = ROOT::VecOps::RVec<int>;

// Lazy execution methods
std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>> createHitPlots( RNode d1){

  int xChip = 448;
  int yChip = 512;
  std::map<int,int> xPixelMin = {{1,-xChip*3},{2,-xChip*3}};
  std::map<int,int> xPixelMax = {{1, xChip*3},{2, xChip*3}};

  std::map<int,int> yPixelMin = {{1,-yChip*3},{2,-yChip*3}};
  std::map<int,int> yPixelMax = {{1, yChip*3},{2, yChip*3}};


  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
    

  d1 = d1.Define("boardID",   [xChip](RVecI xID){ return (xID + 3*xChip) / (2 * xChip); },{"xID"})
         .Define("xChipID",   [xChip](RVecI xID){ return (xID + 3*xChip) / (xChip); },{"xID"})
         .Define("xColumnID", [xChip](RVecI xID){ return (xID + 3*xChip) / 2; },{"xID"})
         .Define("yChipID",   [yChip](RVecI yID){ return (yID + 2.75*yChip) / (yChip); },{"yID"});

  hHists1D["hmoduleID"] = d1.Histo1D({"hmoduleID",  "hmoduleID;Module ID", 2,     1,    3   }, "moduleID");

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
        
    hHists1D[moduleTag+"/"+layerHistName]  = d2.Histo1D({layerHistName,   layerHistName+";Layer ID",  4,     0,    4   }, layerName );
    
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

      
    
      hHists1D[moduleTag+"/"+layerTag+"/"+xHistName]    = d3.Histo1D({xHistName,   xHistName+";x pixel column",  xRange, xMin, xMax   }, xName );
      hHists1D[moduleTag+"/"+layerTag+"/"+yHistName]    = d3.Histo1D({yHistName,   yHistName+";y pixel column",  yRange, yMin, yMax   }, yName );
            
      hHists2D[moduleTag+"/"+layerTag+"/"+xyHistName]   = d3.Histo2D({xyHistName,  xyHistName+";x pixel;y pixel", xRange, xMin, xMax, yRange, yMin, yMax}, xName, yName );
      
      hHists1D[moduleTag+"/"+layerTag+"/"+sizeHistName] = d3.Histo1D({sizeHistName,   sizeHistName+";hits per event",  100, 0, 100   }, layerSizeName );

      // Plot number of hits per boardID for each layer
      TString boardIDHistName = "hBoardID" +moduleTag + layerTag;
      hHists1D[moduleTag+"/"+layerTag+"/"+boardIDHistName] = d3.Histo1D({boardIDHistName, boardIDHistName+";board ID", 3, 0, 3}, boardName);
      
      // Plot number of hits per chipID for each layer
      TString chipIDHistName = "hChipID" +moduleTag + layerTag;
      hHists2D[moduleTag+"/"+layerTag+"/"+chipIDHistName] = d3.Histo2D({chipIDHistName, chipIDHistName+";x Chip;y Chip", 6, 0, 6, 6, 0, 6}, xChipName, yChipName);

      // Plot number of hits per chipID for each layer
      TString xColumnIDHistName = "hxColumnID" +moduleTag + layerTag;
      hHists2D[moduleTag+"/"+layerTag+"/"+xColumnIDHistName] = d3.Histo2D({xColumnIDHistName, xColumnIDHistName+";x Column;y Chip", 3*xChip, 0, 3.0*xChip, 6, 0, 6}, xColumnName, yChipName);
      

    }
  }

  return {hHists1D,hHists2D};
 
}
