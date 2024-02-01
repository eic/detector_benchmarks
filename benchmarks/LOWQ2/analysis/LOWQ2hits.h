#pragma once

#include "functors.h"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using H3ResultPtr = ROOT::RDF::RResultPtr<TH3D>;
using RVecI       = ROOT::VecOps::RVec<int>;
using RVecD       = ROOT::VecOps::RVec<double>;

// Lazy execution methods
std::tuple<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>,std::map<TString,H3ResultPtr>> createHitPlots( RNode d2, double eventWeight ){

  int xChip = 448;
  int yChip = 512;
  int xChips = 6;
  int yChips = 6;
  std::map<int,int> xPixelMin = {{1,-xChip*xChips/2},{2,-xChip*xChips/2}};
  std::map<int,int> xPixelMax = {{1, xChip*xChips/2},{2, xChip*xChips/2}};

  std::map<int,int> yPixelMin = {{1,-yChip*yChips/2},{2,-yChip*yChips/2}};
  std::map<int,int> yPixelMax = {{1, yChip*yChips/2},{2, yChip*yChips/2}};


  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
  std::map<TString,H3ResultPtr> hHists3D;
  
  d2 = d2.Define("boardID",     [xChip](RVecI xID){ return (xID + 3*xChip) / (2 * xChip); },{"xID"})
         .Define("xChipID",     [xChip](RVecI xID){ return (xID + 3*xChip) / (xChip); },{"xID"})
         .Define("xColumnID",   [xChip](RVecI xID){ return (xID + 3*xChip) / 2; },{"xID"})
         .Define("yChipID",     [yChip](RVecI yID){ return (yID + 2.75*yChip) / (yChip); },{"yID"})
         .Define("yHalfChipID", [yChip](RVecI yID){ return (yID + 2.75*yChip) / (0.5*yChip); },{"yID"})
         .Define("eventWeight", [eventWeight](){return eventWeight;},{});


  // Plot number of hits per moduleID
  hHists1D["hmoduleID"] = d2.Histo1D({"hmoduleID",  "hmoduleID;Module ID;Rate [Hz]", 2,     1,    3   }, "moduleID", "eventWeight");
  
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
    

    auto d3 = d2.Define("ModuleFilter",[i](RVecI modID){return modID==i;},{"moduleID"})
      .Define(layerName,"layerID[ModuleFilter]");
        
    
    hHists1D[moduleTag+"/"+layerHistName]  = d3.Histo1D({layerHistName,   layerHistName+";Layer ID;Rate [Hz]",  4,     0,    4   }, layerName, "eventWeight");
    
    // Layer Loop
    for(int j=0; j<=3; j++){
      
      TString layerNum      = std::to_string(j);
      TString layerTag      = "layer"+layerNum;
      
      TString xName         = "xPixel";
      TString xHistName     = "h"+xName;

      TString yName         = "yPixel";
      TString yHistName     = "h"+yName;

      TString xChipName     = "xChip";
      TString yChipName     = "yChip";

      TString yHalfChipName = "yHalfChip";

      TString xyHistName    = "h"+xName+yName;

      TString xColumnName    = "xColumn";

      TString boardName     = "board";

      std::vector<string>  layerSizeInput = {xName.Data()};
      TString layerSizeName = "HitsPerEvent";
      TString sizeHistName  = "h"+layerSizeName;
      
      auto d4 = d3.Define("LayerFilter",[j](RVecI layID){return layID==j;},{"layerID"})
                  .Define(xName,"xID[LayerFilter&&ModuleFilter]")
                  .Define(yName,"yID[LayerFilter&&ModuleFilter]")
                  .Define(boardName,"boardID[LayerFilter&&ModuleFilter]")
                  .Define(xChipName,"xChipID[LayerFilter&&ModuleFilter]")
                  .Define(yChipName,"yChipID[LayerFilter&&ModuleFilter]")
                  .Define(yHalfChipName,"yHalfChipID[LayerFilter&&ModuleFilter]")
                  .Define(xColumnName,"xColumnID[LayerFilter&&ModuleFilter]")
                  .Define(layerSizeName,[](RVecI lay){return lay.size();},layerSizeInput);

      
    
      hHists1D[moduleTag+"/"+layerTag+"/"+xHistName]    = d4.Histo1D({xHistName,   xHistName+";x pixel column;Rate [Hz]",  xRange, xMin, xMax   }, xName, "eventWeight");
      
      hHists1D[moduleTag+"/"+layerTag+"/"+yHistName]    = d4.Histo1D({yHistName,   yHistName+";y pixel column;Rate [Hz]",  yRange, yMin, yMax   }, yName, "eventWeight");
            
      hHists2D[moduleTag+"/"+layerTag+"/"+xyHistName]   = d4.Histo2D({xyHistName,  xyHistName+";x pixel;y pixel;Rate [Hz]", xRange, xMin, xMax, yRange, yMin, yMax}, xName, yName, "eventWeight");
      
      hHists1D[moduleTag+"/"+layerTag+"/"+sizeHistName] = d4.Histo1D({sizeHistName,   sizeHistName+";hits per event",  100, 0, 100   }, layerSizeName );
      
      // Plot number of hits per boardID for each layer
      TString boardIDHistName = "hBoardID";
      hHists1D[moduleTag+"/"+layerTag+"/"+boardIDHistName] = d4.Histo1D({boardIDHistName, boardIDHistName+";board ID;Rate [Hz]", 3, 0, 3}, boardName, "eventWeight");
      
      // Plot number of hits per chipID for each layer
      TString chipIDHistName = "hChipID";
      hHists2D[moduleTag+"/"+layerTag+"/"+chipIDHistName] = d4.Histo2D({chipIDHistName, chipIDHistName+";x Chip;y Chip;Rate [Hz]", 6, 0, 6, 6, 0, 6}, xChipName, yChipName, "eventWeight");
      
      // Plot number of hits per chipID for each layer
      TString xColumnIDHistName = "hxColumnID";
      hHists2D[moduleTag+"/"+layerTag+"/"+xColumnIDHistName] = d4.Histo2D({xColumnIDHistName, xColumnIDHistName+";x Column;y Half Chip;Rate [Hz]", 3*xChip, 0, 3.0*xChip, 12, 0, 12}, xColumnName, yHalfChipName, "eventWeight");
      
    }
  }


  return {hHists1D,hHists2D,hHists3D};
 
}
