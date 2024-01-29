#pragma once

#include "functors.h"
#include "DD4hep/Detector.h"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using RVecI       = ROOT::VecOps::RVec<int>;

// Lazy execution methods
std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>> createHitPlots( RNode d1, double eventRate ){

  int xChip = 448;
  int yChip = 512;
  std::map<int,int> xPixelMin = {{1,-xChip*3},{2,-xChip*3}};
  std::map<int,int> xPixelMax = {{1, xChip*3},{2, xChip*3}};

  std::map<int,int> yPixelMin = {{1,-yChip*3},{2,-yChip*3}};
  std::map<int,int> yPixelMax = {{1, yChip*3},{2, yChip*3}};


  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
        
  int events = *d1.Count();
  double eventWeight = eventRate / events;

  d1 = d1.Define("boardID",   [xChip](RVecI xID){ return (xID + 3*xChip) / (2 * xChip); },{"xID"})
         .Define("xChipID",   [xChip](RVecI xID){ return (xID + 3*xChip) / (xChip); },{"xID"})
         .Define("xColumnID", [xChip](RVecI xID){ return (xID + 3*xChip) / 2; },{"xID"})
         .Define("yChipID",   [yChip](RVecI yID){ return (yID + 2.75*yChip) / (yChip); },{"yID"})
         .Define("yHalfChipID",   [yChip](RVecI yID){ return (yID + 2.75*yChip) / (0.5*yChip); },{"yID"})
         .Define("eventWeight", [eventWeight](){return eventWeight;},{});


  hHists1D["hits/hmoduleID"] = d1.Histo1D({"hmoduleID",  "hmoduleID;Module ID;Counts", 2,     1,    3   }, "moduleID");
  hHists1D["rate/hmoduleID"] = d1.Histo1D({"hmoduleID",  "hmoduleID;Module ID;Rate [Hz]", 2,     1,    3   }, "moduleID", "eventWeight");
  hHists1D["rate/hmoduleID"]->Sumw2(false);

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
        
    hHists1D["hits/"+moduleTag+"/"+layerHistName]  = d2.Histo1D({layerHistName,   layerHistName+";Layer ID;Counts",  4,     0,    4   }, layerName );
    hHists1D["rate/"+moduleTag+"/"+layerHistName]  = d2.Histo1D({layerHistName,   layerHistName+";Layer ID;Rate [Hz]",  4,     0,    4   }, layerName, "eventWeight");
    hHists1D["rate/"+moduleTag+"/"+layerHistName]->Sumw2(false);
    
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

      TString yHalfChipName = "yHalfChipID"+moduleTag+layerTag;

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
                  .Define(yHalfChipName,"yHalfChipID[LayerFilter&&ModuleFilter]")
                  .Define(xColumnName,"xColumnID[LayerFilter&&ModuleFilter]")
                  .Define(layerSizeName,[](RVecI lay){return lay.size();},layerSizeInput);

      
    
      hHists1D["hits/"+moduleTag+"/"+layerTag+"/"+xHistName]    = d3.Histo1D({xHistName,   xHistName+";x pixel column;Counts",  xRange, xMin, xMax   }, xName );
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+xHistName]    = d3.Histo1D({xHistName,   xHistName+";x pixel column;Rate [Hz]",  xRange, xMin, xMax   }, xName, "eventWeight");
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+xHistName]->Sumw2(false);

      hHists1D["hits/"+moduleTag+"/"+layerTag+"/"+yHistName]    = d3.Histo1D({yHistName,   yHistName+";y pixel column;Counts",  yRange, yMin, yMax   }, yName );
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+yHistName]    = d3.Histo1D({yHistName,   yHistName+";y pixel column;Rate [Hz]",  yRange, yMin, yMax   }, yName, "eventWeight");
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+yHistName]->Sumw2(false);
            
      hHists2D["hits/"+moduleTag+"/"+layerTag+"/"+xyHistName]   = d3.Histo2D({xyHistName,  xyHistName+";x pixel;y pixel;Counts", xRange, xMin, xMax, yRange, yMin, yMax}, xName, yName );
      hHists2D["rate/"+moduleTag+"/"+layerTag+"/"+xyHistName]   = d3.Histo2D({xyHistName,  xyHistName+";x pixel;y pixel;Rate [Hz]", xRange, xMin, xMax, yRange, yMin, yMax}, xName, yName, "eventWeight");
      hHists2D["rate/"+moduleTag+"/"+layerTag+"/"+xyHistName]->Sumw2(false);
      
      hHists1D["hits/"+moduleTag+"/"+layerTag+"/"+sizeHistName] = d3.Histo1D({sizeHistName,   sizeHistName+";hits per event",  100, 0, 100   }, layerSizeName );
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+sizeHistName] = d3.Histo1D({sizeHistName,   sizeHistName+";hits per event",  100, 0, 100   }, layerSizeName, "eventWeight");
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+sizeHistName]->Sumw2(false);

      // Plot number of hits per boardID for each layer
      TString boardIDHistName = "hBoardID" +moduleTag + layerTag;
      hHists1D["hits/"+moduleTag+"/"+layerTag+"/"+boardIDHistName] = d3.Histo1D({boardIDHistName, boardIDHistName+";board ID;Counts", 3, 0, 3}, boardName);
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+boardIDHistName] = d3.Histo1D({boardIDHistName, boardIDHistName+";board ID;Rate [Hz]", 3, 0, 3}, boardName, "eventWeight");
      hHists1D["rate/"+moduleTag+"/"+layerTag+"/"+boardIDHistName]->Sumw2(false);
      
      // Plot number of hits per chipID for each layer
      TString chipIDHistName = "hChipID" +moduleTag + layerTag;
      hHists2D["hits/"+moduleTag+"/"+layerTag+"/"+chipIDHistName] = d3.Histo2D({chipIDHistName, chipIDHistName+";x Chip;y Chip;Counts", 6, 0, 6, 6, 0, 6}, xChipName, yChipName);
      hHists2D["rate/"+moduleTag+"/"+layerTag+"/"+chipIDHistName] = d3.Histo2D({chipIDHistName, chipIDHistName+";x Chip;y Chip;Rate [Hz]", 6, 0, 6, 6, 0, 6}, xChipName, yChipName, "eventWeight");
      hHists2D["rate/"+moduleTag+"/"+layerTag+"/"+chipIDHistName]->Sumw2(false);

      // Plot number of hits per chipID for each layer
      TString xColumnIDHistName = "hxColumnID" +moduleTag + layerTag;
      hHists2D["hits/"+moduleTag+"/"+layerTag+"/"+xColumnIDHistName] = d3.Histo2D({xColumnIDHistName, xColumnIDHistName+";x Column;y Half Chip;Counts", 3*xChip, 0, 3.0*xChip, 12, 0, 12}, xColumnName, yHalfChipName);
      hHists2D["rate/"+moduleTag+"/"+layerTag+"/"+xColumnIDHistName] = d3.Histo2D({xColumnIDHistName, xColumnIDHistName+";x Column;y Half Chip;Rate [Hz]", 3*xChip, 0, 3.0*xChip, 12, 0, 12}, xColumnName, yHalfChipName, "eventWeight");
      hHists2D["rate/"+moduleTag+"/"+layerTag+"/"+xColumnIDHistName]->Sumw2(false);
      
    }
  }


  return {hHists1D,hHists2D};
 
}
