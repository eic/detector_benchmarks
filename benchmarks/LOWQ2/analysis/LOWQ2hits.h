#include "functors.h"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using RVecI       = ROOT::VecOps::RVec<int>;

std::map<int,int> xPixelMin = {{1,-1344},{2,-1344}};
std::map<int,int> xPixelMax = {{1, 1344},{2, 1344}};

std::map<int,int> yPixelMin = {{1,-1818},{2,-1364}};
std::map<int,int> yPixelMax = {{1, 1818},{2, 1364}};

// Lazy execution methods
std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>> createHitPlots( RNode d1 ){

  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
  
  string compactName = "/home/simon/EIC/epic/epic_18x275.xml";
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

      TString xyHistName    = "h"+xName+yName;

      std::vector<string>  layerSizeInput = {xName.Data()};
      TString layerSizeName = "HitsPerEvent"+moduleTag+layerTag;
      TString sizeHistName  = "h"+layerSizeName;
      
      auto d3 = d2.Define("LayerFilter",[j](RVecI layID){return layID==j;},{"layerID"})
	.Define(xName,"xID[LayerFilter&&ModuleFilter]")
	.Define(yName,"yID[LayerFilter&&ModuleFilter]")
	.Define(layerSizeName,[](RVecI lay){return lay.size();},layerSizeInput);
      
    
      hHists1D[xHistName]    = d3.Histo1D({xHistName,   xHistName,  xRange, xMin, xMax   }, xName );
      hHists1D[yHistName]    = d3.Histo1D({yHistName,   yHistName,  yRange, yMin, yMax   }, yName );
            
      hHists2D[xyHistName]   = d3.Histo2D({xyHistName,  xyHistName, xRange, xMin, xMax, yRange, yMin, yMax}, xName, yName );
      
      hHists1D[sizeHistName] = d3.Histo1D({sizeHistName,   sizeHistName,  100, 0, 100   }, layerSizeName );
      
    }

  }

  return {hHists1D,hHists2D};
 
}
