#pragma once

#include "functors.h"
#include <DD4hep/Detector.h>
#include <DD4hep/DetElement.h>

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using H3ResultPtr = ROOT::RDF::RResultPtr<TH3D>;
using RVecI       = ROOT::VecOps::RVec<int>;
using RVecD       = ROOT::VecOps::RVec<double>;

// Structure to hold tagger dimensions
struct TaggerDimensions {
    double width_mm;
    double height_mm;
    double pixel_size_mm;
    int n_pixels_x;
    int n_pixels_y;
};

// Function to extract tagger dimensions from detector description
std::map<int, TaggerDimensions> getTaggerDimensions(const dd4hep::Detector& detector) {
    std::map<int, TaggerDimensions> dimensions;
    
    TaggerDimensions tagger1, tagger2;
    
    try {
        // Try to extract tagger dimensions from detector constants
        std::map<std::string, double> constants;
        
        // Attempt to access detector constants
        auto constantsMap = detector.constants();
        
        // Look for tagger dimension constants and pixel size
        bool foundTagger1Width = false, foundTagger1Height = false;
        bool foundTagger2Width = false, foundTagger2Height = false;
        bool foundPixelSize = false;
        double pixel_size_mm = 0.055; // Default 55 μm pixel size

        pixel_size_mm = detector.constant<double>("tracker_pixel_size") * 10.0; // Convert cm to mm
        tagger1.width_mm = detector.constant<double>("Tagger1_Width") * 10.0; // Convert cm to mm
        tagger1.height_mm = detector.constant<double>("Tagger1_Height") * 10.0; // Convert cm to mm
        tagger2.width_mm = detector.constant<double>("Tagger2_Width") * 10.0; // Convert cm to mm
        tagger2.height_mm = detector.constant<double>("Tagger2_Height") * 10.0; // Convert cm to mm
        
        // Set pixel size for both taggers
        tagger1.pixel_size_mm = pixel_size_mm;
        tagger2.pixel_size_mm = pixel_size_mm;
        
        if (!foundPixelSize) {
            std::cout << "Using default pixel size: " << pixel_size_mm << " mm" << std::endl;
        } else {
            std::cout << "Loaded pixel size from detector: " << pixel_size_mm << " mm" << std::endl;
        }
        
        // Calculate pixel counts from dimensions if found, otherwise use defaults
        // if (foundTagger1Width && foundTagger1Height) {
            tagger1.n_pixels_x = static_cast<int>(tagger1.width_mm / pixel_size_mm);
            tagger1.n_pixels_y = static_cast<int>(tagger1.height_mm / pixel_size_mm);
        // } else {
        //     // Default values based on current pixel layout
        //     tagger1.n_pixels_x = 6 * 448;  // 6 chips × 448 pixels
        //     tagger1.n_pixels_y = 6 * 512;  // 6 chips × 512 pixels
        //     tagger1.width_mm = tagger1.n_pixels_x * tagger1.pixel_size_mm;
        //     tagger1.height_mm = tagger1.n_pixels_y * tagger1.pixel_size_mm;
        //     std::cout << "Using default dimensions for Tagger1" << std::endl;
        // }
        
        // if (foundTagger2Width && foundTagger2Height) {
            tagger2.n_pixels_x = static_cast<int>(tagger2.width_mm / pixel_size_mm);
            tagger2.n_pixels_y = static_cast<int>(tagger2.height_mm / pixel_size_mm);
        // } else {
        //     // Default values (same as tagger1)
        //     tagger2.n_pixels_x = tagger1.n_pixels_x;
        //     tagger2.n_pixels_y = tagger1.n_pixels_y;
        //     tagger2.width_mm = tagger1.width_mm;
        //     tagger2.height_mm = tagger1.height_mm;
        //     std::cout << "Using default dimensions for Tagger2" << std::endl;
        // }
        
        dimensions[1] = tagger1;
        dimensions[2] = tagger2;
        
        std::cout << "Loaded tagger dimensions from detector description:" << std::endl;
        std::cout << "  Tagger1: " << tagger1.width_mm << " mm × " << tagger1.height_mm << " mm" << std::endl;
        std::cout << "  Tagger2: " << tagger2.width_mm << " mm × " << tagger2.height_mm << " mm" << std::endl;
        
    } catch (...) {
        // Fallback to default values
        std::cout << "Warning: Could not extract tagger dimensions from XML, using defaults" << std::endl;
        
        tagger1.pixel_size_mm = 0.055;
        tagger1.n_pixels_x = 6 * 448;
        tagger1.n_pixels_y = 6 * 512;
        tagger1.width_mm = tagger1.n_pixels_x * tagger1.pixel_size_mm;
        tagger1.height_mm = tagger1.n_pixels_y * tagger1.pixel_size_mm;
        
        tagger2 = tagger1;
        dimensions[1] = tagger1;
        dimensions[2] = tagger2;
    }
    
    return dimensions;
}

// Lazy execution methods
std::tuple<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>,std::map<TString,H3ResultPtr>> createHitPlots( RNode d2, double eventWeight, const dd4hep::Detector& detector ){

  // Get tagger dimensions from detector description
  auto taggerDims = getTaggerDimensions(detector);
  
  // Use the detector dimensions to set pixel ranges
  std::map<int,int> xPixelMin, xPixelMax, yPixelMin, yPixelMax;
  
  for (auto& [module, dims] : taggerDims) {
      int half_x = dims.n_pixels_x / 2;
      int half_y = dims.n_pixels_y / 2;
      
      xPixelMin[module] = -half_x;
      xPixelMax[module] = half_x;
      yPixelMin[module] = -half_y;
      yPixelMax[module] = half_y;
      
      std::cout << "Module " << module << ": " << dims.width_mm << " mm × " << dims.height_mm 
                << " mm (" << dims.n_pixels_x << " × " << dims.n_pixels_y << " pixels)" << std::endl;
  }

  // Get chip dimensions (assuming consistent across modules)
  int xChip = 448;
  int yChip = 512;
  
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
  hHists1D["hmoduleIDRate"] = d2.Histo1D({"hmoduleIDRate",  "hmoduleID;Module ID;Rate [Hz]", 2,     1,    3   }, "moduleID", "eventWeight");
  
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
    TString layerHistName  = "h"+layerName+"Rate";
    

    auto d3 = d2.Define("ModuleFilter",[i](RVecI modID){return modID==i;},{"moduleID"})
      .Define(layerName,"layerID[ModuleFilter]");
        
    
    hHists1D[moduleTag+"/"+layerHistName]  = d3.Histo1D({layerHistName,   layerHistName+";Layer ID;Rate [Hz]",  4,     0,    4   }, layerName, "eventWeight");
    
    // Layer Loop
    for(int j=0; j<=3; j++){
      
      TString layerNum      = std::to_string(j);
      TString layerTag      = "layer"+layerNum;
      
      TString xName         = "xPixel";
      TString xHistName     = "h"+xName+"Rate";

      TString yName         = "yPixel";
      TString yHistName     = "h"+yName+"Rate";

      TString xChipName     = "xChip";
      TString yChipName     = "yChip";

      TString yHalfChipName = "yHalfChip";

      TString xyHistName    = "h"+xName+yName+"Rate";

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
      TString boardIDHistName = "hBoardIDRate";
      hHists1D[moduleTag+"/"+layerTag+"/"+boardIDHistName] = d4.Histo1D({boardIDHistName, boardIDHistName+";board ID;Rate [Hz]", 3, 0, 3}, boardName, "eventWeight");
      
      // Plot number of hits per chipID for each layer
      TString chipIDHistName = "hChipIDRate";
      hHists2D[moduleTag+"/"+layerTag+"/"+chipIDHistName] = d4.Histo2D({chipIDHistName, chipIDHistName+";x Chip;y Chip;Rate [Hz]", 6, 0, 6, 6, 0, 6}, xChipName, yChipName, "eventWeight");
      
      // Plot number of hits per chipID for each layer
      TString xColumnIDHistName = "hxColumnIDRate";
      hHists2D[moduleTag+"/"+layerTag+"/"+xColumnIDHistName] = d4.Histo2D({xColumnIDHistName, xColumnIDHistName+";x Column;y Half Chip;Rate [Hz]", 3*xChip, 0, 3.0*xChip, 12, 0, 12}, xColumnName, yHalfChipName, "eventWeight");
      
    }
  }

  return {hHists1D,hHists2D,hHists3D};
 
}
