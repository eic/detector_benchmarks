#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "DD4hep/VolumeManager.h"
#include "TFile.h"

using RVecHits = ROOT::VecOps::RVec<edm4hep::SimTrackerHitData>;

//-----------------------------------------------------------------------------------------
// Grab Component functor
//-----------------------------------------------------------------------------------------
struct getSubID{
  getSubID(std::string cname, dd4hep::Detector& det, std::string rname = "BackwardsBeamlineHits") : componentName(cname), detector(det), readoutName(rname){}
  
  ROOT::VecOps::RVec<int> operator()(const RVecHits& evt) {
    auto decoder = detector.readout(readoutName).idSpec().decoder();
    auto indexID = decoder->index(componentName);
    ROOT::VecOps::RVec<int> result;
    for(const auto& h: evt) {
      result.push_back(decoder->get(h.cellID,indexID));      
    }
    return result;    
  };
  
  void SetComponent(std::string cname){
    componentName = cname;
  }
  void SetReadout(std::string rname){
    readoutName = rname;
  }
  
  private: 
  std::string componentName;
  dd4hep::Detector& detector;
  std::string readoutName;
};

//-----------------------------------------------------------------------------------------
// Transform global x,y,z position and momentum into local coordinates
//-----------------------------------------------------------------------------------------
struct globalToLocal{
  globalToLocal(dd4hep::Detector& det) : detector(det){
     volumeManager = dd4hep::VolumeManager::getVolumeManager(det);
  }
  
  ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> operator()(const RVecHits& evt) {

    ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> result;
    ROOT::VecOps::RVec<double> xvec;
    ROOT::VecOps::RVec<double> yvec;
    ROOT::VecOps::RVec<double> zvec;
    ROOT::VecOps::RVec<double> pxvec;
    ROOT::VecOps::RVec<double> pyvec;
    ROOT::VecOps::RVec<double> pzvec;

    for(const auto& h: evt) {
      auto cellID = h.cellID;
      // dd4hep::DetElement detelement = volumeManager.lookupDetElement(cellID);
      // auto detelement = volumeManager.lookupVolumePlacement(cellID);
      auto detelement = volumeManager.lookupDetElement(cellID);
      const TGeoMatrix& transform = detelement.nominal().worldTransformation();
      // transform.Print();
      //auto transform = volumeManager.worldTransformation(cellID);
      //Convert position to local coordinates
      auto pos = h.position;
      Double_t globalPos[3] = {pos.x/10, pos.y/10, pos.z/10};
      Double_t localPos[3];
      transform.MasterToLocal(globalPos, localPos);
      // std::cout << "Global: " << globalPos[0] << " " << globalPos[1] << " " << globalPos[2] << std::endl;
      // std::cout << "Local: " << localPos[0] << " " << localPos[1] << " " << localPos[2] << std::endl;
      //Transform global momentum to local coordinates
      auto mom = h.momentum;
      Double_t globalMom[3] = {mom.x, mom.y, mom.z};
      Double_t localMom[3];
      transform.MasterToLocalVect(globalMom, localMom);
      // std::cout << "Global: " << globalMom[0] << " " << globalMom[1] << " " << globalMom[2] << std::endl;
      // std::cout << "Local: " << localMom[0] << " " << localMom[1] << " " << localMom[2] << std::endl;

      xvec.push_back(localPos[0]);
      yvec.push_back(localPos[1]);
      zvec.push_back(localPos[2]);
      pxvec.push_back(localMom[0]);
      pyvec.push_back(localMom[1]);
      pzvec.push_back(localMom[2]);
      
    }

    result.push_back(xvec);
    result.push_back(yvec);
    result.push_back(zvec);
    result.push_back(pxvec);
    result.push_back(pyvec);
    result.push_back(pzvec);    

    return result;    
  };
  
  private: 
  dd4hep::Detector& detector;
  dd4hep::VolumeManager volumeManager;
};

struct volParams{
  double radius;
  double xPos;
  double yPos;
  double zPos;
  double rotation;
  bool isConeSegment;
};

// Functor to get the volume element parameters from the CellID
struct getVolumeParametersFromCellID {
  getVolumeParametersFromCellID(dd4hep::Detector& det) : detector(det) {
    volumeManager = dd4hep::VolumeManager::getVolumeManager(det);
  }

  ROOT::VecOps::RVec<volParams> operator()(const RVecHits& evt) {
    ROOT::VecOps::RVec<volParams> result;
    // Look up the detector element using the CellID
    for(const auto& h : evt) {
      auto  cellID = h.cellID;
      auto  detelement = volumeManager.lookupDetElement(cellID);
      const TGeoMatrix& transform      = detelement.nominal().worldTransformation();
      const Double_t*   translation    = transform.GetTranslation();
      const Double_t*   rotationMatrix = transform.GetRotationMatrix(); // Compute rotation angle around the Y-axis
      double rotationAngleY = std::atan2(rotationMatrix[2], rotationMatrix[8]); // atan2(r13, r33)
      auto volume = detelement.volume();
      dd4hep::Solid solid = volume.solid();
      bool isCone = solid.isValid() && std::string(solid->IsA()->GetName()) == "TGeoConeSeg";
      double radius = 0.0;
      if (isCone) {
        dd4hep::ConeSegment cone = solid;
        radius = cone.rMax1();
      }
      volParams params{
        radius,
        translation[0],
        translation[1],
        translation[2],
        rotationAngleY,
        isCone
      };
      result.push_back(params);
    }
    return result;
  }

private:
  dd4hep::Detector& detector;
  dd4hep::VolumeManager volumeManager;
};

TH1F* CreateFittedHistogram(const std::string& histName, 
  const std::string& histTitle, 
  const std::map<TString, double>& valueMap, 
  const std::map<TString, double>& errorMap, 
  const std::string& xAxisLabel) {
  // Number of bins corresponds to the number of entries in the value map
  int nPoints = valueMap.size();
  TH1F* hist = new TH1F(histName.c_str(), (";" + xAxisLabel + ";" + histTitle).c_str(), nPoints, 0.5, nPoints + 0.5);

  // Fill the histogram with values and errors, and set custom bin labels
  int binIndex = 1; // Start from bin 1
  for (const auto& [name, value] : valueMap) {
    hist->SetBinContent(binIndex, value); // Set the bin content
    if(errorMap.size()==valueMap.size()){
      hist->SetBinError(binIndex, errorMap.at(name)); // Set the bin error
    }
    else{
      hist->SetBinError(binIndex, 0); // Set the bin error to 0 if not provided
    }
    hist->GetXaxis()->SetBinLabel(binIndex, name); // Set the bin label
    binIndex++;
  }

  // Customize the histogram
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.0);
  hist->SetLineColor(kBlue);
  hist->SetMarkerColor(kRed);

  return hist;
}