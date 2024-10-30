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