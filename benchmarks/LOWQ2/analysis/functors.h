#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/SimCalorimeterHit.h"

//-----------------------------------------------------------------------------------------
// Grab Component functor
//-----------------------------------------------------------------------------------------
struct getSubID{
  getSubID(std::string cname, dd4hep::Detector& det, std::string rname = "TaggerTrackerHits") : componentName(cname), detector(det), readoutName(rname){}
  
  ROOT::VecOps::RVec<int> operator()(const std::vector<edm4hep::SimTrackerHitData>& evt) {
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
