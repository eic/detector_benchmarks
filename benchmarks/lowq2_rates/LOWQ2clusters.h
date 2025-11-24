// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using H3ResultPtr = ROOT::RDF::RResultPtr<TH3D>;
using RVecI       = ROOT::VecOps::RVec<int>;

// Lazy execution methods
std::tuple<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>,std::map<TString,H3ResultPtr>> createClusterPlots( RNode d1 ){

  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
  std::map<TString,H3ResultPtr> hHists3D;
  
  //  auto d2 = d1.Define("ClusterSize",[](RVecI clu){return clu.size();},{"TaggerTrackerClusterPositions_0.index"});
  auto d2 = d1.Define("ClusterSize","TaggerTrackerClusterPositions.rawHits_end-TaggerTrackerClusterPositions.rawHits_begin");
  hHists1D["hClusterSize"] = d2.Histo1D({"hClusterSize","hClusterSize",100,0,100}, "ClusterSize"); 

  //Need associations back to simulated hit positions to get hit resolutions... tricky...

  return {hHists1D,hHists2D,hHists3D};
 
}
