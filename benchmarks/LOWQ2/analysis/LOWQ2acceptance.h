#pragma once

// #include "functors.h"

  #include "ROOT/RDF/RInterface.hxx"
  #include <TH1D.h>
  #include <TH2D.h>
  #include "ROOT/RVec.hxx"

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using RVecI       = ROOT::VecOps::RVec<int>;
using RVecD       = ROOT::VecOps::RVec<double>;
using RVecS       = ROOT::VecOps::RVec<string>;

// Lazy execution methods
std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>> createAcceptancePlots( RNode d1 ){

  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
  
  int ePrimID = 4;
  int pBeamID = 1;
  int eBeamID = 0;

  // d1 = d1.Define("primMom","MCParticles[2].momentum.z");
  d1 = d1.Define("eBeam", [eBeamID]
                            (std::vector<edm4hep::MCParticleData> h)
                            {return ROOT::Math::PxPyPzMVector (h[eBeamID].momentum.x,h[eBeamID].momentum.y,h[eBeamID].momentum.z,h[eBeamID].mass);},
                            {"MCParticles"})
        .Define("pBeam", [pBeamID]
                            (std::vector<edm4hep::MCParticleData> h)
                            {return ROOT::Math::PxPyPzMVector (h[pBeamID].momentum.x,h[pBeamID].momentum.y,h[pBeamID].momentum.z,h[pBeamID].mass);},
                            {"MCParticles"})
        .Define("primMom", [ePrimID]
                            (std::vector<edm4hep::MCParticleData> h)
                            {return ROOT::Math::PxPyPzMVector (h[ePrimID].momentum.x,h[ePrimID].momentum.y,h[ePrimID].momentum.z,h[ePrimID].mass);},
                            {"MCParticles"})
        .Define("primE","primMom.E()")
        .Define("primTheta","1000*(TMath::Pi()-primMom.Theta())")
        .Define("primEta","primMom.Eta()")
        .Define("scatteredV","eBeam-primMom")
        .Define("Q2","-scatteredV.M2()")
        .Define("logQ2","log10(Q2)")
        .Define("x","Q2/(2*scatteredV.Dot(pBeam))")
        .Define("logx","log10(x)")
        .Define("primPhi","primMom.Phi()")
        .Define("W2","(pBeam+eBeam-primMom).M2()");
        //  .Define("primTheta",[](std::vector<edm4hep::MCParticle> part){return part[4].getMomentum().getTheta()},"{MCParticles}");

  RVecS colNames = d1.GetColumnNames();
  std::map<TString, std::string_view> filters;
  filters["All"]  = "true";
  filters["ThetaCut"] = "primTheta<11";
  if(Any(colNames=="TaggerTrackerHits")){
    filters["Any"]  = "TaggerTrackerHits.size() > 0";
    filters["Mod1"] = "moduleID[moduleID==1].size()>0";
    filters["Mod2"] = "moduleID[moduleID==2].size()>0";
  }

  if(Any(colNames=="LowQ2Tracks")){
    filters["Reconstructed"] = "LowQ2Tracks.size()>0";
  }

  //Store the counts for each filter
  RVecI filterCounts;
  int rawCount = 0;
  TString filtersName = "";

  for (auto& filter : filters)
  {
    RNode filterNode = d1.Filter(filter.second);
    TString rawString = filter.first;
    filtersName += rawString + "_";
    TString energyHistName = "h" + rawString + "Energy";
    TString thetaHistName = "h" + rawString + "Theta";
    TString etaHistName = "h" + rawString + "Eta";
    TString logQ2HistName = "h" + rawString + "logQ2";
    TString logxHistName = "h" + rawString + "logx";
    TString Q2HistName = "h" + rawString + "Q2";
    TString xHistName = "h" + rawString + "x";
    TString W2HistName = "h" + rawString + "W2";
    TString Q2xHistName = "h" + rawString + "Q2x";
    TString logQ2logxHistName = "h" + rawString + "logQ2logx";
    TString WlogQ2HistName = "h" + rawString + "WlogQ2";
    TString WlogxHistName = "h" + rawString + "Wlogx";
    TString primElogQ2HistName = "h" + rawString + "primElogQ2";

      hHists1D[rawString+"/"+energyHistName] = filterNode.Histo1D({energyHistName, energyHistName + ";Energy [GeV]", 100, 0, 18}, "primE");
      hHists1D[rawString+"/"+thetaHistName]  = filterNode.Histo1D({thetaHistName, thetaHistName + ";Theta [mrad]", 100, 0, 18}, "primTheta");
      hHists1D[rawString+"/"+etaHistName] = filterNode.Histo1D({etaHistName, etaHistName + ";Eta", 100, -15, 5}, "primEta");
      hHists1D[rawString+"/"+logQ2HistName] = filterNode.Histo1D({logQ2HistName, logQ2HistName + ";log(Q2) [GeV]", 100, -10, 4}, "logQ2");
      hHists1D[rawString+"/"+logxHistName] = filterNode.Histo1D({logxHistName, logxHistName + ";log(x) [GeV]", 100, -12, 0}, "logx");
      hHists1D[rawString+"/"+Q2HistName] = filterNode.Histo1D({Q2HistName, Q2HistName + ";Q2 [GeV]", 100, 0, 500}, "Q2");
      hHists1D[rawString+"/"+xHistName] = filterNode.Histo1D({xHistName, xHistName + ";x [GeV]", 100, 0, 1}, "x");
      hHists1D[rawString+"/"+W2HistName] = filterNode.Histo1D({W2HistName, W2HistName + ";W2 [GeV]", 100, 0, 500}, "W2");

      hHists2D[rawString+"/"+Q2xHistName] = filterNode.Histo2D({Q2xHistName, Q2xHistName + ";Q2 [GeV];x [GeV]", 100, 0, 500, 100, 0, 1}, "Q2", "x");
      hHists2D[rawString+"/"+logQ2logxHistName] = filterNode.Histo2D({logQ2logxHistName, logQ2logxHistName + ";log(Q2) [GeV];log(x) [GeV]", 100, -10, 4, 100, -12, 0}, "logQ2", "logx");
      hHists2D[rawString+"/"+WlogQ2HistName] = filterNode.Histo2D({WlogQ2HistName, WlogQ2HistName + ";W2 [GeV];log(Q2) [GeV]", 100, 0, 500, 100, -10, 4}, "W2", "logQ2");
      hHists2D[rawString+"/"+WlogxHistName] = filterNode.Histo2D({WlogxHistName, WlogxHistName + ";W2 [GeV];log(x) [GeV]", 100, 0, 500, 100, -12, 0}, "W2", "logx");
      hHists2D[rawString+"/"+primElogQ2HistName] = filterNode.Histo2D({primElogQ2HistName, primElogQ2HistName + ";E [GeV];log(Q2) [GeV]", 100, 0, 18, 100, -10, 4}, "primE", "logQ2");
  
  
    int counts = *filterNode.Count();
    filterCounts.push_back(counts);    
    if (filter.first == "All") {
      rawCount = counts;
    }
  }

  filtersName = filtersName(0,filtersName.Length()-1);

  // Number of filters
  int Nfilters = filterCounts.size();
  RVecI entryInt = ROOT::VecOps::Range(Nfilters);

  // Fraction of events passing each filter
  RVecD fractions = filterCounts/static_cast<double>(rawCount);
  
  ROOT::RDataFrame df(Nfilters);
  auto dfWithCounts = df.Define("entry", [entryInt](ULong64_t entry) { return (double)entryInt[(int)entry]; }, {"rdfentry_"})
       .Define("count", [filterCounts](ULong64_t entry) { return (double)filterCounts[(int)entry]; }, {"rdfentry_"})
       .Define("fraction", [fractions](ULong64_t entry) { return (double)fractions[(int)entry]; }, {"rdfentry_"});

  // Use the new column in the histograms
  hHists1D["TotalCounts"] = dfWithCounts.Histo1D({"hTotalCounts", "Total Counts;"+filtersName+";Count", Nfilters, 0,  static_cast<double>(Nfilters)}, "entry","count");
  hHists1D["Fractions"]   = dfWithCounts.Histo1D({"hFractions", "Fractions;"+filtersName+";Fraction", Nfilters, 0,  static_cast<double>(Nfilters)}, "entry","fraction");
  hHists1D["TotalCounts"]->Sumw2(false);
  hHists1D["Fractions"]->Sumw2(false);

  return {hHists1D,hHists2D};
}