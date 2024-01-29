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
  
  int ePrimID = 2;
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

  

    //Histogram names
    TString energyHistName        = "h" + rawString + "Energy";
    TString thetaHistName         = "h" + rawString + "Theta";
    TString etaHistName           = "h" + rawString + "Eta";
    TString logQ2HistName         = "h" + rawString + "logQ2";
    TString logxHistName          = "h" + rawString + "logx";
    TString Q2HistName            = "h" + rawString + "Q2";
    TString xHistName             = "h" + rawString + "x";
    TString W2HistName            = "h" + rawString + "W2";
    TString Q2xHistName           = "h" + rawString + "Q2x";
    TString logQ2logxHistName     = "h" + rawString + "logQ2logx";
    TString WlogQ2HistName        = "h" + rawString + "WlogQ2";
    TString WlogxHistName         = "h" + rawString + "Wlogx";
    TString primElogQ2HistName    = "h" + rawString + "primElogQ2";
    TString primEthetaHistName    = "h" + rawString + "primEtheta";

    //Histogram Ranges
    double energyMin = 0;
    double energyMax = 18;
    double thetaMin  = 0;
    double thetaMax  = 12;
    double etaMin    = -15;
    double etaMax    = 5;
    double logQ2Min  = -8.5;
    double logQ2Max  = 0;
    double logxMin   = -12;
    double logxMax   = 0;
    double Q2Min     = 0;
    double Q2Max     = 500;
    double xMin      = 0;
    double xMax      = 1;
    double W2Min     = 0;
    double W2Max     = 500;

    int nBins1D = 400;
    int nBins2D = 200;

    //Histogram axis names
    TString energyAxisName        = "Energy [GeV]";
    TString thetaAxisName         = "Theta [mrad]";
    TString etaAxisName           = "Eta";
    TString logQ2AxisName         = "log(Q2) [GeV]";
    TString logxAxisName          = "log(x) [GeV]";
    TString Q2AxisName            = "Q2 [GeV]";
    TString xAxisName             = "x [GeV]";
    TString W2AxisName            = "W2 [GeV]";
    TString Q2xAxisName           = Q2AxisName+";"+xAxisName;
    TString logQ2logxAxisName     = logQ2AxisName+";"+logxAxisName;
    TString WlogQ2AxisName        = W2AxisName+";"+logQ2AxisName;
    TString WlogxAxisName         = W2AxisName+";"+logxAxisName;
    TString primElogQ2AxisName    = energyAxisName+";"+logQ2AxisName;
    TString primEthetaAxisName    = energyAxisName+";"+thetaAxisName;


    hHists1D[rawString+"/"+energyHistName] = filterNode.Histo1D({energyHistName, energyHistName + ";Energy [GeV]", nBins1D, energyMin, energyMax}, "primE");
    hHists1D[rawString+"/"+thetaHistName]  = filterNode.Histo1D({thetaHistName, thetaHistName + ";Theta [mrad]", nBins1D, thetaMin, thetaMax}, "primTheta");
    hHists1D[rawString+"/"+etaHistName]    = filterNode.Histo1D({etaHistName, etaHistName + ";Eta", nBins1D, etaMin, etaMax}, "primEta");
    hHists1D[rawString+"/"+logQ2HistName]  = filterNode.Histo1D({logQ2HistName, logQ2HistName + ";log(Q2) [GeV]", nBins1D, logQ2Min, logQ2Max}, "logQ2");
    hHists1D[rawString+"/"+logxHistName]   = filterNode.Histo1D({logxHistName, logxHistName + ";log(x) [GeV]", nBins1D, logxMin, logxMax}, "logx");
    hHists1D[rawString+"/"+Q2HistName]     = filterNode.Histo1D({Q2HistName, Q2HistName + ";Q2 [GeV]", nBins1D, Q2Min, Q2Max}, "Q2");
    hHists1D[rawString+"/"+xHistName]      = filterNode.Histo1D({xHistName, xHistName + ";x [GeV]", nBins1D, xMin, xMax}, "x");
    hHists1D[rawString+"/"+W2HistName]     = filterNode.Histo1D({W2HistName, W2HistName + ";W2 [GeV]", nBins1D, W2Min, W2Max}, "W2");

    hHists2D[rawString+"/"+Q2xHistName]        = filterNode.Histo2D({Q2xHistName, Q2xHistName + ";Q2 [GeV];x [GeV]", nBins2D, Q2Min, Q2Max, nBins2D, xMin, xMax}, "Q2", "x");
    hHists2D[rawString+"/"+logQ2logxHistName]  = filterNode.Histo2D({logQ2logxHistName, logQ2logxHistName + ";log(Q2) [GeV];log(x) [GeV]", nBins2D, logQ2Min, logQ2Max, nBins2D, logxMin, logxMax}, "logQ2", "logx");
    hHists2D[rawString+"/"+WlogQ2HistName]     = filterNode.Histo2D({WlogQ2HistName, WlogQ2HistName + ";W2 [GeV];log(Q2) [GeV]", nBins2D, W2Min, W2Max, nBins2D, logQ2Min, logQ2Max}, "W2", "logQ2");
    hHists2D[rawString+"/"+WlogxHistName]      = filterNode.Histo2D({WlogxHistName, WlogxHistName + ";W2 [GeV];log(x) [GeV]", nBins2D, W2Min, W2Max, nBins2D, logxMin, logxMax}, "W2", "logx");
    hHists2D[rawString+"/"+primElogQ2HistName] = filterNode.Histo2D({primElogQ2HistName, primElogQ2HistName + ";E [GeV];log(Q2) [GeV]", nBins2D, energyMin, energyMax, nBins2D, logQ2Min, logQ2Max}, "primE", "logQ2");
    hHists2D[rawString+"/"+primEthetaHistName] = filterNode.Histo2D({primEthetaHistName, primEthetaHistName + ";E [GeV];Theta [mrad]", nBins2D, energyMin, energyMax, nBins2D, thetaMin, thetaMax}, "primE", "primTheta");


    int counts = *filterNode.Count();
    filterCounts.push_back(counts);    
    if (filter.first == "All") {
      rawCount = counts;
    }
  }

  filtersName = filtersName(0,filtersName.Length()-1);

  // Number of filters
  int Nfilters   = filterCounts.size();
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