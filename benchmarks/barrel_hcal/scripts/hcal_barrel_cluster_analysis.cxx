// ----------------------------------------------------------------------------
// 'hcal_barrel_clusters_analysis.cxx'
// Derek Anderson
// 09.07.2023
//
// ePIC BHCal benchmark macro.
// ----------------------------------------------------------------------------

// c++ utilities
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>
#include <optional>
// root classes
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
// dataframe related classes
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/HistoModels.hxx>
// benchmark utilities
#include "common_bench/mt.h"
#include "common_bench/util.h"
#include "common_bench/plot.h"
#include "common_bench/benchmark.h"
// formatting utilities
#include "fmt/core.h"
#include "fmt/color.h"
// misc
#include "nlohmann/json.hpp"

// make ROOT namespaces implicit
using namespace ROOT;
using namespace ROOT::VecOps; 

// set up aliases
using TH1Def = ROOT::RDF::TH1DModel;
using TH2Def = ROOT::RDF::TH2DModel;



// save canvas ----------------------------------------------------------------

void save_canvas(TCanvas* canvas, std::optional<std::string> label = nullopt) {

  if (label.has_value()) {
    canvas -> SaveAs( fmt::format("results/{}.png", label).c_str() );
    canvas -> SaveAs( fmt::format("results/{}.pdf", label).c_str() );
  } else {
    canvas -> SaveAs( fmt::format("results/{}.png", canvas -> GetName()).c_str() );
    canvas -> SaveAs( fmt::format("results/{}.pdf", canvas -> GetName()).c_str() );
  }

}  // end 'save_canvas(TCanvas*, std::string)'



// bhcal benchmarks -----------------------------------------------------------

int hcal_barrel_clusters_analysis(std::string file) {

  // enable multithreading
  EnableImplicitMT(kNumThreads);

  // open dataframe
  RDataFrame frame("events", file);

  // make sure file isn't empty
  auto events = frame.Count();
  if (events == 0) {
    cerr << "Error: No events found!" << endl;
    assert(events > 0);
  }

  // lambdas for analysis -----------------------------------------------------

  auto getParticleEnergy = [](const RVec<int> &types, const RVec<float> &energies) {
    float    energy = -1;
    uint64_t index  = 0;
    for (const int type : types) {
      if (type == 1) {
        energy = energies.at(index);
        break;
      }
      ++index;
    }
    return energy;
  };

  auto getLeadEnergy = [](const RVec<float> &energies) {
    float lead = -1;
    for (const float energy : energies) {
      if (energy > lead) {
        lead = energy;
      }
    }
    return lead;
  };

  auto getEnergySum = [](const RVec<float> &energies) {
    float sum = 0.;
    for (const float energy : energies) {
      sum += energy;
    }
    return sum;
  };

  auto getPercentDiffVec = [](const float a, const RVec<float> &b) {
    RVec<float> diff;
    for (uint64_t index = 0; index < b.size(); index++) {
      diff.push_back((a - b.at(index)) / a);
    }
    return diff;
  };

  auto getPercentDiff = [](const float a, const float b) {
    return (a - b) / b;
  };

  auto getSumFraction = [](const float a, const float b) {
    return a / (a + b);
  };

  auto getMultiplicity = [](const RVec<float> &collection) {
    return collection.size();
  };

  // define histograms --------------------------------------------------------

  // turn on histogram errors
  TH1::SetDefaultSumw2(true);
  TH2::SetDefaultSumw2(true);

  // histogram titles
  const string eneTitle       = ";E [GeV];counts";
  const string diffTitle      = ";(E_{reco} - E_{par}) / E_{par};counts";
  const string multTitle      = ";multiplicity;counts";
  const string diffCvsPTitle  = ";E_{par} [GeV];E_{clust} [GeV];counts";
  const string diffLvsPTitle  = ";E_{par} [GeV];E_{clust}^{lead} [GeV];counts";
  const string diffSvsPTitle  = ";E_{par} [GeV];#SigmaE_{clust} [GeV];counts";
  const string sumVsFracTitle = ";#SigmaE_{BECal} / (#SigmaE_{BECal} + #SigmaE_{BHCal});#SigmaE_{BECal};counts";

  // histogram binning
  const tuple<int, double, double> eneBins  = make_tuple(200, 0.,   100.);
  const tuple<int, double, double> diffBins = make_tuple(200, -10., 10.);
  const tuple<int, double, double> fracBins = make_tuple(100, 0.,   10.);
  const tuple<int, double, double> multBins = make_tuple(100, 0.,   100.);

  // histogram definitions
  const vector<TH1Def> vecHistDefs1D = {
    TH1Def("hEnePar",         eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneHit",         eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneClust",       eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneLead",        eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hEneSum",         eneTitle.data(),  get<0>(eneBins),  get<1>(eneBins),  get<2>(eneBins)),
    TH1Def("hDiffClust",      diffTitle.data(), get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hDiffLead",       diffTitle.data(), get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hDiffSum",        diffTitle.data(), get<0>(diffBins), get<1>(diffBins), get<2>(diffBins)),
    TH1Def("hMultHit",        multTitle.data(), get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    TH1Def("hMultClust",      multTitle.data(), get<0>(multBins), get<1>(multBins), get<2>(multBins)),
    TH1Def("hMultHitInClust", multTitle.data(), get<0>(multBins), get<2>(multBins), get<2>(multBins))
  };
  const vector<TH2Def> vecHistDefs2D = {
    TH2Def("hSumVsFrac", sumVsFracTitle.data(), get<0>(fracBins), get<1>(fracBins), get<2>(fracBins), get<0>(diffBins), get<1>(diffBins), get<1>(diffBins))
  };

  // run analysis -------------------------------------------------------------

  // define columns
  auto analysis = frame.Define("enePar",        getParticleEnergy, {"GeneratedParticles.type", "GeneratedParticles.energy"})
                       .Define("eneLeadHCal",   getLeadEnergy,     {"HcalBarrelClusters.energy"})
                       .Define("eneSumHCal",    getEnergySum,      {"HcalBarrelClusters.energy"})
                       .Define("eneSumECal",    getEnergySum,      {"EcalBarrelClusters.energy"})
                       .Define("diffClustHCal", getPercentDiffVec, {"enePar", "HcalBarrelClusters.energy"})
                       .Define("diffLeadHCal",  getPercentDiff,    {"enePar", "eneLeadHCal"})
                       .Define("diffSumHCal",   getPercentDiff,    {"enePar", "eneSumHCal"})
                       .Define("fracSumHxECal", getSumFraction,    {"eneSumECal", "eneSumHCal"})
                       .Define("multHit",       getMultiplicity,   {"HcalBarrelRecHits.energy"})
                       .Define("multClust",     getMultiplicity,   {"HcalBarrelClusters.energy"});

  // get 1D histograms
  //   TODO it might be nice to collect these into a vector
  //   and automate histogram operations...
  auto hEnePar         = analysis.Histo1D(vecHistDefs1D[0],  "enePar");
  auto hEneHit         = analysis.Histo1D(vecHistDefs1D[1],  "HcalBarrelRecHits.energy");
  auto hEneClust       = analysis.Histo1D(vecHistDefs1D[2],  "HcalBarrelClusters.energy");
  auto hEneLead        = analysis.Histo1D(vecHistDefs1D[3],  "eneLeadHCal");
  auto hEneSum         = analysis.Histo1D(vecHistDefs1D[4],  "eneSumHCal");
  auto hDiffClust      = analysis.Histo1D(vecHistDefs1D[5],  "diffClustHCal");
  auto hDiffLead       = analysis.Histo1D(vecHistDefs1D[6],  "diffLeadHCal");
  auto hDiffSum        = analysis.Histo1D(vecHistDefs1D[7],  "diffSumHCal");
  auto hMultHit        = analysis.Histo1D(vecHistDefs1D[8],  "multHit");
  auto hMultClust      = analysis.Histo1D(vecHistDefs1D[9],  "multClust");
  auto hMultHitInClust = analysis.Histo1D(vecHistDefs1D[10], "HcalBarrelClusters.nhits");

  // get 2D histograms
  auto hSumVsFrac = analysis.Histo2D(vecHistDefs2D[0], "fracSumHxECal", "eneSumECal");

  // make plots ---------------------------------------------------------------

  // define styles
  const tuple<uint32_t, uint32_t, uint32_t> stylePar     = {1,   1, 2};
  const tuple<uint32_t, uint32_t, uint32_t> styleHit     = {633, 1, 2};
  const tuple<uint32_t, uint32_t, uint32_t> styleClust   = {417, 1, 2};
  const tuple<uint32_t, uint32_t, uint32_t> styleLead    = {601, 1, 2};
  const tuple<uint32_t, uint32_t, uint32_t> styleSum     = {617, 1, 2};
  const tuple<uint32_t, uint32_t, uint32_t> styleFrac    = {1,   1, 2};
  const tuple<uint32_t, uint32_t, uint32_t> styleInClust = {601, 1, 2};

  // set energy styles
  hEnePar   -> SetLineColor(get<0>(stylePar));
  hEnePar   -> SetLineStyle(get<1>(stylePar));
  hEnePar   -> SetLineWidth(get<2>(stylePar));
  hEneHit   -> SetLineColor(get<0>(styleHit));
  hEneHit   -> SetLineStyle(get<1>(styleHit));
  hEneHit   -> SetLineWidth(get<2>(styleHit));
  hEneClust -> SetLineColor(get<0>(styleClust));
  hEneClust -> SetLineStyle(get<1>(styleClust));
  hEneClust -> SetLineWidth(get<2>(styleClust));
  hEneLead  -> SetLineColor(get<0>(styleLead));
  hEneLead  -> SetLineStyle(get<1>(styleLead));
  hEneLead  -> SetLineWidth(get<2>(styleLead));
  hEneSum   -> SetLineColor(get<0>(styleSum));
  hEneSum   -> SetLineStyle(get<1>(styleSum));

  // set difference styles
  hDiffClust -> SetLineColor(get<0>(styleClust));
  hDiffClust -> SetLineStyle(get<1>(styleClust));
  hDiffClust -> SetLineWidth(get<2>(styleClust));
  hDiffLead  -> SetLineColor(get<0>(styleLead));
  hDiffLead  -> SetLineStyle(get<1>(styleLead));
  hDiffLead  -> SetLineWidth(get<2>(styleLead));
  hDiffSum   -> SetLineColor(get<0>(styleSum));
  hDiffSum   -> SetLineStyle(get<1>(styleSum));
  hDiffSum   -> SetLineWidth(get<2>(styleSum));

  // set multiplicity style
  hMultHit        -> SetLineColor(get<0>(styleHit));
  hMultHit        -> SetLineStyle(get<1>(styleHit));
  hMultHit        -> SetLineWidth(get<2>(styleHit));
  hMultClust      -> SetLineColor(get<0>(styleClust));
  hMultClust      -> SetLineStyle(get<1>(styleClust));
  hMultClust      -> SetLineWidth(get<2>(styleClust));
  hMultHitInClust -> SetLineColor(get<0>(styleInClust));
  hMultHitInClust -> SetLineStyle(get<1>(styleInClust));
  hMultHitInClust -> SetLineWidth(get<2>(styleInClust));

  // canvas parameters
  const uint32_t width  = 750;
  const uint32_t height = 750;
  const uint8_t  logY   = 1;
  const uint8_t  logZ   = 1;

  TCanvas *cEne = new TCanvas("cEne", "", width, height);
  cEne      -> cd();
  cEne      -> SetLogy(logY);
  hEnePar   -> Draw();
  hEneHit   -> Draw("same");
  hEneClust -> Draw("same");
  hEneLead  -> Draw("same");
  hEneSum   -> Draw("same");

  TCanvas *cDiff = new TCanvas("cDiff", "", width, height);
  cDiff      -> cd();
  cDiff      -> SetLogy(logY);
  hDiffClust -> Draw();
  hDiffLead  -> Draw("same");
  hDiffSum   -> Draw("same");

  TCanvas *cMult = new TCanvas("cMult", "", width, height);
  cMult           -> cd();
  cMult           -> SetLogy(logY);
  hMultHit        -> Draw();
  hMultClust      -> Draw("same");
  hMultHitInClust -> Draw("same");

  TCanvas *cSumVsFrac = new TCanvas("cSumVsFrac", "", width, height);
  cSumVsFrac -> cd();
  cSumVsFrac -> SetLogz(logZ);
  hSumVsFrac -> Draw("colz");

  // save & exit --------------------------------------------------------------

  save_canvas(cEne);
  save_canvas(cDiff);
  save_canvas(cMult);
  save_canvas(cSumVsFrac);
  cEne       -> Close();
  cDiff      -> Close();
  cMult      -> Close();
  cSumVsFrac -> Close();

  // succesfully exit macro
  return 0;

}  // end 'hcal_barrel_clusters_analysis(std::string)'

// end ------------------------------------------------------------------------
