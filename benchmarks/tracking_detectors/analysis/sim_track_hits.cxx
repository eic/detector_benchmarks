#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TProfile.h"
#include "THStack.h"
#include "Math/Vector4D.h"

#include <cstdlib>
#include <iostream>

R__LOAD_LIBRARY(libeicd.so)
R__LOAD_LIBRARY(libDD4pod.so)

#include <fmt/format.h>

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/TrackerHitCollection.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

auto p_track = [](std::vector<eicd::TrackParametersData> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(std::abs(1.0 / (in[i].qOverP)));
  }
  return result;
};

std::vector<float> pt(std::vector<edm4hep::MCParticleData> const& in) {
  std::vector<float> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(std::sqrt(in[i].momentum.x * in[i].momentum.x + in[i].momentum.y * in[i].momentum.y));
  }
  return result;
}

auto momentum = [](std::vector<ROOT::Math::PxPyPzMVector> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(in[i].E());
  }
  return result;
};
auto theta = [](std::vector<ROOT::Math::PxPyPzMVector> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(in[i].Theta() * 180 / M_PI);
  }
  return result;
};
auto fourvec = [](ROOT::VecOps::RVec<edm4hep::MCParticleData> const& in) {
  std::vector<ROOT::Math::PxPyPzMVector> result;
  ROOT::Math::PxPyPzMVector lv;
  for (size_t i = 0; i < in.size(); ++i) {
    lv.SetCoordinates(in[i].momentum.x, in[i].momentum.y, in[i].momentum.z, in[i].mass);
    result.push_back(lv);
  }
  return result;
};

auto delta_p = [](const std::vector<double>& tracks, const std::vector<double>& thrown) {
  std::vector<double> res;
  for (const auto& p1 : thrown) {
    for (const auto& p2 : tracks) {
      res.push_back(p1 - p2);
    }
  }
  return res;
};

auto delta_p_over_p = [](const std::vector<double>& tracks, const std::vector<double>& thrown) {
  std::vector<double> res;
  for (const auto& p1 : thrown) {
    for (const auto& p2 : tracks) {
      res.push_back((p1 - p2) / p1);
    }
  }
  return res;
};

// Create dataframe with extra definitions for hit multiplicities
ROOT::RDF::RNode add_subsystems(ROOT::RDF::RNode df, std::vector<std::pair<std::string, std::string>> hitcols) {
  if (hitcols.empty()) {
    return df;
  }
  const auto [name, collection] = hitcols.back();
  std::cout << " --> registering subsystem " << collection << "\n";
  auto df2 = df.Define("N_" + name, [](std::vector<edm4hep::SimTrackerHitData> hits) { return hits.size(); }, {collection});
  hitcols.pop_back();
  return add_subsystems(df2, hitcols);
};

int sim_track_hits(const char* fname = "sim_track_hits.edm4hep.root") {

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df("events", fname);

  // detect detector setup
  std::string detector = "default";
  if (const char* env_detector = std::getenv("DETECTOR_VERSION")) {
    if (detector.find("acadia") != std::string::npos) {
      detector = "acadia";
    } else if (detector.find("canyonlands") != std::string::npos) {
      detector = "canyonlands";
    }
  }
  std::cout << "sim_track_hits: detector set to " << detector << std::endl;

  // minimal hit collection setup
  std::vector<std::pair<std::string, std::string>> hitCollections{{"vtx_barrel", "VertexBarrelHits"},
                                                                  {"trk_barrel_outer", "OuterSiBarrelHits"},
                                                                  {"trk_barrel_sagitta", "SagittaSiBarrelHits"},
                                                                  {"trk_endcap_p_inner", "InnerTrackerEndcapPHits"},
                                                                  {"trk_endcap_p_middle", "MiddleTrackerEndcapPHits"},
                                                                  {"trk_endcap_p_outer", "OuterTrackerEndcapPHits"},
                                                                  {"trk_endcap_n_inner", "InnerTrackerEndcapNHits"},
                                                                  {"trk_endcap_n_middle", "MiddleTrackerEndcapNHits"},
                                                                  {"trk_endcap_n_outer", "OuterTrackerEndcapNHits"},
                                                                  {"mpgd_barrel_outer", "OuterMPGDBarrelHits"},
                                                                  {"mpgd_barrel_inner", "InnerMPGDBarrelHits"}
  };

  auto df0 = df.Define("isThrown", "MCParticles.generatorStatus == 1")
                 .Define("thrownParticles", "MCParticles[isThrown]")
                 .Define("thrownP", fourvec, {"thrownParticles"})
                 .Define("p_thrown", momentum, {"thrownP"})
                 .Define("theta_thrown", theta, {"thrownP"})
                 .Define("theta0", "theta_thrown[0]");
  auto df1 = add_subsystems(df0, hitCollections);

  // get 1D and 2D histogram definitions
  std::vector<decltype(df1.Histo2D({}, "", ""))> h2D;
  std::vector<decltype(df1.Histo1D(""))> hnth;
  std::vector<decltype(df1.Histo1D(""))> hn;

  for (const auto& [name, col] : hitCollections) {
    hnth.push_back(
        df1.Histo1D({(name + "_nhits_vs_theta").c_str(), "; #theta [deg.]; ", 20, 0, 180}, "theta0", "N_" + name));
    hn.push_back(df1.Histo1D({(name + "_nhits").c_str(), "; Nhits; #", 20, 0, 20}, "N_" + name));
    h2D.push_back(df1.Histo2D({(name + "_x_vs_y").c_str(), "; x ; y ", 100, -600, 600, 100, -600, 600},
                              col + ".position.x", col + ".position.y"));
    h2D.push_back(df1.Histo2D({(name + "_x_vs_z").c_str(), "; z ; x ", 100, -1200, 1200, 100, -900, 900},
                              col + ".position.z", col + ".position.x"));
  }
  auto hth = df1.Histo1D({"theta0", "; #theta [deg.]", 20, 0, 180}, "theta0");

  // Draw multiplicity versus theta
  {
    TCanvas c;
    THStack hs{"n_hits", "; #theta [deg.] "};
    int idx = 1;
    for (auto& h : hnth) {
      h->Divide(&*hth);
      h->SetLineColor(idx++);
      hs.Add((TH1D*)h->Clone());
    }
    hs.Draw("nostack, hist");
    c.BuildLegend();
    c.SaveAs("results/tracking_detectors/sim_track_hits_n_hits_vs_theta.png");
    c.SaveAs("results/tracking_detectors/sim_track_hits_n_hits_vs_theta.pdf");
  }
  // Draw nhits
  {
    TCanvas c;
    THStack hs{"n_hits", "; NHits "};
    int idx = 1;
    for (auto& h : hn) {
      h->SetLineColor(idx++);
      hs.Add((TH1D*)h->Clone());
    }
    hs.Draw("nostack, hist");
    c.BuildLegend();
    c.SaveAs("results/tracking_detectors/sim_track_hits_n_hits.png");
    c.SaveAs("results/tracking_detectors/sim_track_hits_n_hits.pdf");
  }

  // Draw total multiplicity versus theta
  {
    TCanvas c;
    hth->DrawClone();
    c.BuildLegend();
    c.SaveAs("results/tracking_detectors/sim_track_hits_theta.png");
    c.SaveAs("results/tracking_detectors/sim_track_hits_theta.pdf");
  }

  for (auto& h : h2D) {
    TCanvas c;
    h->DrawClone("colz");
    c.SaveAs(fmt::format("results/tracking_detectors/sim_track_hits_{}.png", h->GetName()).c_str());
    c.SaveAs(fmt::format("results/tracking_detectors/sim_track_hits_{}.pdf", h->GetName()).c_str());
  }

  return 0;
}
