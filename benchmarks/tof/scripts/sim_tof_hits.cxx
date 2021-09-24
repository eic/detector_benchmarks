#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TProfile.h"

#include <iostream>

R__LOAD_LIBRARY(libeicd.so)
R__LOAD_LIBRARY(libDD4pod.so)

#include "dd4pod/TrackerHitCollection.h"
#include "dd4pod/TrackerHitData.h"
#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/TrackParametersCollection.h"
#include "eicd/ClusterCollection.h"
#include "eicd/ClusterData.h"
#include "eicd/TrackerHitCollection.h"

using ROOT::RDataFrame;
using namespace ROOT::VecOps;

auto p_track = [](std::vector<eic::TrackParametersData> const& in) {
  std::vector<double> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(std::abs(1.0/(in[i].qOverP)));
  }
  return result;
};


std::vector<float> pt (std::vector<dd4pod::Geant4ParticleData> const& in){
  std::vector<float> result;
  for (size_t i = 0; i < in.size(); ++i) {
    result.push_back(std::sqrt(in[i].ps.x * in[i].ps.x + in[i].ps.y * in[i].ps.y));
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
   result.push_back(in[i].Theta()*180/M_PI);
  }
  return result;
};
auto fourvec = [](ROOT::VecOps::RVec<dd4pod::Geant4ParticleData> const& in) {
  std::vector<ROOT::Math::PxPyPzMVector> result;
  ROOT::Math::PxPyPzMVector lv;
  for (size_t i = 0; i < in.size(); ++i) {
    lv.SetCoordinates(in[i].ps.x, in[i].ps.y, in[i].ps.z, in[i].mass);
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
      res.push_back((p1 - p2)/p1);
    }
  }
  return res;
};

int sim_tof_hits(const char* fname = "sim_tof_hits.root")
{

  ROOT::EnableImplicitMT();
  ROOT::RDataFrame df("events", fname);

  auto df0 = df.Define("isThrown", "mcparticles.genStatus == 1")
                 .Define("thrownParticles", "mcparticles[isThrown]")
                 .Define("thrownP", fourvec, {"thrownParticles"})
                 .Define("p_thrown", momentum, {"thrownP"})
                 .Define("theta_thrown", theta, {"thrownP"})
                 .Define("theta0", "theta_thrown[0]")
                 .Define("N_VertexBarrelHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"VertexBarrelHits"})
                 .Define("N_VertexEndcapHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"VertexEndcapHits"})
                 .Define("N_BarrelHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"TrackerBarrelHits"})
                 .Define("N_EndcapHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"TrackerEndcapHits"})
                 .Define("N_BarrelTOFHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"BarrelTOFHits"})
                 .Define("N_ForwardTOFHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"ForwardTOFHits"})
                 .Define("N_BackwardTOFHits", [](std::vector<dd4pod::TrackerHitData> hits) { return hits.size();}, {"BackwardTOFHits"})
                 ;

  auto hBarrel_x_vs_y = df0.Histo2D(
      {"hBarrel_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
      "TrackerBarrelHits.position.x", "TrackerBarrelHits.position.y");
  auto hEndcap_x_vs_y = df0.Histo2D(
      {"hEndcap_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
      "TrackerEndcapHits.position.x", "TrackerEndcapHits.position.y");
  auto hVertexBarrel_x_vs_y = df0.Histo2D(
      {"hVertexBarrel_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
      "VertexBarrelHits.position.x", "VertexBarrelHits.position.y");
  auto hVertexEndcap_x_vs_y = df0.Histo2D(
      {"hVertexEndcap_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
      "VertexEndcapHits.position.x", "VertexEndcapHits.position.y");
  auto hBarrelTof_x_vs_y = df0.Histo2D(
       {"hBarrelTof_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
       "BarrelTOFHits.position.x", "BarrelTOFHits.position.y");
  auto hForwardTof_x_vs_y = df0.Histo2D(
       {"hForwardTof_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
       "ForwardTOFHits.position.x", "ForwardTOFHits.position.y");
  auto hBackwardTof_x_vs_y = df0.Histo2D(
       {"hBackwardTof_x_vs_y", "; x ; y ", 100, -600, 600, 100, -600, 600},
       "BackwardTOFHits.position.x", "BackwardTOFHits.position.y");

  auto hBarrel_x_vs_z = df0.Histo2D(
      {"hBarrel_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "TrackerBarrelHits.position.z", "TrackerBarrelHits.position.x");
  auto hEndcap_x_vs_z = df0.Histo2D(
      {"hEndcap_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "TrackerEndcapHits.position.z", "TrackerEndcapHits.position.x");
  auto hVertexBarrel_x_vs_z = df0.Histo2D(
      {"hVertexBarrel_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "VertexBarrelHits.position.z", "VertexBarrelHits.position.x");
  auto hVertexEndcap_x_vs_z = df0.Histo2D(
      {"hVertexEndcap_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "VertexEndcapHits.position.z", "VertexEndcapHits.position.x");
  auto hBarrelTof_x_vs_z = df0.Histo2D(
      {"hBarrelTof_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "BarrelTOFHits.position.z", "BarrelTOFHits.position.x");
  auto hForwrardTof_x_vs_z = df0.Histo2D(
      {"hForwardTof_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "ForwardTOFHits.position.z", "ForwardTOFHits.position.x");
  auto hBackwrardTof_x_vs_z = df0.Histo2D(
      {"hBackwardTof_x_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "BackwardTOFHits.position.z", "BackwardTOFHits.position.x");

  auto hBarrel_y_vs_z = df0.Histo2D(
      {"hBarrel_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "TrackerBarrelHits.position.z", "TrackerBarrelHits.position.y");
  auto hEndcap_y_vs_z = df0.Histo2D(
      {"hEndcap_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "TrackerEndcapHits.position.z", "TrackerEndcapHits.position.y");
  auto hVertexBarrel_y_vs_z = df0.Histo2D(
      {"hVertexBarrel_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "VertexBarrelHits.position.z", "VertexBarrelHits.position.y");
  auto hVertexEndcap_y_vs_z = df0.Histo2D(
      {"hVertexEndcap_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "VertexEndcapHits.position.z", "VertexEndcapHits.position.y");
  auto hBarrelTof_y_vs_z = df0.Histo2D(
      {"hBarrelTof_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
      "BarrelTOFHits.position.z", "BarrelTOFHits.position.y");
    auto hForwardTof_y_vs_z = df0.Histo2D(
        {"hForwardTof_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
        "ForwardTOFHits.position.z", "ForwardTOFHits.position.y");
    auto hBackwardTof_y_vs_z = df0.Histo2D(
        {"hBackwardTof_y_vs_z", "; z ; x ", 200, -1200, 1200, 100, -900, 900},
        "BackwardTOFHits.position.z", "BackwardTOFHits.position.y");
    
  auto hBarrel_N_vs_theta = df0.Histo1D({"hBarrel_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_BarrelHits");
  auto hEndcap_N_vs_theta = df0.Histo1D({"hEndcap_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_EndcapHits");
  auto hVertexBarrel_N_vs_theta = df0.Histo1D({"hVertexBarrel_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_VertexBarrelHits");
  auto hVertexEndcap_N_vs_theta = df0.Histo1D({"hVertexEndcap_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_VertexEndcapHits");
    auto hBarrelTof_N_vs_theta = df0.Histo1D({"hBarrelTof_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_BarrelTOFHits");
    auto hForwardTof_N_vs_theta = df0.Histo1D({"hForwardTof_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_ForwardTOFHits");
    auto hBackwardTof_N_vs_theta = df0.Histo1D({"hBackwardTof_N_vs_theta", "; #theta [deg.]",   20, 0, 180 }, "theta0", "N_BackwardTOFHits");

  auto hBarrel_Nhits  = df0.Histo1D({"hBarrel_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_BarrelHits");
  auto hEndcap_Nhits  = df0.Histo1D({"hEndcap_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_EndcapHits");
  auto hVertexBarrel_Nhits  = df0.Histo1D({"hVertexBarrel_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_VertexBarrelHits");
  auto hVertexEndcap_Nhits  = df0.Histo1D({"hVertexEndcap_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_VertexEndcapHits");
    auto hBarrelTof_Nhits  = df0.Histo1D({"hBarrelTof_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_BarrelTOFHits");
    auto hForwardTof_Nhits  = df0.Histo1D({"hForwardTof_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_ForwardTOFHits");
    auto hBackwardTof_Nhits  = df0.Histo1D({"hBackwardTof_Nhits", "; #theta [deg.]",   20, 0, 20 }, "N_BackwardTOFHits");

  auto hBarrel_Ntheta = df0.Histo1D({"hBarrel_Ntheta", "; #theta [deg.]",   20, 0, 180 }, "theta0");
  auto hEndcap_Ntheta = df0.Histo1D({"hEndcap_Ntheta", "; #theta [deg.]",   20, 0, 180 }, "theta0");
  auto hVertexBarrel_Ntheta = df0.Histo1D({"hVertexBarrel_Ntheta", "; #theta [deg.]",   20, 0, 180 }, "theta0");
  auto hVertexEndcap_Ntheta = df0.Histo1D({"hVertexEndcap_Ntheta", "; #theta [deg.]",   20, 0, 180 }, "theta0");
    auto hBarrelTof_Ntheta = df0.Histo1D({"hBarrelTof_Ntheta", "; #theta [deg.]",   20, 0, 180 }, "theta0");

  auto c = new TCanvas();
  auto hs = new THStack("n_hits","; #theta  ");
  auto h1 = (TH1D*) hBarrel_N_vs_theta->Clone();
  auto h2 = (TH1D*) hBarrel_Ntheta->Clone();
  h1->Divide(h2);
  hs->Add(h1);
  h1 = (TH1D*) hEndcap_N_vs_theta->Clone();
  h2 = (TH1D*) hEndcap_Ntheta->Clone();
  h1->Divide(h2);
  h1->SetLineColor(2);
  hs->Add(h1);
  h1 = (TH1D*) hVertexEndcap_N_vs_theta->Clone();
  h2 = (TH1D*) hVertexEndcap_Ntheta->Clone();
  h1->Divide(h2);
  h1->SetLineColor(4);
  hs->Add(h1);
  h1 = (TH1D*) hVertexBarrel_N_vs_theta->Clone();
  h2 = (TH1D*) hVertexBarrel_Ntheta->Clone();
  h1->Divide(h2);
  h1->SetLineColor(8);
  hs->Add(h1);
  hs->Draw("nostack, hist");
  c->BuildLegend();
  c->SaveAs("results/tof/sim_tof_hits_n_hits_vs_theta.png");
  c->SaveAs("results/tof/sim_tof_hits_n_hits_vs_theta.pdf");

  c  = new TCanvas();
  hs = new THStack("theta","; #theta  ");
  h1 = (TH1D*) hBarrel_N_vs_theta->Clone();
  h2 = (TH1D*) hBarrel_Ntheta->Clone();
  //h1->Divide(h2);
  hs->Add(h2);
  h1 = (TH1D*) hEndcap_N_vs_theta->Clone();
  h2 = (TH1D*) hEndcap_Ntheta->Clone();
  //h1->Divide(h2);
  h1->SetLineColor(2);
  h2->SetLineColor(2);
  hs->Add(h2);
  //h1 = (TH1D*) hVtxBarrel_vs_theta->Clone();
  //h1->SetLineColor(4);
  //h1->SetFillStyle(3001);
  //h1->SetFillColor(4);
  //hs->Add(h1);
  hs->Draw("nostack hist");
  c->BuildLegend();
  c->SaveAs("results/tof/sim_tof_hits_theta.png");
  c->SaveAs("results/tof/sim_tof_hits_theta.pdf");

  c  = new TCanvas();
  hs = new THStack("hits","; hits  ");
  h1 = (TH1D*) hBarrel_Nhits->Clone();
  hs->Add(h1);
  h1 = (TH1D*) hEndcap_Nhits->Clone();
  h1->SetLineColor(2);
  h2->SetLineColor(2);
  hs->Add(h2);
  //h1 = (TH1D*) hVtxBarrel_Nhits->Clone();
  //h1->SetLineColor(4);
  //h1->SetFillStyle(3001);
  //h1->SetFillColor(4);
  //hs->Add(h1);
  //hs->Draw("nostack hist");
  c->BuildLegend();
  c->SaveAs("results/tof/sim_tof_hits_nhits.png");
  c->SaveAs("results/tof/sim_tof_hits_nhits.pdf");

  c = new TCanvas();
  hBarrelTof_x_vs_y->DrawCopy("colz");
  c->SaveAs("results/tof/sim_tof_hits_tofBarrel_xy.png");
  c->SaveAs("results/tof/sim_tof_hits_tofBarrel_xy.pdf");
  c = new TCanvas();
  hForwardTof_x_vs_y->DrawCopy("colz");
  c->SaveAs("results/tof/sim_tof_hits_tofForward_xy.png");
  c->SaveAs("results/tof/sim_tof_hits_tofForward_xy.pdf");
  c = new TCanvas();
  hBackwardTof_x_vs_y->DrawCopy("colz");
  c->SaveAs("results/tof/sim_tof_hits_tofBackward_xy.png");
  c->SaveAs("results/tof/sim_tof_hits_tofBackward_xy.pdf");

  c = new TCanvas();
  hBarrel_x_vs_y->DrawCopy("colz");
  c->SaveAs("results/tof/sim_tof_hits_trkBarrel_xy.png");
  c->SaveAs("results/tof/sim_tof_hits_trkBarrel_xy.pdf");

  c = new TCanvas();
  hEndcap_x_vs_y->DrawCopy("colz");
  hVertexEndcap_x_vs_y->DrawCopy("colz same");
  //hBarrel_x_vs_y->DrawCopy("colz same");
  //hvtxEndcap_x_vs_y->DrawCopy("colz same");

  c->SaveAs("results/tof/sim_tof_hits_Hits_xy.png");
  c->SaveAs("results/tof/sim_tof_hits_Hits_xy.pdf");

  //hAllHits_x_vs_z->DrawCopy("colz");
  hBarrel_x_vs_z->DrawCopy("colz");
  hEndcap_x_vs_z->DrawCopy("colz same");
  c->SaveAs("results/tof/sim_tof_hits_Hits_xz.png");
  c->SaveAs("results/tof/sim_tof_hits_Hits_xz.pdf");

  hBarrel_y_vs_z->DrawCopy("colz");
  hEndcap_y_vs_z->DrawCopy("colz same");
  c->SaveAs("results/tof/sim_tof_hits_Hits_yz.png");
  c->SaveAs("results/tof/sim_tof_hits_Hits_yz.pdf");
  return 0;
}
