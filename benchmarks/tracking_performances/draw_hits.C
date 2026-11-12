// Hits distribution of detectors
// Shyam Kumar:INFN Bari, shyam.kumar@ba.infn.it; shyam055119@gmail.com

#include <TGraphErrors.h>
#include <TF1.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>

#include <podio/Frame.h>
#include <podio/Reader.h>
#include <edm4hep/SimTrackerHitCollection.h>

R__LOAD_LIBRARY(podioIO)
R__LOAD_LIBRARY(edm4hep)

void fill_hits(const TString& filename, const std::vector<std::pair<std::string, std::function<void(const edm4hep::SimTrackerHit&)>>>& fillers)
{
  auto reader = podio::makeReader(filename.Data());
  const unsigned nEvents = reader.getEntries("events");
  for (unsigned i = 0; i < nEvents; ++i) {
    const auto frame = reader.readNextEvent();
    for (const auto& [collName, fill] : fillers) {
      for (const auto& hit : frame.get<edm4hep::SimTrackerHitCollection>(collName)) {
        fill(hit);
      }
    }
  }
}

void draw_hits(TString filename="", TString output_prefix=".")
{

//==========Style of the plot============
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(0);
   gStyle->SetTitleOffset(.85,"X");gStyle->SetTitleOffset(.85,"Y");
   gStyle->SetTitleSize(.04,"X");gStyle->SetTitleSize(.04,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);

//=======Reading the root file DD4HEP===========

 // Timer Start
  TStopwatch timer;
  timer.Start();

  TCanvas *c1 = new TCanvas("c1","c1",1200,1000);
  c1->SetMargin(0.09, 0.03 ,0.1,0.06);

 // X-Y Hits
  Int_t nbins = 320;
  Double_t x= 100., y = 100.;
  TH2D *hitsxy_vtx_si = new TH2D("hitsxy_vtx_si","hitsxy_vtx_si",nbins,-x,x,nbins,-y,y);
  TH2D *hitsxy_barrel_si = new TH2D("hitsxy_barrel_si","hitsxy_barrel_si",nbins,-x,x,nbins,-y,y);
  TH2D *hitsxy_barrel_mm_in = new TH2D("hitsxy_barrel_mm_in","hitsxy_barrel_mm_in",nbins,-x,x,nbins,-y,y);
  TH2D *hitsxy_barrel_tof = new TH2D("hitsxy_barrel_tof","hitsxy_barrel_tof",nbins,-x,x,nbins,-y,y);
  TH2D *hitsxy_barrel_mm_out = new TH2D("hitsxy_barrel_mm_out","hitsxy_barrel_mm_out",nbins,-x,x,nbins,-y,y);

  // Multiply by 0.1 for cm
  auto xy = [](TH2D* h) { return [h](const edm4hep::SimTrackerHit& hit) { h->Fill(hit.getPosition().x*0.1, hit.getPosition().y*0.1); }; };
  fill_hits(filename, {
    {"VertexBarrelHits", xy(hitsxy_vtx_si)},
    {"SiBarrelHits", xy(hitsxy_barrel_si)},
    {"MPGDBarrelHits", xy(hitsxy_barrel_mm_in)},
    {"TOFBarrelHits", xy(hitsxy_barrel_tof)},
    {"OuterMPGDBarrelHits", xy(hitsxy_barrel_mm_out)},
  });

        hitsxy_vtx_si->SetMarkerStyle(31);
        hitsxy_vtx_si->SetTitle("Hit Points (XY)");
        hitsxy_vtx_si->SetMarkerColor(kBlack);
        hitsxy_vtx_si->SetMarkerSize(0.1);
        hitsxy_vtx_si->SetLineColor(kBlack);
        hitsxy_vtx_si->GetXaxis()->SetTitle("X [cm]");
        hitsxy_vtx_si->GetYaxis()->SetTitle("Y [cm]");
        hitsxy_vtx_si->GetXaxis()->CenterTitle();
        hitsxy_vtx_si->GetYaxis()->CenterTitle();

        hitsxy_barrel_si->SetMarkerStyle(20);
        hitsxy_barrel_si->SetMarkerSize(0.1);
        hitsxy_barrel_si->SetMarkerColor(kMagenta);
        hitsxy_barrel_si->SetLineColor(kMagenta);

        hitsxy_barrel_mm_in->SetMarkerStyle(20);
        hitsxy_barrel_mm_in->SetMarkerSize(0.1);
        hitsxy_barrel_mm_in->SetMarkerColor(kBlue);
        hitsxy_barrel_mm_in->SetLineColor(kBlue);

        hitsxy_barrel_tof->SetMarkerStyle(20);
        hitsxy_barrel_tof->SetMarkerSize(0.1);
        hitsxy_barrel_tof->SetMarkerColor(kGreen);
        hitsxy_barrel_tof->SetLineColor(kGreen);hitsxy_barrel_mm_out->SetMarkerStyle(20);
        hitsxy_barrel_mm_out->SetMarkerSize(0.1);
        hitsxy_barrel_mm_out->SetMarkerColor(kBlue-7);
        hitsxy_barrel_mm_out->SetLineColor(kBlue-7);
        
        c1->cd();
        hitsxy_vtx_si->Draw("BOX");
        hitsxy_barrel_si->Draw("BOX same");
        hitsxy_barrel_mm_in->Draw("BOX same");
        hitsxy_barrel_tof->Draw("BOX same");
        hitsxy_barrel_mm_out->Draw("BOX same");

 TLegend *l= new TLegend(0.65,0.85,0.90,1.0);
 l->SetTextSize(0.025);
 l->SetBorderSize(0);
 l->AddEntry(hitsxy_vtx_si,"VertexBarrelHits");
 l->AddEntry(hitsxy_barrel_si,"SiBarrelHits");
 l->AddEntry(hitsxy_barrel_mm_in,"MPGDBarrelHits");
 l->AddEntry(hitsxy_barrel_tof,"TOFBarrelHits");
 l->AddEntry(hitsxy_barrel_mm_out,"OuterMPGDBarrelHits");
 l->Draw();
 c1->SaveAs(Form("%s/hitsxy_dd4hep.png", output_prefix.Data()));
 c1->SaveAs(Form("%s/hitsxy_dd4hep.eps", output_prefix.Data()));
 c1->SaveAs(Form("%s/hitsxy_dd4hep.root", output_prefix.Data()));

  TCanvas *c2 = new TCanvas("c2","c2",1200,1000);
  c2->SetMargin(0.09, 0.03 ,0.1,0.06);
 // Y-Z Hits
 Int_t nbinsx = 400, nbinsy=1800.;
 x= 200.; y = 90.;
 double xmin = -200.;

   TH2D *hitsrz_vtx_si = new TH2D("hitsrz_vtx_si","hitsrz_vtx_si",nbinsx,xmin,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_si = new TH2D("hitsrz_barrel_si","hitsrz_barrel_si",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_mm_in = new TH2D("hitsrz_barrel_mm_in","hitsrz_barrel_mm_in",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_mm_out = new TH2D("hitsrz_barrel_mm_out","hitsrz_barrel_mm_out",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_tof = new TH2D("hitsrz_barrel_tof","hitsrz_barrel_tof",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_disks_si = new TH2D("hitsrz_disks_si","hitsrz_disks_si",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_endcap_tof = new TH2D("hitsrz_endcap_tof","hitsrz_endcap_tof",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_fwd_mpgd = new TH2D("hitsrz_fwd_mpgd","hitsrz_fwd_mpgd",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_bwd_mpgd = new TH2D("hitsrz_bwd_mpgd","hitsrz_bwd_mpgd",nbinsx,-x,x,nbinsy,-y,y);

   TH2D *hitsrz_vtx_si_1 = new TH2D("hitsrz_vtx_si_1","hitsrz_vtx_si_1",nbinsx,xmin,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_si_1 = new TH2D("hitsrz_barrel_si_1","hitsrz_barrel_si_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_mm_in_1 = new TH2D("hitsrz_barrel_mm_in_1","hitsrz_barrel_mm_in_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_mm_out_1 = new TH2D("hitsrz_barrel_mm_out_1","hitsrz_barrel_mm_out_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_barrel_tof_1 = new TH2D("hitsrz_barrel_tof_1","hitsrz_barrel_tof_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_disks_si_1 = new TH2D("hitsrz_disks_si_1","hitsrz_disks_si_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_endcap_tof_1 = new TH2D("hitsrz_endcap_tof_1","hitsrz_endcap_tof_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_fwd_mpgd_1 = new TH2D("hitsrz_fwd_mpgd_1","hitsrz_fwd_mpgd_1",nbinsx,-x,x,nbinsy,-y,y);
   TH2D *hitsrz_bwd_mpgd_1 = new TH2D("hitsrz_bwd_mpgd_1","hitsrz_bwd_mpgd_1",nbinsx,-x,x,nbinsy,-y,y);

    // R signed by y, multiply by 0.1 for cm
    auto rz = [](TH2D* hpos, TH2D* hneg) {
      return [hpos, hneg](const edm4hep::SimTrackerHit& hit) {
        const auto& p = hit.getPosition();
        const double r = std::sqrt(p.x*p.x + p.y*p.y)*0.1;
        if (p.y > 0) hpos->Fill(p.z*0.1, r);
        else if (p.y < 0) hneg->Fill(p.z*0.1, -r);
      };
    };
    fill_hits(filename, {
      {"VertexBarrelHits", rz(hitsrz_vtx_si, hitsrz_vtx_si_1)},
      {"SiBarrelHits", rz(hitsrz_barrel_si, hitsrz_barrel_si_1)},
      {"MPGDBarrelHits", rz(hitsrz_barrel_mm_in, hitsrz_barrel_mm_in_1)},
      {"TOFBarrelHits", rz(hitsrz_barrel_tof, hitsrz_barrel_tof_1)},
      {"TrackerEndcapHits", rz(hitsrz_disks_si, hitsrz_disks_si_1)},
      {"TOFEndcapHits", rz(hitsrz_endcap_tof, hitsrz_endcap_tof_1)},
      {"OuterMPGDBarrelHits", rz(hitsrz_barrel_mm_out, hitsrz_barrel_mm_out_1)},
      {"ForwardMPGDEndcapHits", rz(hitsrz_fwd_mpgd, hitsrz_fwd_mpgd_1)},
      {"BackwardMPGDEndcapHits", rz(hitsrz_bwd_mpgd, hitsrz_bwd_mpgd_1)},
    });



          hitsrz_vtx_si->SetMarkerStyle(31);
          hitsrz_vtx_si->SetTitle("");
          hitsrz_vtx_si->SetMarkerSize(0.1);
          hitsrz_vtx_si->SetLineColor(kBlack);
          hitsrz_vtx_si->SetMarkerColor(kBlack);
          hitsrz_vtx_si->GetXaxis()->SetTitle("Z [cm]");
          hitsrz_vtx_si->GetYaxis()->SetTitle("R [cm]");
          hitsrz_vtx_si->GetXaxis()->CenterTitle();
          hitsrz_vtx_si->GetYaxis()->CenterTitle();

          hitsrz_vtx_si_1->SetMarkerSize(0.1);
          hitsrz_vtx_si_1->SetLineColor(kBlack);
          hitsrz_vtx_si_1->SetMarkerColor(kBlack);

          hitsrz_barrel_si->SetMarkerStyle(20);
          hitsrz_barrel_si->SetMarkerSize(0.1);
          hitsrz_barrel_si->SetMarkerColor(kMagenta);
          hitsrz_barrel_si->SetLineColor(kMagenta);

          hitsrz_barrel_si_1->SetMarkerSize(0.1);
          hitsrz_barrel_si_1->SetLineColor(kMagenta);
          hitsrz_barrel_si_1->SetMarkerColor(kMagenta);

          hitsrz_barrel_mm_in->SetMarkerStyle(20);
          hitsrz_barrel_mm_in->SetMarkerSize(0.1);
          hitsrz_barrel_mm_in->SetLineColor(kBlue);
          hitsrz_barrel_mm_in->SetMarkerColor(kBlue);
          hitsrz_barrel_mm_in_1->SetMarkerSize(0.1);
          hitsrz_barrel_mm_in_1->SetLineColor(kBlue);
          hitsrz_barrel_mm_in_1->SetMarkerColor(kBlue);

          hitsrz_barrel_tof->SetMarkerStyle(20);
          hitsrz_barrel_tof->SetMarkerSize(0.1);
          hitsrz_barrel_tof->SetLineColor(kGreen);
          hitsrz_barrel_tof->SetMarkerColor(kGreen);
          hitsrz_barrel_tof_1->SetMarkerSize(0.1);
          hitsrz_barrel_tof_1->SetLineColor(kGreen);
          hitsrz_barrel_tof_1->SetMarkerColor(kGreen);


          hitsrz_barrel_mm_out->SetMarkerStyle(20);
          hitsrz_barrel_mm_out->SetMarkerSize(0.1);
          hitsrz_barrel_mm_out->SetMarkerColor(kBlue-7);
          hitsrz_barrel_mm_out->SetLineColor(kBlue-7);

          hitsrz_barrel_mm_out_1->SetMarkerSize(0.1);
          hitsrz_barrel_mm_out_1->SetLineColor(kBlue-7);
          hitsrz_barrel_mm_out_1->SetMarkerColor(kBlue-7);
          hitsrz_endcap_tof->SetMarkerStyle(20);
          hitsrz_endcap_tof->SetMarkerSize(0.1);
          hitsrz_endcap_tof->SetMarkerColor(kCyan);
          hitsrz_endcap_tof->SetLineColor(kCyan);

          hitsrz_endcap_tof_1->SetMarkerSize(0.1);
          hitsrz_endcap_tof_1->SetLineColor(kCyan);
          hitsrz_endcap_tof_1->SetMarkerColor(kCyan);

          hitsrz_disks_si->SetMarkerStyle(20);
          hitsrz_disks_si->SetMarkerSize(0.1);
          hitsrz_disks_si->SetMarkerColor(kRed);
          hitsrz_disks_si->SetLineColor(kRed);

          hitsrz_disks_si_1->SetMarkerSize(0.1);
          hitsrz_disks_si_1->SetLineColor(kRed);
          hitsrz_disks_si_1->SetMarkerColor(kRed);

          hitsrz_fwd_mpgd->SetMarkerSize(0.1);
          hitsrz_fwd_mpgd->SetLineColor(kRed-7);
          hitsrz_fwd_mpgd->SetMarkerColor(kRed-7);
          hitsrz_fwd_mpgd_1->SetMarkerSize(0.1);
          hitsrz_fwd_mpgd_1->SetLineColor(kRed-7);
          hitsrz_fwd_mpgd_1->SetMarkerColor(kRed-7);
          hitsrz_bwd_mpgd->SetMarkerSize(0.1);
          hitsrz_bwd_mpgd->SetLineColor(kOrange);
          hitsrz_bwd_mpgd->SetMarkerColor(kOrange);

          hitsrz_bwd_mpgd_1->SetMarkerSize(0.1);
          hitsrz_bwd_mpgd_1->SetLineColor(kOrange);
          hitsrz_bwd_mpgd_1->SetMarkerColor(kOrange);

        c2->cd();
        hitsrz_vtx_si->Draw("BOX");
        hitsrz_vtx_si_1->Draw("BOX same");
        hitsrz_barrel_si->Draw("BOX same");
        hitsrz_barrel_si_1->Draw("BOX same");
        hitsrz_barrel_mm_in->Draw("BOX same");
        hitsrz_barrel_mm_in_1->Draw("BOX same");
        hitsrz_barrel_tof->Draw("BOX same");
        hitsrz_barrel_tof_1->Draw("BOX same");
        hitsrz_barrel_mm_out->Draw("BOX same");
        hitsrz_barrel_mm_out_1->Draw("BOX same");
        hitsrz_endcap_tof->Draw("BOX same");
        hitsrz_endcap_tof_1->Draw("BOX same");
        hitsrz_disks_si->Draw("BOX same");
        hitsrz_disks_si_1->Draw("BOX same");
        hitsrz_fwd_mpgd->Draw("BOX same");
        hitsrz_fwd_mpgd_1->Draw("BOX same");
        hitsrz_bwd_mpgd->Draw("BOX same");
        hitsrz_bwd_mpgd_1->Draw("BOX same");

  TLegend *l1= new TLegend(0.11,0.88,0.95,0.99);
  l1->SetNColumns(3);
  l1->SetTextSize(0.025);
  l1->SetBorderSize(0);
  l1->AddEntry(hitsrz_vtx_si,"VertexBarrelHits");
  l1->AddEntry(hitsrz_barrel_si,"SiBarrelHits");
  l1->AddEntry(hitsrz_barrel_mm_in,"MPGDBarrelHits");
  l1->AddEntry(hitsrz_barrel_tof,"TOFBarrelHits");
  l1->AddEntry(hitsrz_barrel_mm_out,"OuterMPGDBarrelHits");
  l1->AddEntry(hitsrz_disks_si,"TrackerEndcapHits");
  l1->AddEntry(hitsrz_endcap_tof,"TOFEndcapHits");
  l1->AddEntry(hitsrz_fwd_mpgd,"ForwardMPGDEndcapHits");
  l1->AddEntry(hitsrz_bwd_mpgd,"BackwardMPGDEndcapHits");
  l1->Draw();

  c2->SaveAs(Form("%s/hitsrz_dd4hep.png", output_prefix.Data()));
  c2->SaveAs(Form("%s/hitsrz_dd4hep.eps", output_prefix.Data()));
  c2->SaveAs(Form("%s/hitsrz_dd4hep.root", output_prefix.Data()));
 // Timer Stop
  timer.Stop();
  Double_t realtime = timer.RealTime();
  Double_t cputime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",realtime,cputime);

}
