// Code to draw average number of hits vs eta at the generated level
// Shyam Kumar; shyam055119@gmail.com; shyam.kumar@ba.infn.it
#include <podio/Frame.h>
#include <podio/Reader.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHitCollection.h>

R__LOAD_LIBRARY(podioIO)
R__LOAD_LIBRARY(edm4hep)

void NhitsvsEta_ePIC(TString filePath="", TString label="", TString output_prefix=".")
  {
  
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleOffset(.85,"X");gStyle->SetTitleOffset(1.0,"Y");
   gStyle->SetTitleSize(.05,"X");gStyle->SetTitleSize(.05,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetTitleAlign(23);     
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   
   auto reader = podio::makeReader(filePath.Data());
   const unsigned nEvents = reader.getEntries("events");

   const std::vector<std::pair<std::string, std::string>> trackerHits = {
     {"VertexBarrelHits", "SVT IB"},
     {"SiBarrelHits", "SVT OB"},
     {"MPGDBarrelHits", "Inner MPGD"},
     {"OuterMPGDBarrelHits", "Outer MPGD"},
     {"TOFBarrelHits", "BTOF"},
     {"TrackerEndcapHits", "SVT Disks"},
     {"TOFEndcapHits", "ETOF"},
     {"ForwardMPGDEndcapHits", "Forward MPGD"},
     {"BackwardMPGDEndcapHits", "Backward MPGD"},
   };

   TCanvas * c1 = new TCanvas("c1","coutput",1400,1000);
   c1->SetMargin(0.10, 0.05 ,0.1,0.08);
   c1->SetGridx();
   c1->SetGridy();
   
    TProfile* hits = new TProfile("hits","Nhits (#theta)",70,-3.5,3.5);
    hits->GetXaxis()->CenterTitle();
    hits->GetYaxis()->CenterTitle();
    hits->SetMinimum(0.); 
  
    Double_t pmc = 0.;
    
     for (unsigned iEvent = 0; iEvent < nEvents; ++iEvent)
     {
      const auto frame = reader.readNextEvent();
      const bool debug = iEvent < 11;
      if (debug) printf("<--------------------Event No. %d---------------------> \n",iEvent);      
      // Generated primary track
      const auto& mcParticles = frame.get<edm4hep::MCParticleCollection>("MCParticles");
      if (mcParticles.size()>1) continue; // skip event for larger than 1 tracks
      Double_t eta_Track = 100.; // set it ouside from -3.5 to 3.5      
      for (const auto& mc : mcParticles){
      if (mc.getGeneratorStatus() !=1) continue;
      const auto& p = mc.getMomentum();
      pmc = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);   
      eta_Track = -1.0*TMath::Log(TMath::Tan((TMath::ACos(p.z/pmc))/2));
      }
      if (fabs(eta_Track)>3.5) continue; //outside tracker acceptance
      if (debug) printf("Eta of the generated track: %f, Momentum (GeV/c): %f \n",eta_Track,pmc);

      // Associated hits with the primary track (quality==0)
      int nhits = 0;
      for (const auto& [collName, detName] : trackerHits) {
        int nhits_det = 0;
        for (const auto& hit : frame.get<edm4hep::SimTrackerHitCollection>(collName)) {
          if (hit.getQuality()==0) nhits_det++;
        }
        if (debug) printf("%s Associated hits: %d \n", detName.c_str(), nhits_det);
        nhits += nhits_det;
      }
    if (nhits>0) hits->Fill(eta_Track,nhits); 

 }
 
  c1->cd();
  gPad->SetTicks(1,1);
  hits->SetTitle(";#eta_{mc};Nhits");
  hits->SetLineWidth(2);
  hits->Draw("hist");
  TPaveText *pt = new TPaveText(0.1, 0.95, 0.9, 1.0, "NDC");
  pt->AddText(Form("MC Truth-Level Hits; p_{mc} = %1.1f GeV/c",pmc));
  pt->SetBorderSize(0);     
  pt->SetFillStyle(0);      
  pt->SetTextAlign(23); 
  pt->Draw();

  c1->SaveAs(Form("%s/Nhits_vs_eta.png", output_prefix.Data()));
  c1->SaveAs(Form("%s/Nhits_vs_eta.root", output_prefix.Data()));
}
