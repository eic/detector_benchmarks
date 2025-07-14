// Code to extract the Tracking Performances
// Shyam Kumar; INFN Bari, Italy
// shyam.kumar@ba.infn.it; shyam.kumar@cern.ch

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TVector3.h"

#define mpi 0.139  // 1.864 GeV/c^2

void LFHCAL_Performance(TString filename="tracking_output",TString particle="pi-", double mom=0.1, Double_t pTcut = 0.0, TString name = "")
{

  // style of the plot
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleOffset(.85,"X");gStyle->SetTitleOffset(.85,"Y");
   gStyle->SetTitleSize(.05,"X");gStyle->SetTitleSize(.05,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(1);
   
   TString dir = "";
   TString dist_dir_mom = "mom_resol";
      
   bool debug=true;	  
  // Tree with reconstructed tracks
   const int nbins_eta = 5;
   int nfiles = 100; 
   double eta[nbins_eta+1]={1.2,1.5,2,2.5,3,3.5};
   double pt[nbins_eta+1]={0.5,1.0,2.0,5.0,10.0,20.1};
   TH1D *histp[nbins_eta]; 
   
   
   for (int i=0; i<nbins_eta; i++){
   histp[i] = new TH1D(Form("hist_etabin%d",i),Form("hist_etabin%d",i),600,-1,1);
   histp[i]->SetTitle(Form("%1.1f < #eta < %1.1f && p = %1.1f ",eta[i],eta[i+1],mom));
   histp[i]->SetName(Form("hist_mom_%1.1f_%1.1f_pmax_%1.1f",mom,eta[i],eta[i+1]));
   }
   
   TFile* file = TFile::Open(filename.Data());
   if (!file) {printf("file not found !!!"); return;}
   TTreeReader myReader("events", file); // name of tree and file
   if (debug) cout<<"Filename: "<<file->GetName()<<"\t NEvents: "<<myReader.GetEntries()<<endl;
  
   // MC and Reco information 
   TTreeReaderArray<Float_t> charge(myReader, "MCParticles.charge"); 
   TTreeReaderArray<Double_t> vx_mc(myReader, "MCParticles.vertex.x");
   TTreeReaderArray<Double_t> vy_mc(myReader, "MCParticles.vertex.y");
   TTreeReaderArray<Double_t> vz_mc(myReader, "MCParticles.vertex.z");
   TTreeReaderArray<Double_t> px_mc(myReader, "MCParticles.momentum.x");
   TTreeReaderArray<Double_t> py_mc(myReader, "MCParticles.momentum.y");
   TTreeReaderArray<Double_t> pz_mc(myReader, "MCParticles.momentum.z");
   TTreeReaderArray<Int_t> status(myReader, "MCParticles.generatorStatus"); 
   TTreeReaderArray<Int_t> pdg(myReader, "MCParticles.PDG");

   TTreeReaderArray<Float_t> pe_lc(myReader, "LFHCALClusters.energy"); 
   TTreeReaderArray<Float_t> px_lc(myReader, "LFHCALClusters.position.x"); 
   TTreeReaderArray<Float_t> py_lc(myReader, "LFHCALClusters.position.y"); 
   TTreeReaderArray<Float_t> pz_lc(myReader, "LFHCALClusters.position.z"); 
   TTreeReaderArray<Float_t> pe_ec(myReader, "EcalEndcapPClusters.energy"); 
   TTreeReaderArray<Float_t> px_ec(myReader, "EcalEndcapPClusters.position.x"); 
   TTreeReaderArray<Float_t> py_ec(myReader, "EcalEndcapPClusters.position.y"); 
   TTreeReaderArray<Float_t> pz_ec(myReader, "EcalEndcapPClusters.position.z"); 

  int count =0;
  int matchId = 1; // Always matched track assigned the index 0 
  while (myReader.Next()) 
    {
      std::cout << "events = " << count++ << std::endl;
      for (int j = 0; j < pdg.GetSize(); ++j)
	{
	  if (status[j] !=1 && pdg.GetSize()!=1) continue;
	  Double_t pzmc = pz_mc[j];  
	  Double_t ptmc = sqrt(px_mc[j]*px_mc[j]+py_mc[j]*py_mc[j]); 
	  Double_t pmc = sqrt(px_mc[j]*px_mc[j]+py_mc[j]*py_mc[j]+pz_mc[j]*pz_mc[j]); // 1./(q/p); similar to prec
	  Double_t etamc = -1.0*TMath::Log(TMath::Tan((TMath::ACos(pzmc/fabs(pmc)))/2));
	  Double_t phimc = TMath::ATan2(py_mc[j],px_mc[j]);
	  std::cout << "neutron p=" << pmc << " pt=" << ptmc << std::endl;
	  
	  if (fabs(ptmc) < pTcut) continue;

	  float l_px_tot=0;
	  float l_py_tot=0;
	  float l_pz_tot=0;
	  float l_e_tot=0;

	  std::cout << "LFHCAL nclus=" << px_lc.GetSize() << " ECAL nclus=" << pe_ec.GetSize() << std::endl;
	  for (int jl = 0;jl<px_lc.GetSize();jl++)
	    {
	      float e = pe_lc[jl];
	      TVector3 v(px_lc[jl],py_lc[jl],pz_lc[jl]);
	      v.Print();
	      float eta = v.PseudoRapidity();
	      float phi = v.Phi();
	      float pt = e/cosh(eta);
	      std::cout << "LFHCAL clus: e=" << e << " eta=" << eta << " pt=" << pt << std::endl;
	      l_e_tot += e;
	      l_px_tot += pt*cos(phi);
	      l_py_tot += pt*sin(phi);
	      l_pz_tot += pt*sinh(eta);
	    }
	  
	  float e_px_tot=0;
	  float e_py_tot=0;
	  float e_pz_tot=0;
	  float e_e_tot=0;

	  for (int je = 0;je<px_ec.GetSize();je++)
	    {
	      float e = pe_ec[je];
	      TVector3 v(px_ec[je],py_ec[je],pz_ec[je]);
	      float eta = v.PseudoRapidity();
	      float phi = v.Phi();
	      float pt = e/cosh(eta);
	      std::cout << "ECAL clus: e=" << e << " eta=" << eta << " pt=" << pt << std::endl;
	      e_e_tot += e;
	      e_px_tot += pt*cos(phi);
	      e_py_tot += pt*sin(phi);
	      e_pz_tot += pt*sinh(eta);
	    }

	  std::cout << "LFHCAL e=" <<l_e_tot << " ECAL e=" << e_e_tot << std::endl;
	  float px_tot = l_px_tot+e_px_tot;
	  float py_tot = l_py_tot+e_py_tot;
	  float pz_tot = l_pz_tot+e_pz_tot;
	  float e_tot = l_e_tot+e_e_tot;
	  
	  float prec = sqrt(px_tot*px_tot+py_tot*py_tot+pz_tot*pz_tot);
	  float ptrec = sqrt(px_tot*px_tot+py_tot*py_tot);
	  float pzrec = pz_tot;
	  
	  float p_resol = (e_tot-pmc)/pmc;
	  std::cout << "p_resol = " << p_resol << std::endl;
	  for (int ibin=0; ibin<nbins_eta; ++ibin){ 
	    if(etamc>eta[ibin] && etamc<eta[ibin+1]) histp[ibin]->Fill(p_resol); 
	  }
	} // Generated Tracks  
  
    }// event loop ends    
  
   TFile *fout_mom = new TFile(Form("%s/mom/lfhcal_mom_%1.1f_%s_%s.root",particle.Data(),mom,dist_dir_mom.Data(),particle.Data()),"recreate");
   fout_mom->cd();
   for (int ibin=0; ibin<nbins_eta; ++ibin) histp[ibin]->Write();
   fout_mom->Close();

}




