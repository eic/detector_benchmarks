// Macro to draw the pulls of the track parameters
// Shyam Kumar; INFN Bari, shyam.kumar@ba.infn.it

void draw_Pulls(TString particle = "pi-", double etamin = -1.0, double etamax = 1.0){
	
	gStyle->SetPalette(kRainBow);
	gStyle->SetTitleSize(0.045,"XY");	
	gStyle->SetTitleSize(0.04,"XY");	
	gStyle->SetLabelSize(0.04,"XY");	
	gStyle->SetTitleOffset(1.0,"XY");	
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
	gStyle->SetOptTitle(1);
	gStyle->SetGridColor(kBlack);     
	gStyle->SetGridWidth(2);        
	gStyle->SetGridStyle(2);

  const Int_t nfiles = 6;
  double mom[nfiles] ={0.5,1.0,2.0,5.0,10.0,20.0};
  TFile *fmom_real[nfiles];
  
  for (int i =0; i<nfiles; ++i){
      
  TCanvas *can = new TCanvas("can","can",2000,1000);
  can->Divide(3,2);

	fmom_real[i] = TFile::Open(Form("./realseed/pi-/mom/Performances_mom_%1.1f_mom_resol_realseed_%s.root",mom[i],particle.Data()));
	TH1D *hpull_invp = (TH1D*) fmom_real[i]->Get(Form("hpull_invp_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	TH1D *hpull_d0xy = (TH1D*) fmom_real[i]->Get(Form("hpull_d0xy_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	TH1D *hpull_d0z = (TH1D*) fmom_real[i]->Get(Form("hpull_d0z_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	TH1D *hpull_phi = (TH1D*) fmom_real[i]->Get(Form("hpull_phi_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	TH1D *hpull_theta = (TH1D*) fmom_real[i]->Get(Form("hpull_theta_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
	
	can->cd(1);
	hpull_invp->Draw("hist")

	can->cd(2);
	hpull_d0xy->Draw("hist")

	can->cd(3);
	hpull_d0z->Draw("hist")

    can->cd(4);
	hpull_phi->Draw("hist")
	
	can->cd(5);
	hpull_theta->Draw("hist")
	can->SaveAs(Form("Final_Results/%s/mom/hpulls_%1.1f_%1.1f_pmax_%1.1f.png",particle.Data(),etamin,etamax,mom[i]));
	can->SaveAs(Form("Final_Results/%s/mom/hpulls_%1.1f_%1.1f_pmax_%1.1f.root",particle.Data(),etamin,etamax,mom[i]));
	}

}   



















