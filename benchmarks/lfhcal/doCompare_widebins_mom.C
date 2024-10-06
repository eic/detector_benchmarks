// Code to compare the tracking performances: Truth seeding vs real seeding
// Shyam Kumar; shyam.kumar@ba.infn.it; shyam055119@gmail.com

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#define mpi 0.139  // 1.864 GeV/c^2

//void draw_req_Mom(double etamin, double etamax, double xmin=0., double xmax=0.);
void doCompare_widebins_mom(TString particle = "pi-",double etamin=-1.0, double etamax=1.0, double range =0.3, Bool_t drawreq=1, TString extra_legend = "") // name = p, pt for getting p or pt dependence fitted results
{
  
  //=== style of the plot=========
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(1.0,"XY");
  gStyle->SetTitleSize(.04,"XY");
  gStyle->SetLabelSize(.04,"XY");
  gStyle->SetHistLineWidth(2);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  
  const Int_t nfiles = 6;
  double mom[nfiles] ={0.5,1.0,2.0,5.0,10.0,20.0};
  std::vector<double> momV,  momresolV, err_momresolV;
  momV.clear(); momresolV.clear(); err_momresolV.clear();
  TString symbolname = "";
  if (particle == "pi-") symbolname = "#pi^{-}"; 
  else symbolname = particle; 
  ofstream outfile;
  outfile.open ("Mom_resol.txt",ios_base::app);  
  
  TF1 *f1=new TF1("f1","FitMomentumResolution",0.,30.0,2);
  f1->SetParLimits(0,0.,0.1);	
  f1->SetParLimits(1,0.,5.0);	
  
  TCanvas *c_mom = new TCanvas("cmom","cmom",1400,1000);
  c_mom->SetMargin(0.10, 0.05 ,0.1,0.05);
  c_mom->SetGridy();
  
  //Reading the root file
  TFile *fmom[nfiles];
  TGraphErrors *gr_mom;
  TMultiGraph *mgMom; 
  TLegend *lmom; 
  mgMom = new TMultiGraph("mgMom",";p (GeV/c); #sigmap/p %");
  
  lmom = new TLegend(0.65,0.80,0.90,0.93);
  lmom->SetTextSize(0.03);
  lmom->SetBorderSize(0);
  lmom->SetHeader(extra_legend.Data(), "C");
  lmom->AddEntry((TObject*)0, Form("%s, %1.1f < #eta < %1.1f", symbolname.Data(), etamin, etamax), "C");
  
  TF1 *func = new TF1("func","gaus",-0.5,0.5);
  
  for (int i =0; i<nfiles; ++i){
    
    TCanvas *cp = new TCanvas("cp","cp",1400,1000);
    cp->SetMargin(0.10, 0.05 ,0.1,0.07);

    //pi-/mom/lfhcal_mom_20.0_mom_resol_pi-.root
    fmom[i] = TFile::Open(Form("./%s/mom/lfhcal_mom_%1.1f_mom_resol_%s.root",particle.Data(),mom[i],particle.Data()));
    
    TH1D *hist = (TH1D*) fmom[i]->Get(Form("hist_mom_%1.1f_%1.1f_pmax_%1.1f",mom[i],etamin,etamax));
    hist->Rebin(2);
    hist->SetName(Form("hist_mom_%1.1f_%1.1f_eta_%1.1f_%s",mom[i],etamin,etamax,particle.Data()));
    hist->SetTitle(Form("Momentum = %1.1f && %1.1f<#eta<%1.1f;#Delta p/p; Entries(a.u.)",particle.Data(),mom[i],etamin,etamax));
    
    double mu = hist->GetMean(); 
    double sigma = hist->GetStdDev();
    hist->GetXaxis()->SetRangeUser(-1.0*range,1.0*range);
    func->SetRange(mu-2.0*sigma,mu+2.0*sigma); // fit with in 2 sigma range
    hist->Fit(func,"NR+");
    mu = func->GetParameter(1); 
    sigma = func->GetParameter(2);
    func->SetRange(mu-2.0*sigma,mu+2.0*sigma);
    hist->Fit(func,"R+");
    float par2 = func->GetParameter(2)*100;
    float par2_err = func->GetParError(2)*100;
    momV.push_back(mom[i]);
    momresolV.push_back(par2);
    err_momresolV.push_back(par2_err);
    
    cp->cd();
    hist->Draw();
    //cp->SaveAs(Form("Debug_Plots/%s/mom/mom_resol_mom%1.1f_%1.1f_eta_%1.1f.png",particle.Data(),mom[i],etamin,etamax));
  } // all files
  
  const int size = momV.size();
  double p[size], err_p[size], sigma_p[size], err_sigma_p[size]; 
  
  for (int i=0; i<size; i++){
    p[i] = momV.at(i);
    sigma_p[i] = momresolV.at(i);
    err_sigma_p[i] = err_momresolV.at(i);
    err_p[i] = 0.;
  }
  
  
  TFile *fout = new TFile(Form("Final_Results/%s/mom/lfhcal_mom_resol_%1.1f_eta_%1.1f.root",particle.Data(),etamin,etamax),"recreate");
  TGraphErrors *gr1 = new TGraphErrors(size,p,sigma_p,err_p,err_sigma_p);
  gr1->SetName("grseed");
  gr1->SetMarkerStyle(25);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerSize(2.0);
  gr1->SetTitle(";p (GeV/c);#sigmap/p");
  gr1->GetXaxis()->CenterTitle();
  gr1->GetYaxis()->CenterTitle();
  
  
  mgMom->Add(gr1);
  c_mom->cd();
  mgMom->GetXaxis()->SetRangeUser(0.40,20.2);
  mgMom->GetYaxis()->SetRangeUser(0.0,1.50*TMath::MaxElement(gr1->GetN(),gr1->GetY())); // 50% more of the maximum value on yaxis
  mgMom->Draw("AP");
  lmom->AddEntry(gr1,"Nominal");
  lmom->Draw("same");
  //draw_req_Mom(etamin,etamax,0.,mgMom->GetXaxis()->GetXmax());
  c_mom->SaveAs(Form("Final_Results/%s/mom/lfhcal_mom_resol_%1.1f_eta_%1.1f.png",particle.Data(),etamin,etamax));
  
  // Write the numbers in output file for comparisons
  outfile << extra_legend << endl;
  outfile<<"Etamin"<<setw(20)<<"Etamax"<<setw(20)<<"p (GeV/c) \t"<<setw(20)<<"Resol  #mum "<<endl;
  for (Int_t i = 0; i<gr1->GetN(); ++i){
    double x,y;
    gr1->GetPoint(i,x,y);
    outfile<<etamin<<setw(20)<<etamax<<setw(20)<<x<<setw(20)<<y<<endl;
  }
  outfile.close();
  
  fout->cd();
  mgMom->SetName(Form("mom_resol_%1.1f_eta_%1.1f",etamin,etamax));
  mgMom->Write();
  fout->Close();
}

//===Fit Momentum Resolution
float FitMomentumResolution(Double_t *x, Double_t *par)
{
  float func = sqrt(par[0]*par[0]*x[0]*x[0]+par[1]*par[1]);
  return func;
}

//From Yellow report from section 11.2.2

// 1.2,1.5,2,2.5,3,3.5

/*
void draw_req_Mom(double etamin, double etamax, double xmin=0., double xmax=0.)
{

   TF1 *dd4hep_p;
   if (etamin >= -3.5 && etamax <= -2.5) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.1*x)^2+2.0^2)",xmin,xmax);
   else if (etamin >= -2.5 && etamax <= -1.0) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.05*x)^2+1.0^2)",xmin,xmax);
   else if (etamin >= -1.0 && etamax <= 1.0) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.05*x)^2+0.5^2)",xmin,xmax);
   else if (etamin >= 1.0 && etamax <= 2.5) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.05*x)^2+1.0^2)",xmin,xmax);
   else if (etamin >= 2.5 && etamax <= 3.5) dd4hep_p = new TF1("dd4hep_p", "TMath::Sqrt((0.1*x)^2+2.0^2)",xmin,xmax);
   else return;
   dd4hep_p->SetLineStyle(7);
   dd4hep_p->SetLineColor(kMagenta);
   dd4hep_p->SetLineWidth(3.0);
   dd4hep_p->Draw("same");

  TLegend *l= new TLegend(0.70,0.75,0.90,0.80);
  l->SetTextSize(0.03);
  l->SetBorderSize(0);
  l->AddEntry(dd4hep_p,"PWGReq","l");
  l->Draw("same");
 }
*/
