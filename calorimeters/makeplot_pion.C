//////////////////////////////////////
// Read ROOT file and plot variables
//////////////////////////////////////

int makeplot_pion(void)
{
  // Setting figures
  gROOT->SetStyle("Plain");
  gStyle->SetLineWidth(3);
  gStyle->SetOptStat("nem");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetPadLeftMargin(0.14);

  // Input ROOT file
  TFile *f = new TFile("/home/jihee/topside/rec_crystal_pion_output.root","read");
  TTree *t = (TTree *)f->Get("events");

  // Set Branch status and addressed
  t->SetMakeClass(1);
  t->SetBranchStatus("*", 0);

  Int_t EcalClusters_;
  t->SetBranchStatus("EcalClusters", 1);
  t->SetBranchAddress("EcalClusters", &EcalClusters_);

  const Int_t kMaxEcalClusters = 4;
  Double_t cluster_x_pos[kMaxEcalClusters];
  Double_t cluster_y_pos[kMaxEcalClusters];
  Double_t cluster_z_pos[kMaxEcalClusters];
  Float_t cluster_energy[kMaxEcalClusters];
  t->SetBranchStatus("EcalClusters.position.x",1);
  t->SetBranchStatus("EcalClusters.position.y",1);
  t->SetBranchStatus("EcalClusters.position.z",1);
  t->SetBranchStatus("EcalClusters.energy",1);
  t->SetBranchAddress("EcalClusters.position.x",cluster_x_pos);
  t->SetBranchAddress("EcalClusters.position.y",cluster_y_pos);
  t->SetBranchAddress("EcalClusters.position.z",cluster_z_pos);
  t->SetBranchAddress("EcalClusters.energy",cluster_energy);

  // Setting for Canvas
  TCanvas *c1 = new TCanvas("c1","c1", 600, 600);
  TCanvas *c2 = new TCanvas("c2","c2", 600, 600);
  TCanvas *c3 = new TCanvas("c3","c3", 600, 600);
  TCanvas *c4 = new TCanvas("c4","c4", 600, 600);
  
  TCanvas *c6 = new TCanvas("c6","c6", 600, 600);
  TCanvas *c7 = new TCanvas("c7","c7", 600, 600);

  TCanvas *c8 = new TCanvas("c8","c8", 600, 600);
  TCanvas *c9 = new TCanvas("c9","c9", 600, 600);

  // Declare histograms
  TH1D *h1 = new TH1D("Scattering angle","Scattering Angle(#theta)",90,135.0,180.0);
  TH1D *h2 = new TH1D("Pseudo-rapidity","Pseudo-rapidity(#eta)",50,-5.0,0.0);
  TH2D *h3 = new TH2D("E vs #eta","Cluster E vs Pseudo-rapidity",100,0.0,1.0,50,-5.0,0.0);
  TH1D *h4 = new TH1D("Reconstructed E","Reconstructed energy per event",100,0.0,1.0);
  TH1D *h5 = new TH1D("Thrown E","Thrown energy per event",100,0.0,1.0);
  TH2D *h6 = new TH2D("theta vs #eta","Scattering angle(#theta) vs. Pseudo-rapidity",90,135.0,180.0,50,-5.0,0.0);
  TH1D *h7 = new TH1D("Invariant mass","Invariant mass",30,0.0,300.0);

  TH1D *h8 = new TH1D("E1","E1",100,0.0,1000.0);
  TH1D *h9 = new TH1D("E2","E2",100,0.0,1000.0);
  TH1D *h10 = new TH1D("angle", "angle", 100,0.0,180.0);

  // Total number of entries
  Int_t nentries = t->GetEntries();

  // Variables are used in calculation
  Double_t r;                        // Radius [cm]
  Double_t phi;                      // Azimuth [degree]
  Double_t theta;                    // Inclination [degree]
  Double_t eta;                      // Pseudo-rapidity [unitless]
  Float_t  cluster_e;                // Cluster energy [GeV]
  Float_t  total_cluster_e;          // Add up clusters per event [GeV]
  Double_t dot_product_pos_clusters; // dot product of positions of two photons
  Double_t mag_pos2_cluster_1;       // squared magnitude of position
  Double_t mag_pos2_cluster_2;       // squared magnitude of position 
  Double_t cosine_clusters;          // cos(theta_photons)
  Double_t theta_photons;            // angle between two photons
  Double_t invariant_mass;           // M^2 = 2 * p_1 * p_2 * (1 - cos(theta_photons))

  // Loop over event by event
  for (int ievent = 0; ievent < nentries; ievent++)
  {
	t->GetEntry(ievent);

	Int_t ncluster = EcalClusters_;

	total_cluster_e = 0.0;

	// Loop over cluster by cluster
	for (int icluster=0; icluster < ncluster; icluster++)
	{
		r = TMath::Sqrt((cluster_x_pos[icluster]*cluster_x_pos[icluster]) + 
				(cluster_y_pos[icluster]*cluster_y_pos[icluster]) + 
				(cluster_z_pos[icluster]*cluster_z_pos[icluster]));
		phi = TMath::ATan(cluster_y_pos[icluster]/cluster_x_pos[icluster]) * TMath::RadToDeg();
		theta = TMath::ACos(cluster_z_pos[icluster] / r) * TMath::RadToDeg();
		eta = -1.0 * TMath::Log(TMath::Tan((theta*TMath::DegToRad())/2.0));	
		cluster_e = cluster_energy[icluster] / 1.e+3;
		total_cluster_e += cluster_e;

		// Fill histograms
		h1->Fill(theta, 1.0);
		h2->Fill(eta, 1.0);
		h3->Fill(cluster_e, eta, 1.0);
		h6->Fill(theta, eta, 1.0);
	}
	if(ncluster > 0)
		h4->Fill(total_cluster_e, 1.0);
 
	// Find events with 2 clusters
	// To calculate invariant mass
	// M^2 = 2p1p2(1-cos(theta))
	// p1 = E1
	// p2 = E2
	// theta: angle between two photons	
	if(ncluster == 2)
	{
		dot_product_pos_clusters = cluster_x_pos[0]*cluster_x_pos[1] + cluster_y_pos[0]*cluster_y_pos[1];
		mag_pos2_cluster_1 = (cluster_x_pos[0]*cluster_x_pos[0]) + (cluster_y_pos[0]*cluster_y_pos[0]);
		mag_pos2_cluster_2 = (cluster_x_pos[1]*cluster_x_pos[1]) + (cluster_y_pos[1]*cluster_y_pos[1]);
		cosine_clusters = (dot_product_clusters/TMath::Sqrt(mag_cluster_1*mag_cluster_2));
		theta_photons = TMath::Acos(cosine_clusters)*TMath::RadToDeg();
		
		invariant_mass = TMath::Sqrt(2.0*cluster_energy[0]*cluster_energy[1]*(1.0 - cosine_clusters));
		
		// Fill histograms
		h7->Fill(invariant_mass, 1.0);
		h8->Fill(cluster_energy[0], 1.0);
		h9->Fill(cluster_energy[1], 1.0);
		h10->Fill(theta_photons, 1.0);
	}

  }

  // Drawing and Saving figures
  c1->cd();
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2);
  h1->GetXaxis()->SetTitle("#theta [degree]");
  h1->GetYaxis()->SetTitle("events");
  h1->GetYaxis()->SetTitleOffset(1.4);
  gPad->Update();
  h1->Draw();
  //h1->SaveAs("./plots/theta_hist.png");

  c2->cd();
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);
  h2->GetXaxis()->SetTitle("#eta");
  h2->GetYaxis()->SetTitle("events");
  h2->GetYaxis()->SetTitleOffset(1.4);
  h2->Draw();
  //h2->SaveAs("./plots/eta_hist.png");

  c3->cd();
  h3->GetXaxis()->SetTitle("Cluster energy [GeV]");
  h3->GetYaxis()->SetTitle("#eta");
  h3->GetYaxis()->SetTitleOffset(1.4);
  h3->Draw("COLZ");
  //h3->SaveAs("./plots/e_vs_eta_hist.png");

  c4->cd();
  c4->SetLogy(1);
  h4->SetLineColor(kBlue);
  h4->SetLineWidth(2);
  h4->GetXaxis()->SetTitle("reconstructed energy [GeV]");
  h4->GetYaxis()->SetTitle("events");
  h4->GetYaxis()->SetTitleOffset(1.4);
  h4->Draw();
  //h4->SaveAs("./plots/recon_e_hist.png");

  c6->cd();
  h6->GetXaxis()->SetTitle("#theta [degree]");
  h6->GetYaxis()->SetTitle("#eta");
  h6->GetYaxis()->SetTitleOffset(1.4);
  h6->Draw("COLZ");
  //h6->SaveAs("./plots/theta_vs_eta_hist.png");

  c7->cd();
  h7->SetLineColor(kBlue);
  h7->SetLineWidth(2);
  h7->GetXaxis()->SetTitle("Invariant mass [MeV]");
  h7->GetYaxis()->SetTitle("events");
  h7->GetYaxis()->SetTitleOffset(1.4);
  h7->Draw();

  c8->cd();
  h8->SetLineColor(kBlue);
  h8->SetLineWidth(2);
  h8->GetXaxis()->SetTitle("Cluster energy 1 [MeV]");
  h8->GetYaxis()->SetTitle("events");
  h8->GetYaxis()->SetTitleOffset(1.4);
  h8->Draw();

  c9->cd();
  h9->SetLineColor(kBlue);
  h9->SetLineWidth(2);
  h9->GetXaxis()->SetTitle("Cluster energy 2 [MeV]");
  h9->GetYaxis()->SetTitle("events");
  h9->GetYaxis()->SetTitleOffset(1.4);
  h9->Draw();

  c10->cd();
  h10->SetLineColor(kBlue);
  h10->SetLineWidth(2);
  h10->GetXaxis()->SetTitle("angle between two photons [degree]");
  h10->GetYaxis()->SetTitle("events");
  h10->GetYaxis()->SetTitleOffset(1.4);
  h10->Draw();

  return 0;
}
