//------------------
void fwd_neutrons_recon(std::string inputfile, std::string outputfile){

    //Define Style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);

    //Define histograms
    TH2 *h1_neut = new TH2D("h1_neut","Neutron true energy vs. polar angle",100,0.6,2.2,100,0,200);
    h1_neut->GetXaxis()->SetTitle("#theta [deg]"); h1_neut->GetXaxis()->CenterTitle();
    h1_neut->GetYaxis()->SetTitle("E [GeV]"); h1_neut->GetYaxis()->CenterTitle();

    TH2 *h2_neut = new TH2D("h2_neut","Neutron true azimuthal angle vs. polar angle around p axis",100,-0.1,12,100,-200,200);
    h2_neut->GetXaxis()->SetTitle("#theta^{*} [mRad]"); h2_neut->GetXaxis()->CenterTitle();
    h2_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h2_neut->GetYaxis()->CenterTitle();

    TH2 *h2a_neut = new TH2D("h2a_neut","Neutron true azimuthal angle vs. cosine of polar angle around p axis",100,0.99996,1,100,-200,200);
    h2a_neut->GetXaxis()->SetTitle("cos(#theta^{*})"); h2a_neut->GetXaxis()->CenterTitle();
    h2a_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h2a_neut->GetYaxis()->CenterTitle();

    TH1 *h1_ecal_adc = new TH1D("h1_ecal_adc","ECal ADC amplitude spectrum",1000,-0.1,35000);
    h1_ecal_adc->GetXaxis()->SetTitle("ADC Channel"); h1_ecal_adc->GetXaxis()->CenterTitle();
    h1_ecal_adc->SetLineColor(kRed);h1_ecal_adc->SetLineWidth(2);

    TH1 *h1_hcal_adc = new TH1D("h1_hcal_adc","HCal ADC amplitude spectrum",1000,-0.1,35000);
    h1_hcal_adc->GetXaxis()->SetTitle("ADC Channel"); h1_hcal_adc->GetXaxis()->CenterTitle();
    h1_hcal_adc ->SetLineColor(kBlue);h1_hcal_adc->SetLineWidth(2);

    TH1 *h1_ecal = new TH1D("h1_ecal","Total reconstructed hit energy sum",100,-0.1,2);
    h1_ecal->GetXaxis()->SetTitle("E [GeV]"); h1_ecal->GetXaxis()->CenterTitle();
    h1_ecal->SetLineColor(kRed);h1_ecal->SetLineWidth(2);

    TH1 *h1_hcal = new TH1D("h1_hcal","Total reconstructed hit energy sum",100,-0.1,2);
    h1_hcal->GetXaxis()->SetTitle("E [GeV]"); h1_hcal->GetXaxis()->CenterTitle();
    h1_hcal->SetLineColor(kBlue);h1_hcal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h2_ecal = new TH1D("h2_ecal","Total reconstructed hit energy sum",100,-0.1,2);
    h2_ecal->GetXaxis()->SetTitle("E [GeV]"); h2_ecal->GetXaxis()->CenterTitle();
    h2_ecal->SetLineColor(kRed);h2_ecal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h2_hcal = new TH1D("h2_hcal","Total reconstructed hit energy sum",100,-0.1,2);
    h2_hcal->GetXaxis()->SetTitle("E [GeV]"); h2_hcal->GetXaxis()->CenterTitle();
    h2_hcal->SetLineColor(kBlue);h2_hcal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h3_ecal = new TH1D("h3_ecal","Number of reconstructed clusters in ECal",5,0,5);
    h3_ecal->GetXaxis()->SetTitle("Number of clusters"); h3_ecal->GetXaxis()->CenterTitle();
    h3_ecal->GetXaxis()->SetNdivisions(405);h3_ecal->GetXaxis()->CenterLabels();
    h3_ecal->SetLineColor(kRed);h3_ecal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h3_hcal = new TH1D("h3_hcal","Number of reconstructed clusters in HCal",5,0,5);
    h3_hcal->GetXaxis()->SetTitle("Number of clusters"); h3_hcal->GetXaxis()->CenterTitle();
    h3_hcal->GetXaxis()->SetNdivisions(405);h3_hcal->GetXaxis()->CenterLabels();
    h3_hcal->SetLineColor(kBlue);h3_hcal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h4_hcal = new TH1D("h4_hcal","HCal cluster energy: Events w/ 1 HCal Cluster",100,-0.1,120);
    h4_hcal->GetXaxis()->SetTitle("E [GeV]"); h4_hcal->GetXaxis()->CenterTitle();
    h4_hcal->SetLineColor(kBlue);h4_hcal->SetLineWidth(2);

    //Cut on generated neutron theta* + 1 HCal Cluster
    TH1 *h1_neut_rec = new TH1D("h1_neut_rec","Reconstructed Neutron Energy",100,-0.1,120);
    h1_neut_rec->GetXaxis()->SetTitle("E [GeV]"); h1_neut_rec->GetXaxis()->CenterTitle();
    h1_neut_rec->SetLineColor(kBlue);h1_neut_rec->SetLineWidth(2);

    //Cut on generated neutron theta* + 1 HCal Cluster
    TH1 *h2_neut_rec = new TH2D("h2_neut_rec","Reconstructed Neutron local hit position at ZDC HCal front face",100,-200,200,100,-200,200);
    h2_neut_rec->GetXaxis()->SetTitle("x [mm]"); h2_neut_rec->GetXaxis()->CenterTitle();
    h2_neut_rec->GetYaxis()->SetTitle("y [mm]"); h2_neut_rec->GetYaxis()->CenterTitle();

    //Read ROOT file
    TFile* file = new TFile(inputfile.c_str());
    TTree *tree = (TTree*) file->Get("events"); 
 
    cout<<"Total number of events to analyze is "<<tree->GetEntries()<<endl;

    //Create Array Reader
    TTreeReader tr(tree);

    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<double> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<double> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<double> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass");
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");

    //ZDC LYSO ECal
    TTreeReaderArray<int> ecal_adc(tr,"EcalFarForwardZDCRawHits.amplitude");
    TTreeReaderArray<float> ecal_hit_e(tr,"EcalFarForwardZDCRecHits.energy");
    TTreeReaderArray<float> ecal_cluster_energy(tr, "EcalFarForwardZDCClusters.energy");

    //ZDC SiPM-on-tile HCal
    TTreeReaderArray<int> hcal_adc(tr,"HcalFarForwardZDCRawHits.amplitude");
    TTreeReaderArray<float> hcal_hit_e(tr, "HcalFarForwardZDCRecHits.energy");
    TTreeReaderArray<float> hcal_cluster_energy(tr, "HcalFarForwardZDCClusters.energy");
    
    //Reconstructed neutron quantity
    TTreeReaderArray<float> rec_neutron_energy(tr,"ReconstructedFarForwardZDCNeutrals.energy");
    TTreeReaderArray<float> rec_neutron_px(tr,"ReconstructedFarForwardZDCNeutrals.momentum.x");
    TTreeReaderArray<float> rec_neutron_py(tr,"ReconstructedFarForwardZDCNeutrals.momentum.y");
    TTreeReaderArray<float> rec_neutron_pz(tr,"ReconstructedFarForwardZDCNeutrals.momentum.z");
    TTreeReaderArray<int> rec_neutron_PDG(tr,"ReconstructedFarForwardZDCNeutrals.PDG");
    //Other variables
    int counter(0);
    
    TLorentzVector neut_true; //True neutron in lab coordinates
    TLorentzVector neut_true_rot; //True neutron wrt proton beam direction

    float ecal_e_tot(0);
    float hcal_e_tot(0);
    
    //Loop over events
    while (tr.Next()) {
	
	if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	counter++;

       //Reset variables
       ecal_e_tot = 0;
       hcal_e_tot = 0;

       //Loop over generated particles, select primary neutron
       for(int igen=0;igen<gen_status.GetSize();igen++){
        	if(gen_status[igen]==1 && gen_pid[igen]==2112){

                     neut_true.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                     h1_neut->Fill(neut_true.Theta()*TMath::RadToDeg(),neut_true.E());

			//Wrt proton beam direction
			neut_true_rot = neut_true;
			neut_true_rot.RotateY(0.025);
			h2_neut->Fill(neut_true_rot.Theta()*1000.,neut_true_rot.Phi()*TMath::RadToDeg()); //Theta* in mRad
                        h2a_neut->Fill(std::cos(neut_true_rot.Theta()),neut_true_rot.Phi()*TMath::RadToDeg());

              }
       } //End loop over generated particles

       //Loop over ECal ADC hits (Raw hits may have different length than Rec hits due to zero supression)
       for(int iadc=0;iadc<ecal_adc.GetSize();iadc++){
               h1_ecal_adc->Fill(ecal_adc[iadc]);  
       }

       //Loop over ECal hits
       for(int ihit=0;ihit<ecal_hit_e.GetSize();ihit++){
               ecal_e_tot +=  ecal_hit_e[ihit];           
       }

       //Loop over ECal ADC hits (Raw hits may have different length than Rec hits due to zero supression)
       for(int iadc=0;iadc<hcal_adc.GetSize();iadc++){
               h1_hcal_adc->Fill(hcal_adc[iadc]);            
       }

       //Loop over HCal hits
       for(int ihit=0;ihit<hcal_hit_e.GetSize();ihit++){
               hcal_e_tot +=  hcal_hit_e[ihit];           
       }

       //ECal cluster size
       int ecal_clus_size = ecal_cluster_energy.GetSize();

       //HCal cluster size
       int hcal_clus_size = hcal_cluster_energy.GetSize();

       //Fill histograms for total energy and clusters
       h1_ecal->Fill(ecal_e_tot);
       h1_hcal->Fill(hcal_e_tot);

       if( neut_true_rot.Theta()*1000. < 3.5 ){
                h2_ecal->Fill(ecal_e_tot);
                h2_hcal->Fill(hcal_e_tot);

                h3_ecal->Fill(ecal_clus_size);
                h3_hcal->Fill(hcal_clus_size);

                //HCal cluster energy -- 1 cluster events
                if(hcal_clus_size==1) h4_hcal->Fill(hcal_cluster_energy[0]);
       }

       //Reconstructed neutron(s)
       for(int irec=0;irec<rec_neutron_energy.GetSize();irec++){
	 if (rec_neutron_PDG[irec] != 2112)
	   continue;
                if(neut_true_rot.Theta()*1000. < 3.5 && hcal_clus_size==1){
                        h1_neut_rec->Fill(rec_neutron_energy[irec]);

                        TVector3 neut_rec(rec_neutron_px[irec],rec_neutron_py[irec],rec_neutron_pz[irec]);
                        
                        //Wrt proton beam direction
			TVector3 neut_rec_rot = neut_rec;
			neut_rec_rot.RotateY(0.025);

                        //Projection of reconstructed neutron to front face of HCal
                        auto proj_x = 35.8 * 1000. * tan(neut_rec_rot.Theta()) * cos(neut_rec_rot.Phi()) ;
                        auto proj_y = 35.8 * 1000. * tan(neut_rec_rot.Theta()) * sin(neut_rec_rot.Phi()) ;
                        h2_neut_rec->Fill(proj_x,proj_y);
                }          
       } //End loop over reconstructed neutrons

    } //End loop over events

    //Make plots
    TCanvas *c1 = new TCanvas("c1");
    h1_neut->Draw("colz");

    TCanvas *c2 = new TCanvas("c2");
    h2_neut->Draw("colz");

    TCanvas *c2a = new TCanvas("c2a");
    h2a_neut->Draw("colz");

    TCanvas *c3a = new TCanvas("c3a");
    c3a->SetLogy();
    h1_ecal_adc->Draw();

    TCanvas *c3b = new TCanvas("c3b");
    c3b->SetLogy();
    h1_hcal_adc->Draw();

    TCanvas *c4 = new TCanvas("c4");
    h1_ecal->Draw("");
    h1_hcal->Draw("same");

    TLegend *leg4 = new TLegend(0.6,0.6,0.85,0.8);
    leg4->SetBorderSize(0);leg4->SetFillStyle(0);
    leg4->AddEntry(h1_ecal,"Sum of digitized ZDC Ecal hit energies","l");
    leg4->AddEntry(h1_hcal,"Sum of digitized ZDC Hcal hit energies","l");
    leg4->Draw();

    TCanvas *c5 = new TCanvas("c5");
    h2_ecal->Draw("");
    h2_hcal->Draw("same");

    TLegend *leg5 = new TLegend(0.5,0.6,0.9,0.8);
    leg5->SetBorderSize(0);leg5->SetFillStyle(0);
    leg5->SetHeader("Require neutron #theta^{*} (wrt proton beam) < 3.5 mRad");
    leg5->AddEntry(h1_ecal,"Sum of digitized ZDC Ecal hit energies","l");
    leg5->AddEntry(h1_hcal,"Sum of digitized ZDC Hcal hit energies","l");
    leg5->Draw();

    TCanvas *c6a = new TCanvas("c6a");
    h3_ecal->Draw();

    TLegend *leg6 = new TLegend(0.5,0.7,0.9,0.9);
    leg6->SetBorderSize(0);leg6->SetFillStyle(0);
    leg6->SetHeader("Require neutron #theta^{*} (wrt proton beam) < 3.5 mRad");
    leg6->Draw();

    TCanvas *c6b = new TCanvas("c6b");
    h3_hcal->Draw();
    leg6->Draw();

    TCanvas *c7 = new TCanvas("c7");
    h4_hcal->Draw();

    TLegend *leg7 = new TLegend(0.15,0.7,0.5,0.9);
    leg7->SetBorderSize(0);leg7->SetFillStyle(0);
    leg7->SetHeader("Require neutron #theta^{*} (wrt proton beam) < 3.5 mRad");
    leg7->Draw();

    TCanvas *c8 = new TCanvas("c8");
    h1_neut_rec->Draw();

    TLegend *leg8 = new TLegend(0.15,0.7,0.7,0.9);
    leg8->SetBorderSize(0);leg8->SetFillStyle(0);
    leg8->SetHeader("Generated neutron #theta^{*} < 3.5 mRad + 1 HCal cluster");
    leg8->Draw();

    TCanvas *c9 = new TCanvas("c9");
    h2_neut_rec->Draw("colz");

    TLatex *tex9 = new TLatex(-150,150,"Generated neutron #theta^{*} < 3.5 mRad + 1 HCal cluster");
    tex9->SetTextSize(0.03);
    tex9->Draw();

    //Print plots to file
    c1->Print(Form("%s[",outputfile.c_str()));
    c1->Print(Form("%s",outputfile.c_str()));
    c2->Print(Form("%s",outputfile.c_str()));
    c2a->Print(Form("%s",outputfile.c_str()));
    c3a->Print(Form("%s",outputfile.c_str()));
    c3b->Print(Form("%s",outputfile.c_str()));
    c4->Print(Form("%s",outputfile.c_str()));
    c5->Print(Form("%s",outputfile.c_str()));
    c6a->Print(Form("%s",outputfile.c_str()));
    c6b->Print(Form("%s",outputfile.c_str()));
    c7->Print(Form("%s",outputfile.c_str()));
    c8->Print(Form("%s",outputfile.c_str()));
    c9->Print(Form("%s",outputfile.c_str()));
    c9->Print(Form("%s]",outputfile.c_str())); 

}
