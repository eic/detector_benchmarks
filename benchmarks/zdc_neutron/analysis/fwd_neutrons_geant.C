//------------------
int get_layer_number(TVector3 pos){

       // Get Layer position wrt proton beam direction
	pos.RotateY(0.025);

       auto local_z = pos.Z() - 35.8*1000.; //in mm
       
       // Convert to layer number
       // local_z = 22.25 + (layer_number - 1)*24.9
       int layer_number = (int)std::round( (local_z - 22.25)/24.9 ) + 1;

       return layer_number;
}

//------------------
void fwd_neutrons_geant(std::string inputfile, std::string outputfile){

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

    //Require HCal hit energy sum > 1 GeV
    TH2 *h3_neut = new TH2D("h3_neut","Neutron true azimuthal angle vs. polar angle around p axis",100,-0.1,12,100,-200,200);
    h3_neut->GetXaxis()->SetTitle("#theta^{*} [mRad]"); h3_neut->GetXaxis()->CenterTitle();
    h3_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h3_neut->GetYaxis()->CenterTitle();

    TH2 *h4_neut = new TH2D("h4_neut","Neutron local hit position at ZDC HCal front face",100,-300,300,100,-300,300);
    h4_neut->GetXaxis()->SetTitle("x [mm]"); h4_neut->GetXaxis()->CenterTitle();
    h4_neut->GetYaxis()->SetTitle("y [mm]"); h4_neut->GetYaxis()->CenterTitle();

    TH1 *h1_ecal = new TH1D("h1_ecal","Total true hit energy sum",100,-0.1,2);
    h1_ecal->GetXaxis()->SetTitle("E [GeV]"); h1_ecal->GetXaxis()->CenterTitle();
    h1_ecal->SetLineColor(kRed);h1_ecal->SetLineWidth(2);

    TH1 *h1_hcal = new TH1D("h1_hcal","Total true hit energy sum",100,-0.1,2);
    h1_hcal->GetXaxis()->SetTitle("E [GeV]"); h1_hcal->GetXaxis()->CenterTitle();
    h1_hcal->SetLineColor(kBlue);h1_hcal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h2_ecal = new TH1D("h2_ecal","Total true hit energy sum",100,-0.1,2);
    h2_ecal->GetXaxis()->SetTitle("E [GeV]"); h2_ecal->GetXaxis()->CenterTitle();
    h2_ecal->SetLineColor(kRed);h2_ecal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h2_hcal = new TH1D("h2_hcal","Total true hit energy sum",100,-0.1,2);
    h2_hcal->GetXaxis()->SetTitle("E [GeV]"); h2_hcal->GetXaxis()->CenterTitle();
    h2_hcal->SetLineColor(kBlue);h2_hcal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h2a_ecal = new TH1D("h2a_ecal","Total true hit energy sum",200,-0.1,30);
    h2a_ecal->GetXaxis()->SetTitle("E [GeV]"); h2a_ecal->GetXaxis()->CenterTitle();
    h2a_ecal->SetLineColor(kRed);h2a_ecal->SetLineWidth(2);

    //Cut on theta*
    TH1 *h2a_hcal = new TH1D("h2a_hcal","Total true hit energy sum",200,-0.1,30);
    h2a_hcal->GetXaxis()->SetTitle("E [GeV]"); h2a_hcal->GetXaxis()->CenterTitle();
    h2a_hcal->SetLineColor(kBlue);h2a_hcal->SetLineWidth(2);

    //Cut on theta*
    TH2 *h2_cal_both = new TH2D("h2_cal_both","Total true hit energy sum: ECal vs. HCal",100,-0.1,2.5,100,-0.1,30);
    h2_cal_both->GetXaxis()->SetTitle("E_{HCal.} [GeV]"); h2_cal_both->GetXaxis()->CenterTitle();
    h2_cal_both->GetYaxis()->SetTitle("E_{ECal.} [GeV]"); h2_cal_both->GetYaxis()->CenterTitle();

    TH1 *h3a_hcal = new TH1D("h3a_hcal","Cells w/ non-zero energy deposit ",100,0,64);
    h3a_hcal->GetXaxis()->SetTitle("Layer Number in ZDC HCal"); h3a_hcal->GetXaxis()->CenterTitle();
    h3a_hcal->GetYaxis()->SetTitle("Number of cells hit per event");h3a_hcal->GetYaxis()->CenterTitle();
    h3a_hcal->SetLineColor(kBlue);h3a_hcal->SetLineWidth(3);

    //Cut on theta* < 3.5mRad
    TH1 *h3b_hcal = new TH1D("h3b_hcal","Cells w/ non-zero energy deposit ",100,0,64);
    h3b_hcal->GetXaxis()->SetTitle("Layer Number in ZDC HCal"); h3b_hcal->GetXaxis()->CenterTitle();
    h3b_hcal->GetYaxis()->SetTitle("Number of cells hit per event");h3b_hcal->GetYaxis()->CenterTitle();
    h3b_hcal->SetLineColor(kBlue);h3b_hcal->SetLineWidth(3);

    //Cut on 6mRad < theta* < 8mRad
    TH1 *h3c_hcal = new TH1D("h3c_hcal","Cells w/ non-zero energy deposit ",100,0,64);
    h3c_hcal->GetXaxis()->SetTitle("Layer Number in ZDC HCal"); h3c_hcal->GetXaxis()->CenterTitle();
    h3c_hcal->GetYaxis()->SetTitle("Number of cells hit per event");h3c_hcal->GetYaxis()->CenterTitle();
    h3c_hcal->SetLineColor(kBlue);h3c_hcal->SetLineWidth(3);

    //Read ROOT file
    TFile* file = new TFile(inputfile.c_str());
    TTree *tree = (TTree*) file->Get("events");

    //Set cut
    float edep_cut = 1.; //In GeV
    
    cout<<"Total number of events to analyze is "<<tree->GetEntries()<<endl;

    //Create Array Reader
    TTreeReader tr(tree);

    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<float> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<float> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<float> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass");
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");

    //ZDC LYSO ECal
    TTreeReaderArray<float> ecal_hit_e(tr,"EcalFarForwardZDCHits.energy");

    //ZDC SiPM-on-tile HCal
    TTreeReaderArray<float> hcal_hit_e(tr, "HcalFarForwardZDCHits.energy");
    TTreeReaderArray<float> hcal_hit_x(tr, "HcalFarForwardZDCHits.position.x");
    TTreeReaderArray<float> hcal_hit_y(tr, "HcalFarForwardZDCHits.position.y");
    TTreeReaderArray<float> hcal_hit_z(tr, "HcalFarForwardZDCHits.position.z");

    //Other variables
    int counter(0),counter_3p5(0),counter_6to8(0);
    
    TLorentzVector neut_true; //True neutron in lab coordinates
    TLorentzVector neut_true_rot; //True neutron wrt proton beam direction

    float ecal_e_tot(0);
    float hcal_e_tot(0);
    
    //Loop over events
    while (tr.Next()) {
	
	if(counter%1000==0) cout<<"Analyzing event "<<counter<<endl;
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

       //Loop over ECal hits
       for(int ihit=0;ihit<ecal_hit_e.GetSize();ihit++){
               ecal_e_tot +=  ecal_hit_e[ihit];           
       }

       //Loop over HCal hits
       for(int ihit=0;ihit<hcal_hit_e.GetSize();ihit++){
              hcal_e_tot +=  hcal_hit_e[ihit];

              //Cell hit position in global coordinates
              double hit_x = (double) hcal_hit_x[ihit];
              double hit_y = (double) hcal_hit_y[ihit];
              double hit_z = (double) hcal_hit_z[ihit];
              
              TVector3 hit_pos_det(hit_x,hit_y,hit_z);
              auto layer_number = get_layer_number(hit_pos_det);
              
              h3a_hcal->Fill(layer_number);
              if( neut_true_rot.Theta()*1000. < 3.5 ) h3b_hcal->Fill(layer_number);
              if( neut_true_rot.Theta()*1000. > 6 &&  neut_true_rot.Theta()*1000 < 8) h3c_hcal->Fill(layer_number);
              
       } //End loop over HCal hits

       //Fill histograms
       h1_ecal->Fill(ecal_e_tot);
       h1_hcal->Fill(hcal_e_tot);

       if(hcal_e_tot > edep_cut){ 
              h3_neut->Fill(neut_true_rot.Theta()*1000.,neut_true_rot.Phi()*TMath::RadToDeg()); //Theta* in mRad

              //Projection to front face of HCal
              auto proj_x = 35.8 * 1000. * tan(neut_true_rot.Theta()) * cos(neut_true_rot.Phi()) ;
              auto proj_y = 35.8 * 1000. * tan(neut_true_rot.Theta()) * sin(neut_true_rot.Phi()) ;

              h4_neut->Fill(proj_x,proj_y);
       }

       if( neut_true_rot.Theta()*1000. < 3.5 ){
              h2_ecal->Fill(ecal_e_tot);
              h2_hcal->Fill(hcal_e_tot);
              h2a_ecal->Fill(ecal_e_tot);
              h2a_hcal->Fill(hcal_e_tot);

              h2_cal_both->Fill(hcal_e_tot,ecal_e_tot);

              counter_3p5++;
       }
       
       if( neut_true_rot.Theta()*1000. > 6 &&  neut_true_rot.Theta()*1000 < 8){
              counter_6to8++;       
       }

    } //End loop over events

    //Scale histograms
    h3a_hcal->Scale(1./counter);
    h3b_hcal->Scale(1./counter_3p5);
    h3c_hcal->Scale(1./counter_6to8);

    //Make plots
    TCanvas *c1 = new TCanvas("c1");
    h1_neut->Draw("colz");

    TCanvas *c2 = new TCanvas("c2");
    h2_neut->Draw("colz");

    TCanvas *c2a = new TCanvas("c2a");
    h2a_neut->Draw("colz");

    TCanvas *c3 = new TCanvas("c3");
    h1_ecal->Draw("");
    h1_hcal->Draw("same");

    TLegend *leg3 = new TLegend(0.6,0.6,0.85,0.8);
    leg3->SetBorderSize(0);leg3->SetFillStyle(0);
    leg3->AddEntry(h1_ecal,"Sum of ZDC Ecal hit energies","l");
    leg3->AddEntry(h1_hcal,"Sum of ZDC Hcal hit energies","l");
    leg3->Draw();

    TCanvas *c4 = new TCanvas("c4");
    h3_neut->Draw("colz");

    TLatex *tex4 = new TLatex(7,150,Form("Hcal hit energy sum > %.1f GeV",edep_cut));
    tex4->SetTextSize(0.035);
    tex4->Draw();

    TCanvas *c4a = new TCanvas("c4a");
    h4_neut->Draw("colz");

    TLatex *tex4a = new TLatex(50,250,Form("Hcal hit energy sum > %.1f GeV",edep_cut));
    tex4a->SetTextSize(0.035);
    tex4a->Draw();

    TCanvas *c5 = new TCanvas("c5");
    h2_ecal->Draw("");
    h2_hcal->Draw("same");

    TLegend *leg5 = new TLegend(0.5,0.6,0.9,0.8);
    leg5->SetBorderSize(0);leg5->SetFillStyle(0);
    leg5->SetHeader("Require neutron #theta^{*} (wrt proton beam) < 3.5 mRad");
    leg5->AddEntry(h1_ecal,"Sum of ZDC Ecal hit energies","l");
    leg5->AddEntry(h1_hcal,"Sum of ZDC Hcal hit energies","l");
    leg5->Draw();

    TCanvas *c5a = new TCanvas("c5a");
    c5a->SetLogy();
    h2a_ecal->Draw("");
    h2a_hcal->Draw("same");

    leg5->Draw();

    TCanvas *c5b = new TCanvas("c5b");
    c5b->SetLogz();
    h2_cal_both->Draw("colz");

    TLatex *tex5b = new TLatex(0,28,"Require neutron #theta^{*} (wrt proton beam) < 3.5 mRad");
    tex5b->SetTextSize(0.03);
    tex5b->Draw();

    TCanvas *c6a = new TCanvas("c6a");
    h3a_hcal->Draw();

    TLegend *leg6a = new TLegend(0.5,0.6,0.9,0.8);
    leg6a->SetBorderSize(0);leg6a->SetFillStyle(0);
    leg6a->SetHeader("Average over all events");
    leg6a->Draw();

    TCanvas *c6b = new TCanvas("c6b");
    h3b_hcal->Draw();

    TLegend *leg6b = new TLegend(0.5,0.6,0.9,0.8);
    leg6b->SetBorderSize(0);leg6b->SetFillStyle(0);
    leg6b->SetHeader("Average over events w/ neutron #theta^{*} < 3.5 mRad");
    leg6b->Draw();

    TCanvas *c6c = new TCanvas("c6c");
    h3c_hcal->Draw();

    TLegend *leg6c = new TLegend(0.5,0.6,0.9,0.8);
    leg6c->SetBorderSize(0);leg6c->SetFillStyle(0);
    leg6c->SetHeader("Average over events w/ neutron 6 < #theta^{*} < 8 mRad");
    leg6c->Draw();

    //Print plots to file
    c1->Print(Form("%s[",outputfile.c_str()));
    c1->Print(Form("%s",outputfile.c_str()));
    c2->Print(Form("%s",outputfile.c_str()));
    c2a->Print(Form("%s",outputfile.c_str()));
    c3->Print(Form("%s",outputfile.c_str()));
    c4->Print(Form("%s",outputfile.c_str()));
    c4a->Print(Form("%s",outputfile.c_str()));
    c5->Print(Form("%s",outputfile.c_str()));
    c5a->Print(Form("%s",outputfile.c_str()));
    c5b->Print(Form("%s",outputfile.c_str()));
    c6a->Print(Form("%s",outputfile.c_str()));
    c6b->Print(Form("%s",outputfile.c_str()));
    c6c->Print(Form("%s",outputfile.c_str()));
    c6c->Print(Form("%s]",outputfile.c_str()));

}
