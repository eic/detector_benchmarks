#include <fstream>
#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/RotationY.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

void trk_dis_analysis(const std::string& config_name)
{
    // Read our configuration
    std::ifstream  config_file{config_name};
    nlohmann::json config;
    config_file >> config;
  
    const std::string rec_file      = config["rec_file"];
    const std::string detector      = config["detector"];
    const std::string output_prefix = config["output_prefix"];
    const int         ebeam         = config["ebeam"];
    const int         pbeam         = config["pbeam"];
    const int         Q2_min        = config["Min_Q2"];
    
    fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
                "Running DIS tracking analysis...\n");
    fmt::print(" - Detector package: {}\n", detector);
    fmt::print(" - input file: {}\n", rec_file);
    fmt::print(" - output prefix for histograms: {}\n", output_prefix);
    fmt::print(" - ebeam: {}\n", ebeam);
    fmt::print(" - pbeam: {}\n", pbeam);
    fmt::print(" - Minimum Q2: {}\n", Q2_min);
    
    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Set output file for the histograms
    std::string output_name_hists = fmt::format("{}.root", output_prefix);
    cout << "Output file for histograms = " << output_name_hists << endl;
    TFile* ofile = new TFile(output_name_hists.c_str(), "RECREATE");

    //--------------------------------------------------------------------------------------------------------------------------------------------
 
    // Set up input file chain
    TChain *mychain = new TChain("events");
    mychain->Add(rec_file.c_str());

    //--------------------------------------------------------------------------------------------------------------------------------------------

    // Initialize reader
    TTreeReader tr(mychain);

    // Generated Particle Information
    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<float> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<float> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<float> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass"); //Not important here
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");

    // Reconstructed real-seeded tracks (charged particles)
    TTreeReaderArray<float> rec_px(tr, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> rec_py(tr, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> rec_pz(tr, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> rec_mass(tr, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<int> rec_type(tr, "ReconstructedChargedParticles.type"); //Type 0: successful eta/phi match to generated particle
    TTreeReaderArray<int> rec_pdg(tr, "ReconstructedChargedParticles.PDG"); //Uses PID lookup table information

    // Reconstructed truth-seeded tracks (charged particles)
    TTreeReaderArray<float> rec_ts_px(tr, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> rec_ts_py(tr, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> rec_ts_pz(tr, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TTreeReaderArray<float> rec_ts_mass(tr, "ReconstructedTruthSeededChargedParticles.mass");
    TTreeReaderArray<int> rec_ts_type(tr, "ReconstructedTruthSeededChargedParticles.type"); //Type 0: successful eta/phi match to generated particle
    TTreeReaderArray<int> rec_ts_pdg(tr, "ReconstructedTruthSeededChargedParticles.PDG"); //Uses PID lookup table information

    //-------------------------------------------------------------------------------------------------------------------------------------------- 
    // Define Histograms

    //Eta distribution of generated charged particles
    TH1 *h1a = new TH1D("h1a","Generated Charged Particles",100,-4,4);
    h1a->GetXaxis()->SetTitle("#eta_{gen.}");h1a->GetXaxis()->CenterTitle();
    h1a->SetLineColor(kTeal);h1a->SetLineWidth(2);

    TH1 *h1a1 = new TH1D("h1a1","Generated Charged Particles",100,-4,4);  //Minimum momentum cut of Pt > 200 MeV/c
    h1a1->GetXaxis()->SetTitle("#eta_{gen.}");h1a1->GetXaxis()->CenterTitle();
    h1a1->SetLineColor(kRed);h1a1->SetLineWidth(2);

    TH1 *h1a2 = new TH1D("h1a2","Generated Charged Particles",100,-4,4);  //Minimum momentum cut of Pt > 500 MeV/c
    h1a2->GetXaxis()->SetTitle("#eta_{gen.}");h1a2->GetXaxis()->CenterTitle();
    h1a2->SetLineColor(kBlack);h1a2->SetLineWidth(2);
    h1a2->SetFillColor(kBlack);h1a2->SetFillStyle(3244);

    //Eta distribution of reconstructed real-seeded tracks (charged particles)
    TH1 *h1b = new TH1D("h1b","Reconstructed Real-seeded tracks",100,-4,4);
    h1b->GetXaxis()->SetTitle("#eta_{rec.}");h1b->GetXaxis()->CenterTitle();
    h1b->SetLineColor(kGreen);h1b->SetLineWidth(2);

    TH1 *h1b1 = new TH1D("h1b1","Reconstructed Real-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 200 MeV/c
    h1b1->GetXaxis()->SetTitle("#eta_{rec.}");h1b1->GetXaxis()->CenterTitle();
    h1b1->SetLineColor(kBlue);h1b1->SetLineWidth(2);

    TH1 *h1b2 = new TH1D("h1b2","Reconstructed Real-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 500 MeV/c
    h1b2->GetXaxis()->SetTitle("#eta_{rec.}");h1b2->GetXaxis()->CenterTitle();
    h1b2->SetLineColor(kRed);h1b2->SetLineWidth(2);
    h1b2->SetMarkerColor(kRed);h1b2->SetMarkerStyle(kFullCrossX);

    //Eta distribution of reconstructed truth-seeded tracks (charged particles)
    TH1 *h1c = new TH1D("h1c","Reconstructed Truth-seeded tracks",100,-4,4);
    h1c->GetXaxis()->SetTitle("#eta_{rec.}");h1c->GetXaxis()->CenterTitle();
    h1c->SetLineColor(kRed);h1c->SetLineWidth(2);

    TH1 *h1c1 = new TH1D("h1c1","Reconstructed Truth-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 200 MeV/c
    h1c1->GetXaxis()->SetTitle("#eta_{rec.}");h1c1->GetXaxis()->CenterTitle();
    h1c1->SetLineColor(kOrange);h1c1->SetLineWidth(2);

    TH1 *h1c2 = new TH1D("h1c2","Reconstructed Truth-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 500 MeV/c
    h1c2->GetXaxis()->SetTitle("#eta_{rec.}");h1c2->GetXaxis()->CenterTitle();
    h1c2->SetLineColor(kMagenta);h1c2->SetLineWidth(2);
    h1c2->SetMarkerColor(kMagenta);h1c2->SetMarkerStyle(kFullCrossX);

    // Ratio histograms
    TH1 *h1rb1 = new TH1D("h1rb1","",100,-4,4); //Real-seeded tracks (Pt > 200 MeV/c cut)
    TH1 *h1rc1 = new TH1D("h1rc1","",100,-4,4); //Truth-seeded tracks (Pt > 200 MeV/c cut)
    TH1 *h1rb2 = new TH1D("h1rb2","",100,-4,4); //Real-seeded tracks (Pt > 500 MeV/c cut)
    TH1 *h1rc2 = new TH1D("h1rc2","",100,-4,4); //Truth-seeded tracks (Pt > 500 MeV/c cut)

    //Define additional variables
    TLorentzVector gen_vec;
    TVector3 gen_vertex;

    TLorentzVector rec_vec;
    TVector3 track_vec; //Reconstructed track momentum vector
    int counter(0);

    //Loop over events
    std::cout<<"Analyzing "<<mychain->GetEntries()<<" events!"<<std::endl;
    while (tr.Next()) {

	    if(counter%100==0) std::cout<<"Analyzing event "<<counter<<std::endl;
	    counter++;

        //Loop over generated particles
        for(size_t igen=0;igen<gen_status.GetSize();igen++){
	        
            auto charge = gen_charge[igen];
            auto status = gen_status[igen];

            //Require final-state, charged particle (no secondaries)
            if(status==1 && fabs(charge) > 0.01 ){
                gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
                
                //Fill eta histogram
                h1a->Fill(gen_vec.Eta());
                if( gen_vec.Pt()>0.2 ) h1a1->Fill(gen_vec.Eta());
                if( gen_vec.Pt()>0.5 ) h1a2->Fill(gen_vec.Eta());
            }
        } //End loop over generated particles
        
        //Loop over reconstructed real-seeded charged particles (copy of tracks with PID info)
        size_t rec_mult = rec_type.GetSize();

        for(size_t irec=0;irec<rec_mult;irec++){

            rec_vec.SetXYZM(rec_px[irec],rec_py[irec],rec_pz[irec],rec_mass[irec]);
            
            //Fill histograms
            h1b->Fill(rec_vec.Eta());
            if( rec_vec.Pt() > 0.2 ) h1b1->Fill(rec_vec.Eta());
            if( rec_vec.Pt() > 0.5 ) h1b2->Fill(rec_vec.Eta());

        } //End loop over reconstructed particles

        //Loop over reconstructed truth-seeded charged particles (copy of tracks with PID info)
        size_t rec_ts_mult = rec_ts_type.GetSize();

        for(size_t irec=0;irec<rec_ts_mult;irec++){

            rec_vec.SetXYZM(rec_ts_px[irec],rec_ts_py[irec],rec_ts_pz[irec],rec_ts_mass[irec]);
            
            //Fill histograms
            h1c->Fill(rec_vec.Eta());
            if( rec_vec.Pt() > 0.2 ) h1c1->Fill(rec_vec.Eta());
            if( rec_vec.Pt() > 0.5 ) h1c2->Fill(rec_vec.Eta());

        } //End loop over reconstructed particles

    } //End loop over events

    //Make ratio histograms
    h1rb1 = (TH1*) h1b1->Clone("h1rb1");
    h1rb1->Divide(h1a1);
    h1rb1->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 200 MeV/c");
    h1rb1->GetXaxis()->SetTitle("#eta");h1rb1->GetXaxis()->CenterTitle();
    h1rb1->GetYaxis()->SetTitle("Ratio");h1rb1->GetYaxis()->CenterTitle();

    h1rc1 = (TH1*) h1c1->Clone("h1rc1");
    h1rc1->Divide(h1a1);
    h1rc1->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 200 MeV/c");
    h1rc1->GetXaxis()->SetTitle("#eta");h1rc1->GetXaxis()->CenterTitle();
    h1rc1->GetYaxis()->SetTitle("Ratio");h1rc1->GetYaxis()->CenterTitle();

    h1rb2 = (TH1*) h1b2->Clone("h1rb2");
    h1rb2->Divide(h1a2);
    h1rb2->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 500 MeV/c");
    h1rb2->GetXaxis()->SetTitle("#eta");h1rb2->GetXaxis()->CenterTitle();
    h1rb2->GetYaxis()->SetTitle("Ratio");h1rb2->GetYaxis()->CenterTitle();

    h1rc2 = (TH1*) h1c2->Clone("h1rc2");
    h1rc2->Divide(h1a2);
    h1rc2->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 500 MeV/c");
    h1rc2->GetXaxis()->SetTitle("#eta");h1rc2->GetXaxis()->CenterTitle();
    h1rc2->GetYaxis()->SetTitle("Ratio");h1rc2->GetYaxis()->CenterTitle();

     //--------------------------------------------------------------------------------------------------------------------------------------------

    ofile->Write(); // Write histograms to file
    ofile->Close(); // Close output file

}