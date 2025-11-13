//-------------------------
//
// Simple analysis code to analyze EPIC simulation output 
//
//
// Author: Alex Jentsch
//
//
//
//------------------------


#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TString.h>
#include <TVector3.h>

using namespace std;

#include "detectorResolution.hpp"
#include "render.hpp"

void analyze_DVCS_eicrecon(TString fileList = "inputFileList.list", TString outputDir = "."){
	
	cout << "Input FileList: " << fileList << endl;
	string fileName;
	TFile * inputRootFile;
	TTree * rootTree;
	TString outputFileName = outputDir + "/output.root";
	cout << "Output file: " << outputFileName << endl;


	ifstream fileListStream;
	fileListStream.open(fileList);
	if(!fileListStream) { cout << "NO_LIST_FILE " << fileList << endl; return;}

	//--------------------------------------------------------------------------

	//histograms -- only a few for now
	
	//MC information
	TH1D* h_eta_MC = new TH1D("h_eta",";Pseudorapidity, #eta",100,-10.0,10.0);
	TH1D* h_px_MC = new TH1D("px_MC", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_MC = new TH1D("py_MC", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_MC = new TH1D("pt_MC", ";p_{t} [GeV/c]", 150, 0.0, 2.0);
	TH1D* h_pz_MC = new TH1D("pz_MC", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	TH1D* h_theta_MC = new TH1D("theta_MC", ";#theta [mrad]", 100, 0.0, 25.0);
	TH1D* h_phi_MC = new TH1D("phi_MC", ";#phi [rad]", 100, -3.2, 3.2);
	
	TH1D* h_pt_squared_MC = new TH1D("pt_squared_MC", "; Momentum Transfer, -t [GeV^{2}]", 100, 0.0, 2.0);
	
	//TRUTH information
	TH1D* h_eta_TRUTH = new TH1D("h_eta_TRUTH",";Pseudorapidity, #eta",100,-10.0,10.0);
	TH1D* h_px_TRUTH = new TH1D("px_TRUTH", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_TRUTH = new TH1D("py_TRUTH", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_TRUTH = new TH1D("pt_TRUTH", ";p_{t} [GeV/c]", 150, 0.0, 2.0);
	TH1D* h_pz_TRUTH = new TH1D("pz_TRUTH", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	TH1D* h_theta_TRUTH = new TH1D("theta_TRUTH", ";#theta [mrad]", 100, 0.0, 25.0);
	TH1D* h_phi_TRUTH = new TH1D("phi_TRUTH", ";#phi [rad]", 100, -3.2, 3.2);
	
	TH1D* h_pt_squared_TRUTH = new TH1D("pt_squared_TRUTH", "; Momentum Transfer, -t [GeV^{2}]", 100, 0.0, 2.0);
	
	//ACCEPTANCE ONLY information
	TH1D* h_eta_accep = new TH1D("h_eta_accep",";Pseudorapidity, #eta",100,-10.0,10.0);
	TH1D* h_px_accep = new TH1D("px_accep", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_accep = new TH1D("py_accep", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_accep = new TH1D("pt_accep", ";p_{t} [GeV/c]", 150, 0.0, 2.0);
	TH1D* h_pz_accep = new TH1D("pz_accep", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	TH1D* h_theta_accep = new TH1D("theta_accep", ";#theta [mrad]", 100, 0.0, 25.0);
	TH1D* h_phi_accep = new TH1D("phi_accep", ";#phi [rad]", 100, -3.2, 3.2);
	
	TH1D* h_pt_squared_accep = new TH1D("pt_squared_accep", "; Momentum Transfer, -t [GeV^{2}]", 100, 0.0, 2.0);
	
	//Roman pots
	TH1D* h_px_RomanPots = new TH1D("px_RomanPots", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_RomanPots = new TH1D("py_RomanPots", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_RomanPots = new TH1D("pt_RomanPots", ";p_{t} [GeV/c]", 150, 0.0, 2.0);
	TH1D* h_pz_RomanPots = new TH1D("pz_RomanPots", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	TH1D* h_pt_squared_RomanPots = new TH1D("pt_squared_RomanPots", "; Momentum Transfer, -t [GeV^{2}]", 100, 0.0, 2.0);
	TH2D* h_rp_GLOBAL_occupancy_map = new TH2D("Roman_pots_GLOBAL_occupancy_map", ";hit x [mm];hit y [mm]", 100, -1300, -900, 100, -70, 70);
	TH2D* h_rp_LOCAL_occupancy_map = new TH2D("Roman_pots_LOCAL_occupancy_map", ";hit x [mm];hit y [mm]", 100, -150, 150, 100, -80, 80);
	//100, -150, 150, 100, -70, 70);

	//OMD
    TH1D* h_px_OMD = new TH1D("px_OMD", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
    TH1D* h_py_OMD = new TH1D("py_OMD", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
    TH1D* h_pt_OMD = new TH1D("pt_OMD", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
    TH1D* h_pz_OMD = new TH1D("pz_OMD", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
    TH2D* h_omd_occupancy_map = new TH2D("OMD_occupancy_map", ";hit x [mm];hit y [mm]", 100, -150, 150, 100, -70, 70);	


	//B0 tracker hits
	TH2D* h_B0_occupancy_map_layer_0 = new TH2D("B0_occupancy_map_0", "B0_occupancy_map_0", 100, -400, 0, 100, -170, 170);
	TH2D* h_B0_occupancy_map_layer_1 = new TH2D("B0_occupancy_map_1", "B0_occupancy_map_1", 100, -400, 0, 100, -170, 170);
	TH2D* h_B0_occupancy_map_layer_2 = new TH2D("B0_occupancy_map_2", "B0_occupancy_map_2", 100, -400, 0, 100, -170, 170);
	TH2D* h_B0_occupancy_map_layer_3 = new TH2D("B0_occupancy_map_3", "B0_occupancy_map_3", 100, -400, 0, 100, -170, 170);
	TH1D* h_B0_hit_energy_deposit = new TH1D("B0_tracker_hit_energy_deposit", ";Deposited Energy [keV]", 100, 0.0, 500.0);
	
	//B0 EMCAL clusters
	TH2D* h_B0_emcal_occupancy_map = new TH2D("B0_emcal_occupancy_map", ";hit x [mm];hit y [mm]", 100, -400, 0, 100, -170, 170);
	TH1D* h_B0_emcal_cluster_energy = new TH1D("B0_emcal_cluster_energy", ";Cluster Energy [GeV]", 100, 0.0, 100.0);
	
	//Reconstructed tracks (for usage with B0 too!!)
	TH1D* h_px_reco_track = new TH1D("px_reco_track", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_py_reco_track = new TH1D("py_reco_track", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
	TH1D* h_pt_reco_track = new TH1D("pt_reco_track", ";p_{t} [GeV/c]", 150, 0.0, 2.0);
	TH1D* h_pz_reco_track = new TH1D("pz_reco_track", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
	TH1D* h_pt_squared_B0 = new TH1D("pt_squared_B0", "; Momentum Transfer, -t [GeV^{2}]", 100, 0.0, 2.0);	

	//ZDC EMCAL clusters
	TH2D* h_ZDC_emcal_occupancy_map = new TH2D("ZDC_emcal_occupancy_map", "ZDC_emcal_occupancy_map", 100, -1150, -1050, 100, -60, 60);
	TH1D* h_ZDC_emcal_cluster_energy = new TH1D("ZDC_emcal_cluster_energy", ";Cluster Energy [GeV]", 100, 0.0, 100.0);
	
	//B0 momentum resolution

	TH1D* h_b0_pt_resolution = new TH1D("b0_pt_resolution", ";#Delta p_{T} [GeV/c]", 100, -2.0, 2.0);
    TH2D* h_b0_pt_resolution_percent = new TH2D("b0_deltaPt_over_pt_vs_pt", ";P_{T, MC} [GeV/c]; #Delta p_{T}/p_{T, MC} [percent/100]", 100, 0.0, 2.5, 100, -1.0, 1.0);	
	TH2D* h_b0_p_resolution_percent = new TH2D("b0_deltaP_over_p_vs_p", ";Three-Momentum,  p_{MC} [GeV/c]; #Delta p/p_{MC} [percent/100]", 100, 0.0, 20.0, 100, -1.0, 1.0);	
    TH2D* h_b0_px_resolution_percent = new TH2D("b0_deltaPx_over_px_vs_px", ";  P_{x, MC} [GeV/c]; #Delta p_{x}/p_{x, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0);	
	TH2D* h_b0_py_resolution_percent = new TH2D("b0_deltaPy_over_py_vs_py", ";  p_{y, MC} [GeV/c]; #Delta p_{y}/p_{y, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0);	
	TH2D* h_b0_pz_resolution_percent = new TH2D("b0_deltaPz_over_pz_vs_pz", ";  p_{z, MC} [GeV/c]; #Delta p_{z}/p_{z, MC} [percent/100]", 100, 0.0, 320.0, 100, -1.0, 1.0);	
	
	TH1D* h_b0_extracted_pt_resolution;
	TH1D* h_b0_extracted_p_resolution;
	
	TH1D* h_extracted_px_resolution;
	TH1D* h_extracted_py_resolution;
	TH1D* h_extracted_pz_resolution;
	
	//RP momentum resolution
	

	TH1D* h_RP_pt_resolution = new TH1D("RP_pt_resolution", ";#Delta p_{T} [GeV/c]", 100, -2.0, 2.0);
    TH2D* h_RP_pt_resolution_percent = new TH2D("RP_deltaPt_over_pt_vs_pt", ";P_{T, MC} [GeV/c]; #Delta p_{T}/p_{T, MC} [percent/100]", 100, 0.0, 2.5, 100, -1.0, 1.0);	
	TH2D* h_RP_p_resolution_percent = new TH2D("RP_deltaP_over_p_vs_p", ";Three-Momentum,  p_{MC} [GeV/c]; #Delta p/p_{MC} [percent/100]", 100, 0.0, 320.0, 100, -1.0, 1.0);	
    TH2D* h_RP_px_resolution_percent = new TH2D("RP_deltaPx_over_px_vs_px", ";  P_{x, MC} [GeV/c]; #Delta p_{x}/p_{x, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0);	
	TH2D* h_RP_py_resolution_percent = new TH2D("RP_deltaPy_over_py_vs_py", ";  p_{y, MC} [GeV/c]; #Delta p_{y}/p_{y, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0);	
	TH2D* h_RP_pz_resolution_percent = new TH2D("RP_deltaPz_over_pz_vs_pz", ";  p_{z, MC} [GeV/c]; #Delta p_{z}/p_{z, MC} [percent/100]", 100, 0.0, 320.0, 100, -1.0, 1.0);	
	
	TH1D* h_RP_extracted_pt_resolution;
	TH1D* h_RP_extracted_p_resolution;
	
	TH1D* h_extracted_RP_px_resolution;
	TH1D* h_extracted_RP_py_resolution;
	TH1D* h_extracted_RP_pz_resolution;

	//forward ECAL information

	TH2D* h_forward_emcal_occupancy_map = new TH2D("forward_emcal_occupancy_map", "forward_emcal_occupancy_map", 100, -1000, 1000, 100, -1000, 1000);
	TH1D* h_num_forward_emcal_clusters = new TH1D("number_of_forward_EMCAL_clusters_per_event", "number_of_forward_EMCAL_clusters_per_event", 50, 0, 50);
	
	//barrel ECAL information

	//TH2D* h_forward_emcal_occupancy_map = new TH2D("forward_emcal_occupancy_map", "forward_emcal_occupancy_map", 100, -1000, 1000, 100, -1000, 1000);
	//TH1D* h_num_forward_emcal_clusters = new TH1D("number_of_forward_EMCAL_clusters_per_event", "number_of_forward_EMCAL_clusters_per_event", 50, 0, 50);
	

	int fileCounter = 0;
	int iEvent = 0;

	while(getline(fileListStream, fileName) && iEvent < 150000){

	    TString tmp = fileName;

	    cout << "Input file " << fileCounter << ": " << fileName << endl;

	    inputRootFile = new TFile(tmp);
	    if(inputRootFile->IsZombie()){ cout << "MISSING_ROOT_FILE"<< fileName << endl; continue;}
		
		fileCounter++;

		TTree * evtTree = (TTree*)inputRootFile->Get("events");

		int numEvents = evtTree->GetEntries();

    	TTreeReader tree_reader(evtTree);       // !the tree reader

		//MC particles
    
    	TTreeReaderArray<double> mc_px_array = {tree_reader, "MCParticles.momentum.x"};
    	TTreeReaderArray<double> mc_py_array = {tree_reader, "MCParticles.momentum.y"};
    	TTreeReaderArray<double> mc_pz_array = {tree_reader, "MCParticles.momentum.z"};
    	TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
    	TTreeReaderArray<int> mc_pdg_array = {tree_reader, "MCParticles.PDG"};
		TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};

		TTreeReaderArray<double> truth_px_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.x"};
        TTreeReaderArray<double> truth_py_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.y"};
        TTreeReaderArray<double> truth_pz_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.z"};
        TTreeReaderArray<double> truth_mass_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.mass"};
        TTreeReaderArray<int> truth_pdg_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.PDG"};
        TTreeReaderArray<int> truth_genStatus_array = {tree_reader, "MCParticlesHeadOnFrameNoBeamFX.generatorStatus"};	
		
	
		//Roman pots -- momentum vector
   	 	TTreeReaderArray<float> reco_RP_px = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    	TTreeReaderArray<float> reco_RP_py = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    	TTreeReaderArray<float> reco_RP_pz = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
	
		//Off-Momentum -- momentum vector
		TTreeReaderArray<float> reco_OMD_px = {tree_reader, "ForwardOffMRecParticles.momentum.x"};
        TTreeReaderArray<float> reco_OMD_py = {tree_reader, "ForwardOffMRecParticles.momentum.y"};
        TTreeReaderArray<float> reco_OMD_pz = {tree_reader, "ForwardOffMRecParticles.momentum.z"};
	
		//hit locations (for debugging)
   	 	TTreeReaderArray<float> global_hit_RP_x = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.x"};
    	TTreeReaderArray<float> global_hit_RP_y = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.y"};
    	TTreeReaderArray<float> global_hit_RP_z = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.z"};

		//hit locations (for debugging)
        TTreeReaderArray<float> global_hit_OMD_x = {tree_reader, "ForwardOffMRecParticles.referencePoint.x"};
        TTreeReaderArray<float> global_hit_OMD_y = {tree_reader, "ForwardOffMRecParticles.referencePoint.y"};
        TTreeReaderArray<float> global_hit_OMD_z = {tree_reader, "ForwardOffMRecParticles.referencePoint.z"};
		
		//b0 tracker hits
		TTreeReaderArray<float> b0_hits_x = {tree_reader, "B0TrackerRecHits.position.x"};
    	TTreeReaderArray<float> b0_hits_y = {tree_reader, "B0TrackerRecHits.position.y"};
    	TTreeReaderArray<float> b0_hits_z = {tree_reader, "B0TrackerRecHits.position.z"};
		TTreeReaderArray<float> b0_MC_momentum_x = {tree_reader, "B0TrackerHits.momentum.x"};
        TTreeReaderArray<float> b0_MC_momentum_y = {tree_reader, "B0TrackerHits.momentum.y"};
        TTreeReaderArray<float> b0_MC_momentum_z = {tree_reader, "B0TrackerHits.momentum.z"};
		TTreeReaderArray<float> b0_hits_eDep = {tree_reader, "B0TrackerRecHits.edep"}; //deposited energy per hit
	
	
		//b0 EMCAL
		TTreeReaderArray<float> b0_cluster_x = {tree_reader, "B0ECalClusters.position.x"};
    	TTreeReaderArray<float> b0_cluster_y = {tree_reader, "B0ECalClusters.position.y"};
    	TTreeReaderArray<float> b0_cluster_z = {tree_reader, "B0ECalClusters.position.z"};
		TTreeReaderArray<float>  b0_cluster_energy = {tree_reader, "B0ECalClusters.energy"}; //deposited energy in cluster
		
		//reco tracks (where b0 tracks live!!!)
		//TTreeReaderArray<float> reco_track_x = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
    	//TTreeReaderArray<float> reco_track_y = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
    	//TTreeReaderArray<float> reco_track_z = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
		//TTreeReaderArray<int> reco_track_PDG = {tree_reader, "ReconstructedChargedParticles.PDG"};
		
		TTreeReaderArray<float> reco_track_x = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    	TTreeReaderArray<float> reco_track_y = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    	TTreeReaderArray<float> reco_track_z = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
		TTreeReaderArray<int> reco_track_PDG = {tree_reader, "ReconstructedTruthSeededChargedParticles.PDG"};
		
		//forward EMCAL endcap
		TTreeReaderArray<float> pECAL_cluster_x = {tree_reader, "EcalEndcapPTruthClusters.position.x"};
        TTreeReaderArray<float> pECAL_cluster_y = {tree_reader, "EcalEndcapPTruthClusters.position.y"};
        TTreeReaderArray<float> pECAL_cluster_z = {tree_reader, "EcalEndcapPTruthClusters.position.z"};
        TTreeReaderArray<float>  pECAL_cluster_energy = {tree_reader, "EcalEndcapPTruthClusters.energy"};

		//ZDC EMCAL
		//TTreeReaderArray<float> zdc_ecal_cluster_x = {tree_reader, "ZDCEcalClusters.position.x"};
    	//TTreeReaderArray<float> zdc_ecal_cluster_y = {tree_reader, "ZDCEcalClusters.position.y"};
    	//TTreeReaderArray<float> zdc_ecal_cluster_z = {tree_reader, "ZDCEcalClusters.position.z"};
		//TTreeReaderArray<float>  zdc_ecal_cluster_energy = {tree_reader, "ZDCEcalClusters.energy"}; //deposited energy in cluster
		

		cout << "file has " << evtTree->GetEntries() <<  " events..." << endl;

		tree_reader.SetEntriesRange(0, evtTree->GetEntries());
		while (tree_reader.Next()) {

			//cout << "Reading event: " << iEvent << endl;

	    	//MCParticles
	        
			TVector3 mctrk;
			TVector3 rptrk;
			TVector3 final_mc_track;
		
			double protonMomentum = 0.0;
	
			double maxPt=-99.;			
			
			//loop over TRUTH particles where the afterburner effects + crossing angle are compeltely removed
			for(int imc = 0; imc < truth_px_array.GetSize(); imc++){
				mctrk.SetXYZ(truth_px_array[imc], truth_py_array[imc], truth_pz_array[imc]);
				final_mc_track.SetXYZ(-999, -999, -999);

				if(truth_pdg_array[imc] == 2212 && truth_genStatus_array[imc] == 1){ //only checking final-state protons here
				
					if(mctrk.Mag() < 0.0){ continue; }

					h_eta_TRUTH->Fill(mctrk.Eta()); 
			    	h_px_TRUTH->Fill(mctrk.Px()); 
			    	h_py_TRUTH->Fill(mctrk.Py()); 
			    	h_pt_TRUTH->Fill(mctrk.Perp()); 
			    	h_pz_TRUTH->Fill(mctrk.Pz()); 
			    	h_theta_TRUTH->Fill(1000*mctrk.Theta()); 
			    	h_phi_TRUTH->Fill(mctrk.Phi()); 
				
					h_pt_squared_TRUTH->Fill(mctrk.Perp()*mctrk.Perp());
					
					//if(reco_RP_px.GetSize() > 0){h_pt_RomanPots->Fill(mctrk.Perp());} //for acceptance only plot
					
					//final_mc_track = mctrk;
					//break;
				}
			}
						
			//MC Particles loop
	    	for(int imc=0;imc<mc_px_array.GetSize();imc++){
	    		mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);	
								
				//He4 -- 1000020040
	    		if(mc_pdg_array[imc] == 2212 && mc_genStatus_array[imc] == 1 ){ //only checking for protons here -- change as desired
					mctrk.RotateY(0.025);
				
					protonMomentum = mctrk.Mag();
	
					h_eta_MC->Fill(mctrk.Eta());
					
					h_px_MC->Fill(mctrk.Px());
					h_py_MC->Fill(mctrk.Py());
					h_pt_MC->Fill(mctrk.Perp());
					h_pz_MC->Fill(mctrk.Pz());
					h_theta_MC->Fill(mctrk.Theta()*1000);
					h_phi_MC->Fill(mctrk.Phi());
					
					h_pt_squared_MC->Fill(mctrk.Perp()*mctrk.Perp());
					
					//if(reco_RP_px.GetSize() > 0){h_pt_RomanPots->Fill(mctrk.Perp());}//for acceptance only plot
					final_mc_track = mctrk;					

				}
				
	    	}			
	    	
			//roman pots reco tracks
			for(int iRPPart = 0; iRPPart < reco_RP_px.GetSize(); iRPPart++){
    		
				TVector3 prec_romanpots(reco_RP_px[iRPPart], reco_RP_py[iRPPart], reco_RP_pz[iRPPart]);	
    		
				h_px_RomanPots->Fill(prec_romanpots.Px());
				h_py_RomanPots->Fill(prec_romanpots.Py());
				h_pt_RomanPots->Fill(prec_romanpots.Perp());
				//h_pt_RomanPots->Fill(final_mc_track.Perp());
				h_pz_RomanPots->Fill(prec_romanpots.Pz());
				h_pt_squared_RomanPots->Fill(prec_romanpots.Perp()*prec_romanpots.Perp());
				
				
				h_eta_accep->Fill(final_mc_track.Eta());
				h_px_accep->Fill(final_mc_track.Px());
				h_py_accep->Fill(final_mc_track.Py());
				h_pt_accep->Fill(final_mc_track.Perp());
				h_pz_accep->Fill(final_mc_track.Pz());
				h_theta_accep->Fill(final_mc_track.Theta()*1000);
				h_phi_accep->Fill(final_mc_track.Phi());
				h_pt_squared_accep->Fill(final_mc_track.Perp()*final_mc_track.Perp());
				
				double delPt = prec_romanpots.Perp() - final_mc_track.Perp();
				double delP = prec_romanpots.Mag() - final_mc_track.Mag();

				double delPt_percent = delPt/final_mc_track.Perp();
				double delP_percent = delP/final_mc_track.Mag();

				h_RP_pt_resolution->Fill(delPt_percent);
				h_RP_pt_resolution_percent->Fill(final_mc_track.Perp(), delPt_percent);
				h_RP_p_resolution_percent->Fill(final_mc_track.Mag(), delP_percent);
				
				double delPx = prec_romanpots.Px() - final_mc_track.Px();
				double delPy = prec_romanpots.Py() - final_mc_track.Py();
				double delPz = prec_romanpots.Pz() - final_mc_track.Pz();

				double delPx_percent = delPx/final_mc_track.Px();
				double delPy_percent = delPy/final_mc_track.Py();
				double delPz_percent = delPz/final_mc_track.Pz();
				
				h_RP_px_resolution_percent->Fill(final_mc_track.Px(), delPx_percent);
				h_RP_py_resolution_percent->Fill(final_mc_track.Py(), delPy_percent);
				h_RP_pz_resolution_percent->Fill(final_mc_track.Pz(), delPz_percent);
			
				if(global_hit_RP_z[iRPPart] < 32555.0){ //only the first sensor plane
	
					h_rp_GLOBAL_occupancy_map->Fill(global_hit_RP_x[iRPPart], global_hit_RP_y[iRPPart]);
					h_rp_LOCAL_occupancy_map->Fill(global_hit_RP_x[iRPPart], global_hit_RP_y[iRPPart]);
				}
			}
			
			//OMD reco tracks
            for(int iOMDPart = 0; iOMDPart < reco_OMD_px.GetSize(); iOMDPart++){

                TVector3 prec_omd(reco_OMD_px[iOMDPart], reco_OMD_py[iOMDPart], reco_OMD_pz[iOMDPart]);

                h_px_OMD->Fill(prec_omd.Px());
                h_py_OMD->Fill(prec_omd.Py());
                h_pt_OMD->Fill(prec_omd.Perp());
                h_pz_OMD->Fill(prec_omd.Pz());

                h_omd_occupancy_map->Fill(global_hit_OMD_x[iOMDPart], global_hit_OMD_y[iOMDPart]);
            }
		
		
			
			double hit_x = -9999.;
			double hit_y = -9999.;
			double hit_z = -9999.;
			double hit_deposited_energy = -9999.;
			
			for(int b0cluster = 0; b0cluster < b0_cluster_x.GetSize(); b0cluster++){
			
				hit_x = b0_cluster_x[b0cluster];
				hit_y = b0_cluster_y[b0cluster];
				hit_z = b0_cluster_z[b0cluster];
				hit_deposited_energy = b0_cluster_energy[b0cluster]*1.246; //poor man's calibration constant, for now
							
				h_B0_emcal_occupancy_map->Fill(hit_x, hit_y);
				h_B0_emcal_cluster_energy->Fill(hit_deposited_energy);
				
			}


			//b0 tracker hits -- for debugging or external tracking
			for(int b0hit = 0; b0hit < b0_hits_x.GetSize(); b0hit++){
    		
				hit_x = b0_hits_x[b0hit];
				hit_y = b0_hits_y[b0hit];
				hit_z = b0_hits_z[b0hit];
				hit_deposited_energy = b0_hits_eDep[b0hit]*1e6; //convert GeV --> keV
				
				h_B0_hit_energy_deposit->Fill(hit_deposited_energy);
				
				//if(hit_deposited_energy < 10.0){ continue; } //threshold value -- 10 keV, arbitrary for now
								
				//ACLGAD layout
    			if(hit_z > 5700 && hit_z < 5990){ h_B0_occupancy_map_layer_0->Fill(hit_x, hit_y); }
				if(hit_z > 6100 && hit_z < 6200){ h_B0_occupancy_map_layer_1->Fill(hit_x, hit_y); }
				if(hit_z > 6400 && hit_z < 6500){ h_B0_occupancy_map_layer_2->Fill(hit_x, hit_y); }
				if(hit_z > 6700 && hit_z < 6750){ h_B0_occupancy_map_layer_3->Fill(hit_x, hit_y); }
				
			}
		
			
			//reconstructed tracks with ACTS -- used for B0
			//cout << "Event has " << reco_track_x.GetSize() << " reco tracks from ACTS..." << endl;
			
			bool hasB0Track = false;
            bool hasHit1 = false;
            bool hasHit2 = false;
            bool hasHit3 = false;
            bool hasHit4 = false;

			int b0Hit_layer_1 = 0;
			

			for(int iRecoTrk = 0; iRecoTrk < reco_track_x.GetSize(); iRecoTrk++){
    		
				TVector3 prec_reco_tracks(reco_track_x[iRecoTrk], reco_track_y[iRecoTrk], reco_track_z[iRecoTrk]);
				
				prec_reco_tracks.RotateY(0.025); //remove crossing angle

				//looking for B0 protons, only -- eta cut is cheap way to do it without associations
				if(reco_track_PDG[iRecoTrk] != 0 || prec_reco_tracks.Eta() < 4.0){ continue; } 
				
				h_px_reco_track->Fill(prec_reco_tracks.Px());
				h_py_reco_track->Fill(prec_reco_tracks.Py());
				h_pt_reco_track->Fill(prec_reco_tracks.Perp());
				h_pz_reco_track->Fill(prec_reco_tracks.Pz());
				h_pt_squared_B0->Fill(prec_reco_tracks.Perp()*prec_reco_tracks.Perp());
			

				double delPt = prec_reco_tracks.Perp() - final_mc_track.Perp();
				double delP = prec_reco_tracks.Mag() - final_mc_track.Mag();

				double delPt_percent = delPt/final_mc_track.Perp();
				double delP_percent = delP/final_mc_track.Mag();

				h_b0_pt_resolution->Fill(delPt_percent);
				h_b0_pt_resolution_percent->Fill(final_mc_track.Perp(), delPt_percent);
				h_b0_p_resolution_percent->Fill(final_mc_track.Mag(), delP_percent);
				
				double delPx = prec_reco_tracks.Px() - final_mc_track.Px();
				double delPy = prec_reco_tracks.Py() - final_mc_track.Py();
				double delPz = prec_reco_tracks.Pz() - final_mc_track.Pz();

				double delPx_percent = delPx/final_mc_track.Px();
				double delPy_percent = delPy/final_mc_track.Py();
				double delPz_percent = delPz/final_mc_track.Pz();
				
				h_b0_px_resolution_percent->Fill(final_mc_track.Px(), delPx_percent);
				h_b0_py_resolution_percent->Fill(final_mc_track.Py(), delPy_percent);
				h_b0_pz_resolution_percent->Fill(final_mc_track.Pz(), delPz_percent);
			}
				
		
				
			h_num_forward_emcal_clusters->Fill(pECAL_cluster_x.GetSize());

			for(int iForwardECALCluster = 0; iForwardECALCluster < pECAL_cluster_x.GetSize(); iForwardECALCluster++){

				h_forward_emcal_occupancy_map->Fill(pECAL_cluster_x[iForwardECALCluster], pECAL_cluster_y[iForwardECALCluster]);


			}

			iEvent++;
			
		}// event loop
		
		inputRootFile->Close();
		
	}// input file loop
	
	
	//Calculate detector resolutions here (comment out if not needed)
	
	
    h_b0_extracted_pt_resolution = extractResolution("b0_extracted_pt_resolution", h_b0_pt_resolution_percent, false);
    h_b0_extracted_pt_resolution->GetXaxis()->SetTitle("p_{T, MC} [GeV/c]");
    h_b0_extracted_pt_resolution->GetYaxis()->SetTitle("#Delta p_{T}/p_{T, MC} [percent/100]");

    h_b0_extracted_p_resolution = extractResolution("b0_extracted_p_resolution", h_b0_p_resolution_percent, false);
    h_b0_extracted_p_resolution->GetXaxis()->SetTitle("p_{MC} [GeV/c]");
    h_b0_extracted_p_resolution->GetYaxis()->SetTitle("#Delta p/p_{MC} [percent/100]");

	h_extracted_px_resolution = extractResolution("track_extracted_px_resolution", h_b0_px_resolution_percent, false);
	h_extracted_py_resolution = extractResolution("track_extracted_py_resolution", h_b0_py_resolution_percent, false);
	h_extracted_pz_resolution = extractResolution("track_extracted_pz_resolution", h_b0_pz_resolution_percent, false);
	
	h_extracted_px_resolution->GetXaxis()->SetTitle("p_{x, MC} [GeV/c]");
	h_extracted_px_resolution->GetYaxis()->SetTitle("#Delta p_{x}/p_{x, MC} [percent/100]");
	
	h_extracted_py_resolution->GetXaxis()->SetTitle("p_{y, MC} [GeV/c]");
	h_extracted_py_resolution->GetYaxis()->SetTitle("#Delta p_{y}/p_{y, MC} [percent/100]");
	
	h_extracted_pz_resolution->GetXaxis()->SetTitle("p_{z, MC} [GeV/c]");
	h_extracted_pz_resolution->GetYaxis()->SetTitle("#Delta p_{z}/p_{z, MC} [percent/100]");
	
    h_RP_extracted_pt_resolution = extractResolution("RP_extracted_pt_resolution", h_RP_pt_resolution_percent, false);
    h_RP_extracted_pt_resolution->GetXaxis()->SetTitle("p_{T, MC} [GeV/c]");
    h_RP_extracted_pt_resolution->GetYaxis()->SetTitle("#Delta p_{T}/p_{T, MC} [percent/100]");

    h_RP_extracted_p_resolution = extractResolution("RP_extracted_p_resolution", h_RP_p_resolution_percent, false);
    h_RP_extracted_p_resolution->GetXaxis()->SetTitle("p_{MC} [GeV/c]");
    h_RP_extracted_p_resolution->GetYaxis()->SetTitle("#Delta p/p_{MC} [percent/100]");

	h_extracted_RP_px_resolution = extractResolution("track_extracted_RP_px_resolution", h_RP_px_resolution_percent, false);
	h_extracted_RP_py_resolution = extractResolution("track_extracted_RP_py_resolution", h_RP_py_resolution_percent, false);
	h_extracted_RP_pz_resolution = extractResolution("track_extracted_RP_pz_resolution", h_RP_pz_resolution_percent, false);
	
	h_extracted_RP_px_resolution->GetXaxis()->SetTitle("p_{x, MC} [GeV/c]");
	h_extracted_RP_px_resolution->GetYaxis()->SetTitle("#Delta p_{x}/p_{x, MC} [percent/100]");
	
	h_extracted_RP_py_resolution->GetXaxis()->SetTitle("p_{y, MC} [GeV/c]");
	h_extracted_RP_py_resolution->GetYaxis()->SetTitle("#Delta p_{y}/p_{y, MC} [percent/100]");
	
	h_extracted_RP_pz_resolution->GetXaxis()->SetTitle("p_{z, MC} [GeV/c]");
	h_extracted_RP_pz_resolution->GetYaxis()->SetTitle("#Delta p_{z}/p_{z, MC} [percent/100]");

	//for quick checks - just print plots here

	/*

	TLegend * leg2 = new TLegend(0.55, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_pt_squared_MC, "Burned MC, X'ing angle removed", "p");
    leg2->AddEntry(h_pt_squared_TRUTH, "Generator Truth", "p");
	leg2->AddEntry(h_pt_squared_RomanPots, "Roman Pots FULL reconstruction", "p");
    leg2->AddEntry(h_pt_squared_B0, "B0 FULL reconstruction", "p");


    TCanvas * cannyCompare2 = new TCanvas("can2", "can2", 600, 500);

	//cannyCompare2->Divide(2,1);

	cannyCompare2->cd(1)->SetLogy();

    h_pt_squared_MC->SetLineColor(kBlack);
    h_pt_squared_MC->SetMarkerColor(kBlack);
    h_pt_squared_MC->SetMarkerStyle(20);
    h_pt_squared_MC->SetStats(0);
    h_pt_squared_MC->GetXaxis()->SetTitle("Momentum Transfer, -t [GeV^{2}]");

    h_pt_squared_TRUTH->SetLineColor(kRed);
    h_pt_squared_TRUTH->SetMarkerColor(kRed);
    h_pt_squared_TRUTH->SetMarkerStyle(33);
    h_pt_squared_TRUTH->SetStats(0);
    h_pt_squared_TRUTH->GetXaxis()->SetTitle("Momentum Transfer, -t [GeV^{2}]");

    h_pt_squared_RomanPots->SetLineColor(kGreen + 3);
    h_pt_squared_RomanPots->SetMarkerColor(kGreen + 3);
    h_pt_squared_RomanPots->SetMarkerStyle(23);
    h_pt_squared_RomanPots->SetStats(0);
    h_pt_squared_RomanPots->GetXaxis()->SetTitle("Momentum Transfer, -t [GeV^{2}]");

	h_pt_squared_B0->SetLineColor(kBlue + 3);
    h_pt_squared_B0->SetMarkerColor(kBlue + 3);
    h_pt_squared_B0->SetMarkerStyle(23);
    h_pt_squared_B0->SetStats(0);
    h_pt_squared_B0->GetXaxis()->SetTitle("Momentum Transfer, -t [GeV^{2}]");

    h_pt_squared_MC->Draw("EP");
    h_pt_squared_TRUTH->Draw("SAME EP");
    h_pt_squared_RomanPots->Draw("SAME EP");
	h_pt_squared_B0->Draw("SAME EP");

	leg2->Draw("SAME");

	TBox * RP_outline = new TBox(-128.0, -48.0, 128.0, 48.0);
	RP_outline->SetLineColor(kBlack);
	RP_outline->SetLineWidth(2);
	RP_outline->SetFillStyle(0);

	
	cannyCompare2->cd(2)->SetLogy();

	h_pz_MC->SetLineColor(kBlack);
    h_pz_MC->SetMarkerColor(kBlack);
    h_pz_MC->SetMarkerStyle(20);
    h_pz_MC->SetStats(0);
    h_pz_MC->GetXaxis()->SetTitle("Longitudinal Momentum, p_{z} [GeV/c]");

    h_pz_RomanPots->SetLineColor(kRed);
    h_pz_RomanPots->SetMarkerColor(kRed);
    h_pz_RomanPots->SetMarkerStyle(33);
    h_pz_RomanPots->SetStats(0);
    h_pz_RomanPots->GetXaxis()->SetTitle("Longitudinal Momentum, p_{z} [GeV/c]");

    h_pz_reco_track->SetLineColor(kGreen + 3);
    h_pz_reco_track->SetMarkerColor(kGreen + 3);
    h_pz_reco_track->SetMarkerStyle(23);
    h_pz_reco_track->SetStats(0);
    h_pz_reco_track->GetXaxis()->SetTitle("Longitudinal Momentum, p_{z} [GeV/c]");

    h_pz_MC->Draw("EP");
    h_pz_RomanPots->Draw("SAME EP");
    h_pz_reco_track->Draw("SAME EP");
	

	TString outputPtPlots = "pt_roman_pots_B0_comparisons_";
	outputPtPlots = outputPtPlots + outputName + date + run + fileType_PDF;

	cannyCompare2->SaveAs(outputPtPlots);

	//TCanvas * RP_occupancy_canvas = new TCanvas("canvas_rp", "canvas_rp", 500, 500);
	
	//h_rp_occupancy_map->Draw("COLZ");
	//RP_outline->Draw("SAME");

	*/

	SaveHistogramsFromList(*gDirectory->GetList(), outputDir + "/far_forward_dvcs_");

	TFile * outputFile = new TFile(outputFileName, "RECREATE");

	h_eta_MC->Write();
	h_px_MC->Write();
	h_py_MC->Write();
	h_pt_MC->Write();
	h_pz_MC->Write();
	h_theta_MC->Write();
	h_phi_MC->Write();
	h_pt_squared_MC->Write();
	
	h_eta_TRUTH->Write(); 
    h_px_TRUTH->Write(); 
    h_py_TRUTH->Write(); 
    h_pt_TRUTH->Write(); 
    h_pz_TRUTH->Write(); 
    h_theta_TRUTH->Write(); 
    h_phi_TRUTH->Write(); 
	h_pt_squared_TRUTH->Write();
	
	h_eta_accep->Write();
	h_px_accep->Write();
	h_py_accep->Write();
	h_pt_accep->Write();
	h_pz_accep->Write();
	h_theta_accep->Write();
	h_phi_accep->Write();
	h_pt_squared_accep->Write();
	
	h_px_RomanPots->Write();
	h_py_RomanPots->Write();
	h_pt_RomanPots->Write();
	h_pz_RomanPots->Write();
	h_pt_squared_RomanPots->Write();
	h_rp_GLOBAL_occupancy_map->Write();
    h_rp_LOCAL_occupancy_map->Write();

	h_px_OMD->Write();
    h_py_OMD->Write();
    h_pt_OMD->Write();
    h_pz_OMD->Write();
    h_omd_occupancy_map->Write(); 
	
	h_B0_occupancy_map_layer_0->Write();
	h_B0_occupancy_map_layer_1->Write();
	h_B0_occupancy_map_layer_2->Write();
	h_B0_occupancy_map_layer_3->Write();
	h_B0_hit_energy_deposit->Write();
	
	h_B0_emcal_occupancy_map->Write();
	h_B0_emcal_cluster_energy->Write();

	h_px_reco_track->Write();
	h_py_reco_track->Write();
	h_pt_reco_track->Write();
	h_pz_reco_track->Write();
	h_pt_squared_B0->Write();

	h_b0_pt_resolution->Write();
	h_b0_pt_resolution_percent->Write();
	h_b0_extracted_pt_resolution->Write();
	h_b0_p_resolution_percent->Write();
	h_b0_extracted_p_resolution->Write();
	
	h_extracted_px_resolution->Write();
	h_extracted_py_resolution->Write();
	h_extracted_pz_resolution->Write();
	
	h_RP_pt_resolution->Write();
	h_RP_pt_resolution_percent->Write();
	h_RP_extracted_pt_resolution->Write();
	h_RP_p_resolution_percent->Write();
	h_RP_extracted_p_resolution->Write();
	
	h_extracted_RP_px_resolution->Write();
	h_extracted_RP_py_resolution->Write();
	h_extracted_RP_pz_resolution->Write();

	h_ZDC_emcal_occupancy_map->Write();
	h_ZDC_emcal_cluster_energy->Write();

	h_forward_emcal_occupancy_map->Write();
	h_num_forward_emcal_clusters->Write();

	outputFile->Close();

	

    return;

}

