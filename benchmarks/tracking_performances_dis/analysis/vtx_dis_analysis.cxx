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

void vtx_dis_analysis(const std::string& config_name)
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

    // TTreeReader
    TTreeReader tree_reader(mychain);
  
    // Reco Vertex
    TTreeReaderArray<int> recoVtxType = {tree_reader, "CentralTrackVertices.type"};
    TTreeReaderArray<float> recoVtxX = {tree_reader, "CentralTrackVertices.position.x"};
    TTreeReaderArray<float> recoVtxY = {tree_reader, "CentralTrackVertices.position.y"};
    TTreeReaderArray<float> recoVtxZ = {tree_reader, "CentralTrackVertices.position.z"};
  
    TTreeReaderArray<unsigned int> assoPartBegin = {tree_reader, "CentralTrackVertices.associatedParticles_begin"};
    TTreeReaderArray<unsigned int> assoPartEnd = {tree_reader, "CentralTrackVertices.associatedParticles_end"};
    TTreeReaderArray<int> assoPartIndex = {tree_reader, "_CentralTrackVertices_associatedParticles.index"};

    // MC
    TTreeReaderArray<int> mcGenStat = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int> mcPDG = {tree_reader, "MCParticles.PDG"};
    TTreeReaderArray<float> mcCharge = {tree_reader, "MCParticles.charge"};
    TTreeReaderArray<float> mcMomX = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<float> mcMomY = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<float> mcMomZ = {tree_reader, "MCParticles.momentum.z"};

    TTreeReaderArray<double> mcVtxX = {tree_reader, "MCParticles.vertex.x"};
    TTreeReaderArray<double> mcVtxY = {tree_reader, "MCParticles.vertex.y"};
    TTreeReaderArray<double> mcVtxZ = {tree_reader, "MCParticles.vertex.z"};

    TTreeReaderArray<unsigned int> mcParentBegin = {tree_reader, "MCParticles.parents_begin"};
    TTreeReaderArray<unsigned int> mcParentEnd = {tree_reader, "MCParticles.parents_end"};
    TTreeReaderArray<int> mcParentIndex = {tree_reader, "_MCParticles_parents.index"};

    // Reco
    TTreeReaderArray<int> recoType = {tree_reader, "ReconstructedChargedParticles.type"};


    //-------------------------------------------------------------------------------------------------------------------------------------------- 
    // Define Histograms
    
    // Gen
    TH1D *numGenTracksHist = new TH1D("numGenTracks","",31,-0.5,30.5);
    
    TH2D *genVtxYvsXHist = new TH2D("genVtxYvsXHist","",200,-1.,1.,200,-1,1);
    genVtxYvsXHist->Draw("COLZ");
    genVtxYvsXHist->SetTitle("MC Vertex: v_{x} versus v_{y}");
    genVtxYvsXHist->GetXaxis()->SetTitle("x-coordinate (in mm)");
    genVtxYvsXHist->GetYaxis()->SetTitle("y-coordinate (in mm)");
    genVtxYvsXHist->GetXaxis()->CenterTitle(); genVtxYvsXHist->GetYaxis()->CenterTitle(); 
    
    TH2D *genVtxRvsZHist = new TH2D("genVtxRvsZHist","",200,-100.,100.,100,0,0.5);
    genVtxRvsZHist->Draw("COLZ");
    genVtxRvsZHist->SetTitle("MC Vertex: v_{r} versus v_{z}");
    genVtxRvsZHist->GetXaxis()->SetTitle("z-coordinate (in mm)");
    genVtxRvsZHist->GetYaxis()->SetTitle("#sqrt{x^{2} + y^{2}} (in mm)");
    genVtxRvsZHist->GetXaxis()->CenterTitle(); genVtxRvsZHist->GetYaxis()->CenterTitle(); 
  
    // Reco
    TH1D *numRecoTracksHist = new TH1D("numRecoTracks","",31,-0.5,30.5);
    
    TH1D *recoVtxEffHist = new TH1D("recoVtxEff","",4,-0.5,3.5);
    recoVtxEffHist->Draw("P");
    recoVtxEffHist->SetLineColor(kRed); recoVtxEffHist->SetMarkerColor(kRed);
    recoVtxEffHist->SetMarkerStyle(8); recoVtxEffHist->SetMarkerSize(1.2);
    recoVtxEffHist->SetTitle("Vertexing Efficiency");
    recoVtxEffHist->GetXaxis()->SetTitle("Vertex Count"); recoVtxEffHist->GetYaxis()->SetTitle("nEvents/total_events %");
    recoVtxEffHist->GetXaxis()->CenterTitle(); recoVtxEffHist->GetYaxis()->CenterTitle(); 

    TH2D *recoVtxYvsXHist = new TH2D("recoVtxYvsX","",200,-1.,1.,200,-1,1);
    recoVtxYvsXHist->Draw("COLZ");
    recoVtxYvsXHist->SetTitle("Reconstructed Vertex: v_{x} versus v_{y}");
    recoVtxYvsXHist->GetXaxis()->SetTitle("x-coordinate (in mm)"); recoVtxYvsXHist->GetYaxis()->SetTitle("y-coordinate (in mm)");
    recoVtxYvsXHist->GetXaxis()->CenterTitle(); recoVtxYvsXHist->GetYaxis()->CenterTitle();
    
    TH2D *recoVtxRvsZHist = new TH2D("recoVtxRvsZ","",200,-100.,100.,100,0.,0.8);
    recoVtxRvsZHist->Draw("COLZ");
    recoVtxRvsZHist->SetTitle("Reconstructed Vertex: v_{r} versus v_{z}");
    recoVtxRvsZHist->GetXaxis()->SetTitle("z-coordinate (in mm)");
    recoVtxRvsZHist->GetYaxis()->SetTitle("#sqrt{x^{2} + y^{2}} (in mm)");
    recoVtxRvsZHist->GetXaxis()->CenterTitle(); recoVtxRvsZHist->GetYaxis()->CenterTitle();

    //Reco Vs MC
    TH2D *recoVsMCTracksHist = new TH2D("recoVsMCTracks","",31,-0.5,30.5,31,-0.5,30.5);
    recoVsMCTracksHist->Draw("COLZ");
    recoVsMCTracksHist->SetTitle("Number of particles associated with vertex");
    recoVsMCTracksHist->GetXaxis()->SetTitle("N_{MC}"); recoVsMCTracksHist->GetYaxis()->SetTitle("N_{RC}");
    recoVsMCTracksHist->GetXaxis()->CenterTitle(); recoVsMCTracksHist->GetYaxis()->CenterTitle(); 
    
    // Resolution
    TH2D *vtxResXvsGenTrkHist = new TH2D("vtxResXvsGenTrk","",31,-0.5,30.5,200,-1,1);
    vtxResXvsGenTrkHist->Draw("COLZ");
    vtxResXvsGenTrkHist->SetTitle("Vertex Resolution X vs MC Tracks");
    vtxResXvsGenTrkHist->GetXaxis()->SetTitle("N_{MC}"); vtxResXvsGenTrkHist->GetYaxis()->SetTitle("recVtx_x - mcVtx_x (in mm)");
    vtxResXvsGenTrkHist->GetXaxis()->CenterTitle(); vtxResXvsGenTrkHist->GetYaxis()->CenterTitle();
    
    TH2D *vtxResYvsGenTrkHist = new TH2D("vtxResYvsGenTrk","",31,-0.5,30.5,200,-1,1);
    vtxResYvsGenTrkHist->Draw("COLZ");
    vtxResYvsGenTrkHist->SetTitle("Vertex Resolution Y vs MC Tracks");
    vtxResYvsGenTrkHist->GetXaxis()->SetTitle("N_{MC}");
    vtxResYvsGenTrkHist->GetYaxis()->SetTitle("recVtx_y - mcVtx_y (in mm)");
    vtxResYvsGenTrkHist->GetXaxis()->CenterTitle(); vtxResYvsGenTrkHist->GetYaxis()->CenterTitle();
    
    TH2D *vtxResZvsGenTrkHist = new TH2D("vtxResZvsGenTrk","",31,-0.5,30.5,200,-1,1);
    vtxResZvsGenTrkHist->Draw("COLZ");
    vtxResZvsGenTrkHist->SetTitle("Vertex Resolution Z vs MC Tracks");
    vtxResZvsGenTrkHist->GetXaxis()->SetTitle("N_{MC}");
    vtxResZvsGenTrkHist->GetYaxis()->SetTitle("recVtx_z - mcVtx_z (in mm)");
    vtxResZvsGenTrkHist->GetXaxis()->CenterTitle(); vtxResZvsGenTrkHist->GetYaxis()->CenterTitle();
  
    TH2D *vtxResXvsRecoTrkHist = new TH2D("vtxResXvsRecoTrk","",31,-0.5,30.5,200,-1,1);
    vtxResXvsRecoTrkHist->Draw("COLZ");
    vtxResXvsRecoTrkHist->SetTitle("Vertex Resolution X vs Reconstructed Tracks");
    vtxResXvsRecoTrkHist->GetXaxis()->SetTitle("N_{RC}");
    vtxResXvsRecoTrkHist->GetYaxis()->SetTitle("recVtx_x - mcVtx_x (in mm)");
    vtxResXvsRecoTrkHist->GetXaxis()->CenterTitle(); vtxResXvsRecoTrkHist->GetYaxis()->CenterTitle();
    
    TH2D *vtxResYvsRecoTrkHist = new TH2D("vtxResYvsRecoTrk","",31,-0.5,30.5,200,-1,1);
    vtxResYvsRecoTrkHist->Draw("COLZ");
    vtxResYvsRecoTrkHist->SetTitle("Vertex Resolution Y vs Reconstructed Tracks");
    vtxResYvsRecoTrkHist->GetXaxis()->SetTitle("N_{RC}");
    vtxResYvsRecoTrkHist->GetYaxis()->SetTitle("recVtx_y - mcVtx_y (in mm)");
    vtxResYvsRecoTrkHist->GetXaxis()->CenterTitle(); vtxResYvsRecoTrkHist->GetYaxis()->CenterTitle();
  
    TH2D *vtxResZvsRecoTrkHist = new TH2D("vtxResZvsRecoTrk","",31,-0.5,30.5,200,-1,1);
    vtxResZvsRecoTrkHist->Draw("COLZ");
    vtxResZvsRecoTrkHist->SetTitle("Vertex Resolution Z vs Reconstructed Tracks");
    vtxResZvsRecoTrkHist->GetXaxis()->SetTitle("N_{RC}");
    vtxResZvsRecoTrkHist->GetYaxis()->SetTitle("recVtx_z - mcVtx_z (in mm)");
    vtxResZvsRecoTrkHist->GetXaxis()->CenterTitle(); vtxResZvsRecoTrkHist->GetYaxis()->CenterTitle();
  
    TH1D *numGenTrkswithVtxHist = new TH1D("numGenTrkswithVtx","",31,-0.5,30.5);
    TH1D *numRecoTrkswithVtxHist = new TH1D("numRecoTrkswithVtx","",31,-0.5,30.5);


    int counter(0);

    //Loop over events
    std::cout<<"Analyzing "<<mychain->GetEntries()<<" events!"<<std::endl;
    while (tree_reader.Next()) {

	    if(counter%100==0) std::cout<<"Analyzing event "<<counter<<std::endl;
	    
	    counter++;

    	//////////////////////////////////////////////////////////////////////////
    	////////////////////////  Analyze MC Tracks  /////////////////////////
    	//////////////////////////////////////////////////////////////////////////
    
    	//Finding MC vertex using scattered electron
    	TVector3 mcEvtVtx(-999., -999., -999.);    
    	for(unsigned int i=0; i<mcGenStat.GetSize(); i++)
    	{
		if(mcGenStat[i] != 1) continue;
		if(mcPDG[i] != 11) continue;

		bool scatEfound = false;
		// mcParentBegin and mcParentEnd specify the entries from _MCParticles_parents.index 
		// _MCParticles_parents.index stores the MCParticle index
		for(unsigned int j=mcParentBegin[i]; j<mcParentEnd[i]; j++)
		{
          	int parentPDG = mcPDG[mcParentIndex[j]];
          	if(parentPDG == 11) scatEfound = true;
        	}
       
		if(scatEfound == false) continue;
		//Scattered electron found
		double vtx_mc_x =  mcVtxX[i];
		double vtx_mc_y =  mcVtxY[i];
		double vtx_mc_z =  mcVtxZ[i];
		mcEvtVtx = TVector3(vtx_mc_x, vtx_mc_y, vtx_mc_z);
    	}
    	genVtxYvsXHist->Fill(mcEvtVtx.x(), mcEvtVtx.y());
    	TVector3 mcRadius(mcEvtVtx.x(), mcEvtVtx.y(), 0);
    	genVtxRvsZHist->Fill(mcEvtVtx.z(), mcRadius.Mag());
    
    	//Filtering MC Tracks
    	int numMCTracks=0;
    	for(unsigned int i=0; i<mcGenStat.GetSize(); i++)
    	{
		if(mcGenStat[i] != 1) continue;
		if(mcCharge[i] == 0) continue;
	
		TVector3 mcPartVtx(mcVtxX[i], mcVtxY[i], mcVtxZ[i]);
		TVector3 vtx_diff = mcPartVtx - mcEvtVtx;
		if(vtx_diff.Mag() > 1e-4) continue;
	
		TVector3 mcPartMom(mcMomX[i], mcMomY[i], mcMomZ[i]);
		if(fabs(mcPartMom.Eta()) > 3.5) continue;
	
		numMCTracks++;
    	}
    	numGenTracksHist->Fill(numMCTracks);

    	//////////////////////////////////////////////////////////////////////////
    	//////////////////////  Analyze Reconstructed Tracks  //////////////////////
    	//////////////////////////////////////////////////////////////////////////

    	numRecoTracksHist->Fill(recoType.GetSize());
	  
    	//Finding Reconstructed Vertex and Vertexing Efficiency
    	int nVtx=0;
    	float diff=999.;
    	int nAssoPart=0;
    	TVector3 recoEvtVtx(-999., -999., -999.);    
    	for(unsigned int i=0; i<recoVtxType.GetSize(); i++)
    	{
		nVtx++;
	
		TVector3 recoVtx(recoVtxX[i], recoVtxY[i], recoVtxZ[i]);
	
		//Finding the reconstructed vertex closer to the MC vertex
		TVector3 vtx_diff = recoVtx - mcEvtVtx;
		if(vtx_diff.Mag() < diff)
		{
	    	diff = vtx_diff.Mag();
	    	recoEvtVtx = recoVtx;
	    
	    	for(unsigned int j=assoPartBegin[i]; j<assoPartEnd[i]; j++)
	    	{
              	nAssoPart = j;
            	}
		}
    	}
    
    	recoVtxEffHist->Fill(nVtx);
    	recoVtxYvsXHist->Fill(recoEvtVtx.x(), recoEvtVtx.y());
    	TVector3 recoRadius(recoEvtVtx.x(), recoEvtVtx.y(), 0);
    	recoVtxRvsZHist->Fill(recoEvtVtx.z(), recoRadius.Mag());
    	
    	vtxResXvsGenTrkHist->Fill(numMCTracks, recoEvtVtx.x() - mcEvtVtx.x());
    	vtxResYvsGenTrkHist->Fill(numMCTracks, recoEvtVtx.y() - mcEvtVtx.y());
    	vtxResZvsGenTrkHist->Fill(numMCTracks, recoEvtVtx.z() - mcEvtVtx.z());
    	
    	vtxResXvsRecoTrkHist->Fill(nAssoPart, recoEvtVtx.x() - mcEvtVtx.x());
    	vtxResYvsRecoTrkHist->Fill(nAssoPart, recoEvtVtx.y() - mcEvtVtx.y());
    	vtxResZvsRecoTrkHist->Fill(nAssoPart, recoEvtVtx.z() - mcEvtVtx.z());
    	
    	if(nVtx !=0) {
    	numGenTrkswithVtxHist->Fill(numMCTracks);
    	numRecoTrkswithVtxHist->Fill(recoType.GetSize());
    	
    	recoVsMCTracksHist->Fill(numMCTracks, nAssoPart);} 
		  
  	}
  
    //--------------------------------------------------------------------------------------------------------------------------------------------
    recoVtxEffHist->Scale(100/counter);
    ofile->Write(); // Write histograms to file
    ofile->Close(); // Close output file

}
  

