// Plots the resolution of the reconstructed particles

#include <iostream>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "TCanvas.h"
#include "TStyle.h"

void reconstructionAnalysis(TString inFile             = "/home/simong/EIC/scripts/tagger_inference_new5.root",
                            float   beamEnergy         = 10.0,
                            TString momentumCanvasName = "momentum_resolution.png",
                            TString energyThetaPhiCanvasName = "energy_theta_phi_resolution.png",
                            TString relationCanvasName = "relation_resolution.png",
                            TString resolutionGraphsCanvasName = "resolution_graphs.png") {

    //Set ROOT style    
    gStyle->SetPadLeftMargin(0.1);  // Set left margin
    gStyle->SetPadRightMargin(0.0); // Set right margin
    gStyle->SetPadTopMargin(0.0);   // Set top margin
    gStyle->SetPadBottomMargin(0.1);// Set bottom margin
    gStyle->SetTitleAlign(13);
    gStyle->SetTitleX(0.12);          // Place the title on the top right
    gStyle->SetTitleY(0.985);         // Place the title on the top right
    gStyle->SetTitleSize(0.08, "t");
    gStyle->SetTitleXOffset(1.0);    // Adjust y-axis title offset
    gStyle->SetTitleYOffset(1.0);    // Adjust y-axis title offset
    gStyle->SetOptStat(0);

    ROOT::RDataFrame d0("events",inFile, {"MCParticles","TaggerTrackerReconstructedParticles"});

    auto filterDF = d0.Define("SimParticles", "MCParticles[MCParticles.generatorStatus==1 && MCParticles.PDG==11]")
                      .Filter("SimParticles.size()==1")
                      .Filter("TaggerTrackerReconstructedParticles.size()==1");

    // Plot x,y,z momentum resolution as a function of the MCParticle momentum
    auto momentumDF = filterDF
        .Define("reco_momentum", "TaggerTrackerReconstructedParticles[0].momentum")
        .Define("mc_momentum", "SimParticles[0].momentum")
        .Define("px_rec", "reco_momentum.x")
        .Define("py_rec", "reco_momentum.y")
        .Define("pz_rec", "reco_momentum.z")
        .Define("px_mc", "mc_momentum.x")
        .Define("py_mc", "mc_momentum.y")
        .Define("pz_mc", "mc_momentum.z")
        // Calculate theta and phi for reco and MC
        .Define("theta_rec", "std::atan2(std::sqrt(px_rec*px_rec + py_rec*py_rec), pz_rec)")
        .Define("theta_mc", "std::atan2(std::sqrt(px_mc*px_mc + py_mc*py_mc), pz_mc)")
        .Define("phi_rec_rad", "std::atan2(py_rec, px_rec)")
        .Define("phi_mc_rad", "std::atan2(py_mc, px_mc)")
        .Define("phi_rec", "TMath::RadToDeg()*phi_rec_rad")
        .Define("phi_mc", "TMath::RadToDeg()*phi_mc_rad")
        // Calculate resolutions
        .Define("px_diff", "(px_rec - px_mc)")
        .Define("py_diff", "(py_rec - py_mc)")
        .Define("pz_res", "(pz_rec - pz_mc)/pz_mc")
        .Define("E_rec", "TaggerTrackerReconstructedParticles[0].energy")
        .Define("E_mc", "std::sqrt(px_mc*px_mc + py_mc*py_mc + pz_mc*pz_mc + SimParticles[0].mass*SimParticles[0].mass)")
        .Define("E_res", "(E_rec - E_mc)/E_mc")
        .Define("theta_diff", "(theta_rec - theta_mc)")
        .Define("phi_diff", "TMath::RadToDeg()*ROOT::VecOps::DeltaPhi(phi_rec_rad, phi_mc_rad)");

    //Print the size of the original DataFrame
    std::cout << "Original DataFrame size: " << d0.Count().GetValue() << std::endl;
    //Print the size of the filtered DataFrame
    std::cout << "Filtered DataFrame size: " << filterDF.Count().GetValue() << std::endl;

    int   momentumBins      = 100;
    float momentumXRange[2] = {-0.1, 0.1};
    float momentumYRange[2] = {-0.1, 0.1};
    float momentumZRange[2] = {-beamEnergy, 0};

    int   momentumResolutionBins      = 100;
    float momentumDifferenceXRange[2] = {-0.05, 0.05};
    float momentumDifferenceYRange[2] = {-0.05, 0.05};
    float momentumResolutionZRange[2] = {-0.1, 0.1};

    int   energyBins      = 100;
    int   thetaBins       = 100;
    int   phiBins         = 100;
    float energyRange[2]  = {3.5, beamEnergy}; // GeV
    float thetaRange[2]   = {3.134, TMath::Pi()}; // radians from 3.1 to pi
    float phiRange[2]     = {-180, 180}; // degrees from -180 to 180

    int   resolutionBins           = 100;
    float energyResolutionRange[2] = {-0.1, 0.1};
    float thetaResolutionRange[2]  = {-0.003, 0.003};
    float phiResolutionRange[2]    = {-90, 90}; // degrees from -90 to 90

    // Plot reconstructed vs montecarlo momentum components
    auto px_Hist = momentumDF.Histo2D({"px_vs_px", "Reconstructed vs MC px; px reconstructed [GeV]; px MC [GeV]", momentumBins, momentumXRange[0], momentumXRange[1], momentumBins, momentumXRange[0], momentumXRange[1]}, "px_rec", "px_mc");
    auto py_Hist = momentumDF.Histo2D({"py_vs_py", "Reconstructed vs MC py; py reconstructed [GeV]; py MC [GeV]", momentumBins, momentumYRange[0], momentumYRange[1], momentumBins, momentumYRange[0], momentumYRange[1]}, "py_rec", "py_mc");
    auto pz_Hist = momentumDF.Histo2D({"pz_vs_pz", "Reconstructed vs MC pz; pz reconstructed [GeV]; pz MC [GeV]", momentumBins, momentumZRange[0], momentumZRange[1], momentumBins, momentumZRange[0], momentumZRange[1]}, "pz_rec", "pz_mc");

    // Plot individual momentum resolutions
    auto px_diff_Hist = momentumDF.Histo1D({"px_diff", "px difference; px difference [GeV]; Entries", momentumResolutionBins, momentumDifferenceXRange[0], momentumDifferenceXRange[1]}, "px_diff");
    auto py_diff_Hist = momentumDF.Histo1D({"py_diff", "py difference; py difference [GeV]; Entries", momentumResolutionBins, momentumDifferenceYRange[0], momentumDifferenceYRange[1]}, "py_diff");
    auto pz_res_Hist  = momentumDF.Histo1D({"pz_res",  "pz resolution; pz resolution [GeV]; Entries", momentumResolutionBins, momentumResolutionZRange[0], momentumResolutionZRange[1]}, "pz_res");

    // Plot reconstructed vs montecarlo energy, theta and phi
    auto E_Hist     = momentumDF.Histo2D({"E_vs_E",         "Reconstructed vs MC energy; E reconstructed [GeV]; E MC [GeV]",        energyBins, energyRange[0], energyRange[1], energyBins, energyRange[0], energyRange[1]}, "E_rec", "E_mc");
    auto theta_Hist = momentumDF.Histo2D({"theta_vs_theta", "Reconstructed vs MC theta; theta reconstructed [rad]; theta MC [rad]", thetaBins, thetaRange[0], thetaRange[1], thetaBins, thetaRange[0], thetaRange[1]}, "theta_rec", "theta_mc");
    auto phi_Hist   = momentumDF.Histo2D({"phi_vs_phi",     "Reconstructed vs MC phi; phi reconstructed [deg]; phi MC [deg]",       phiBins, phiRange[0], phiRange[1], phiBins, phiRange[0], phiRange[1]}, "phi_rec", "phi_mc");

    auto E_res_Hist      = momentumDF.Histo1D({"E_res", "E resolution; E resolution [GeV]; Entries", resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "E_res");
    auto theta_diff_Hist = momentumDF.Histo1D({"theta_diff", "theta difference; theta difference [rad]; Entries", resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "theta_diff");
    auto phi_diff_Hist   = momentumDF.Histo1D({"phi_diff", "phi difference; phi difference [deg]; Entries", resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "phi_diff");

    // Plot Reconstructed energy, theta and phi resolutions as a function of each mc value of energy, thata and phi
    auto E_res_vs_E_Hist          = momentumDF.Histo2D({"E_res_vs_E", "E resolution vs E MC; E MC [GeV]; E resolution [GeV]", energyBins, energyRange[0], energyRange[1], resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "E_mc", "E_res");
    auto E_res_vs_theta_Hist      = momentumDF.Histo2D({"E_res_vs_theta", "E resolution vs theta MC; theta MC [rad]; E resolution [GeV]", thetaBins, thetaRange[0], thetaRange[1], resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "theta_mc", "E_res");
    auto E_res_vs_phi_Hist        = momentumDF.Histo2D({"E_res_vs_phi", "E resolution vs phi MC; phi MC [deg]; E resolution [GeV]", phiBins, phiRange[0], phiRange[1], resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "phi_mc", "E_res");
    auto theta_diff_vs_E_Hist     = momentumDF.Histo2D({"theta_diff_vs_E", "theta difference vs E MC; E MC [GeV]; theta difference [rad]", energyBins, energyRange[0], energyRange[1], resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "E_mc", "theta_diff");
    auto theta_diff_vs_theta_Hist = momentumDF.Histo2D({"theta_diff_vs_theta", "theta difference vs theta MC; theta MC [rad]; theta difference [rad]", thetaBins, thetaRange[0], thetaRange[1], resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "theta_mc", "theta_diff");
    auto theta_diff_vs_phi_Hist   = momentumDF.Histo2D({"theta_diff_vs_phi", "theta difference vs phi MC; phi MC [deg]; theta difference [rad]", phiBins, phiRange[0], phiRange[1], resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "phi_mc", "theta_diff");
    auto phi_diff_vs_E_Hist       = momentumDF.Histo2D({"phi_diff_vs_E", "phi difference vs E MC; E MC [GeV]; phi difference [rad]", energyBins, energyRange[0], energyRange[1], resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "E_mc", "phi_diff");
    auto phi_diff_vs_theta_Hist   = momentumDF.Histo2D({"phi_diff_vs_theta", "phi difference vs theta MC; theta MC [rad]; phi difference [deg]", thetaBins, thetaRange[0], thetaRange[1], resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "theta_mc", "phi_diff");
    auto phi_diff_vs_phi_Hist     = momentumDF.Histo2D({"phi_diff_vs_phi", "phi difference vs phi MC; phi MC [deg]; phi difference [deg]", phiBins, phiRange[0], phiRange[1], resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "phi_mc", "phi_diff");

    // Create canvas for momentum component plots
    TCanvas *cMomentum = new TCanvas("momentum_canvas", "Momentum Resolution", 3000, 1600);
    cMomentum->Divide(3, 2);
    cMomentum->cd(1);
    px_Hist->Draw("colz");
    cMomentum->cd(2);
    py_Hist->Draw("colz");
    cMomentum->cd(3);
    pz_Hist->Draw("colz");
    cMomentum->cd(4);
    px_diff_Hist->Draw();
    cMomentum->cd(5);
    py_diff_Hist->Draw();
    cMomentum->cd(6);
    pz_res_Hist->Draw();
    cMomentum->SetGrid();
    cMomentum->Update();
    // Save the canvas as a PNG file
    cMomentum->SaveAs(momentumCanvasName);

    // Create canvas for energy, theta and phi resolution plots
    TCanvas *cEnergyThetaPhi = new TCanvas("energy_theta_phi_canvas", "Energy, Theta and Phi Resolution", 3000, 1600);
    cEnergyThetaPhi->Divide(3, 2);
    cEnergyThetaPhi->cd(1);
    E_Hist->Draw("colz");
    cEnergyThetaPhi->cd(2);
    theta_Hist->Draw("colz");
    cEnergyThetaPhi->cd(3);
    phi_Hist->Draw("colz");
    cEnergyThetaPhi->cd(4);
    E_res_Hist->Draw();
    cEnergyThetaPhi->cd(5);
    theta_diff_Hist->Draw();
    cEnergyThetaPhi->cd(6);
    phi_diff_Hist->Draw();
    cEnergyThetaPhi->SetGrid();
    cEnergyThetaPhi->Update();
    // Save the canvas as a PNG file
    cEnergyThetaPhi->SaveAs(energyThetaPhiCanvasName);
    
    // Create canvas for resolution vs MC values
    TCanvas *cResolutionVsMC = new TCanvas("resolution_vs_mc_canvas", "Resolution vs MC Values", 3000, 1600);
    cResolutionVsMC->Divide(3, 3);
    cResolutionVsMC->cd(1);
    E_res_vs_E_Hist->Draw("colz");
    cResolutionVsMC->cd(2);
    E_res_vs_theta_Hist->Draw("colz");
    cResolutionVsMC->cd(3);
    E_res_vs_phi_Hist->Draw("colz");
    cResolutionVsMC->cd(4);
    theta_diff_vs_E_Hist->Draw("colz");
    cResolutionVsMC->cd(5);
    theta_diff_vs_theta_Hist->Draw("colz");
    cResolutionVsMC->cd(6);
    theta_diff_vs_phi_Hist->Draw("colz");
    cResolutionVsMC->cd(7);
    phi_diff_vs_E_Hist->Draw("colz");
    cResolutionVsMC->cd(8);
    phi_diff_vs_theta_Hist->Draw("colz");
    cResolutionVsMC->cd(9);
    phi_diff_vs_phi_Hist->Draw("colz");
    cResolutionVsMC->SetGrid();
    cResolutionVsMC->Update();
    // Save the canvas as a PNG file
    cResolutionVsMC->SaveAs(relationCanvasName);

    // Fit Gaussians to the E vs E histogram slices
    TObjArray* fitE_vs_E = new TObjArray();
    TObjArray* fitTheta_vs_E = new TObjArray();
    TObjArray* fitPhi_vs_theta = new TObjArray();
    E_res_vs_E_Hist->FitSlicesY(nullptr, 1, -1, 0, "Q", fitE_vs_E);
    theta_diff_vs_E_Hist->FitSlicesY(nullptr, 1, -1, 0, "Q", fitTheta_vs_E);
    phi_diff_vs_theta_Hist->FitSlicesY(nullptr, 1, -1, 0, "Q", fitPhi_vs_theta);

    // Create graphs containing the fitted means and standard deviations
    TH1* hE_vs_E_mean         = (TH1*)fitE_vs_E->At(1);     // mean values (index 1)
    TH1* hE_vs_E_stddev       = (TH1*)fitE_vs_E->At(2);     // stddev values (index 2)
    TH1* hTheta_vs_E_mean     = (TH1*)fitTheta_vs_E->At(1); // mean values (index 1)
    TH1* hTheta_vs_E_stddev   = (TH1*)fitTheta_vs_E->At(2); // stddev values (index 2)
    TH1* hPhi_vs_theta_mean   = (TH1*)fitPhi_vs_theta->At(1);   // mean values (index 1)
    TH1* hPhi_vs_theta_stddev = (TH1*)fitPhi_vs_theta->At(2);   // stddev values (index 2)

    // Create a canvas for the resolution graphs
    TCanvas *cResolutionGraphs = new TCanvas("resolution_graphs_canvas", "Resolution Graphs", 1200, 800);
    cResolutionGraphs->Divide(3, 2);
    cResolutionGraphs->cd(1);
    hE_vs_E_mean->SetTitle("Mean Energy Offset vs E MC; Energy MC [GeV]; Mean Energy Offset [GeV]");
    hE_vs_E_mean->SetMarkerStyle(20);
    hE_vs_E_mean->SetMarkerColor(kBlue);
    hE_vs_E_mean->SetMaximum(0.01); // Adjust maximum for better visibility
    hE_vs_E_mean->SetMinimum(-0.01); // Adjust minimum for better visibility
    hE_vs_E_mean->Draw();
    cResolutionGraphs->cd(4);
    hE_vs_E_stddev->SetTitle("Energy Resolution vs E MC; Energy MC [GeV]; Energy Resolution [GeV]");
    hE_vs_E_stddev->SetMarkerStyle(20);
    hE_vs_E_stddev->SetMarkerColor(kRed);
    hE_vs_E_stddev->SetMaximum(0.04); // Adjust maximum for better visibility
    hE_vs_E_stddev->SetMinimum(0); // Adjust minimum for better visibility
    hE_vs_E_stddev->Draw();
    cResolutionGraphs->cd(2);
    hTheta_vs_E_mean->SetTitle("Mean Theta Offset vs E MC; Energy MC [GeV]; Mean Theta Offset [rad]");
    hTheta_vs_E_mean->SetMarkerStyle(20);
    hTheta_vs_E_mean->SetMarkerColor(kBlue);
    hTheta_vs_E_mean->SetMaximum(0.0003); // Adjust maximum for better visibility
    hTheta_vs_E_mean->SetMinimum(-0.0003); // Adjust minimum for better visibility
    hTheta_vs_E_mean->Draw();
    cResolutionGraphs->cd(5);
    hTheta_vs_E_stddev->SetTitle("Std Dev Theta Offset vs E MC; Energy MC [GeV]; Std Dev Theta Offset [rad]");
    hTheta_vs_E_stddev->SetMarkerStyle(20);
    hTheta_vs_E_stddev->SetMarkerColor(kRed);
    hTheta_vs_E_stddev->SetMaximum(0.0005); // Adjust maximum for better visibility
    hTheta_vs_E_stddev->SetMinimum(0); // Adjust minimum for better visibility
    hTheta_vs_E_stddev->Draw();
    cResolutionGraphs->cd(3);
    hPhi_vs_theta_mean->SetTitle("Mean Phi Offset vs theta MC; theta MC [rad]; Mean Phi Offset [rad]");
    hPhi_vs_theta_mean->SetMarkerStyle(20);
    hPhi_vs_theta_mean->SetMarkerColor(kBlue);
    hPhi_vs_theta_mean->SetMaximum(5); // Adjust maximum for better visibility
    hPhi_vs_theta_mean->SetMinimum(-5); // Adjust minimum for better visibility
    hPhi_vs_theta_mean->Draw();
    cResolutionGraphs->cd(6);
    hPhi_vs_theta_stddev->SetTitle("Std Dev Phi Offset vs theta MC; theta MC [rad]; Std Dev Phi Offset [rad]");
    hPhi_vs_theta_stddev->SetMarkerStyle(20);
    hPhi_vs_theta_stddev->SetMarkerColor(kRed);
    hPhi_vs_theta_stddev->SetMaximum(60); // Adjust maximum for better visibility
    hPhi_vs_theta_stddev->SetMinimum(0); // Adjust minimum for better visibility
    hPhi_vs_theta_stddev->Draw();
    cResolutionGraphs->SetGrid();
    cResolutionGraphs->Update();
    // Save the canvas as a PNG file
    cResolutionGraphs->SaveAs(resolutionGraphsCanvasName);

    // Check to see if resolutions pass tests.


}

