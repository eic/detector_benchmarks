// Plots the resolution of the reconstructed particles

#include <iostream>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

void reconstructionAnalysis( TString inFile                       = "/home/simong/EIC/detector_benchmarks_anl/sim_output/lowq2_reconstruction/analysis/Low-Q2_retrained_Particles_new.eicrecon.edm4hep.root",
                             float   electronEnergy               = 18.0,
                             float   protonEnergy                 = 275.0,
                             TString outFile                      = "reconstruction_results.root",
                             TString momentumCanvasName           = "momentum_resolution.png",
                             TString energyThetaPhiCanvasName     = "energy_theta_phi_resolution.png",
                             TString relationCanvasName           = "relation_resolution.png",
                             TString resolutionGraphsCanvasName   = "resolution_graphs.png",
                             TString acceptanceCanvasName         = "acceptance_canvas.png",
                             TString energyQ2acceptanceCanvasName = "energy_Q2_acceptance_canvas.png",
                             TString xQ2acceptanceCanvasName      = "x_Q2_acceptance_canvas.png",
                             TString Q2xResolutionCanvasName      = "Q2_x_resolution.png",
                             std::string particleCollectionName   = "TaggerTrackerReconstructedParticles") {

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

    ROOT::RDataFrame d0("events",inFile, {"MCParticles",particleCollectionName});

    // Base MC electron definitions prior to reconstruction filter
    auto baseDF = d0.Define("SimParticles", "MCParticles[MCParticles.generatorStatus==1 && MCParticles.PDG==11]")
                    .Define("hasElectron", "SimParticles.size()==1")
                    .Define("mc_momentum", "hasElectron ? SimParticles[0].momentum : MCParticles[0].momentum")
                    .Define("px_mc", "mc_momentum.x")
                    .Define("py_mc", "mc_momentum.y")
                    .Define("pz_mc", "mc_momentum.z")
                    .Define("E_mc", "sqrt(px_mc*px_mc + py_mc*py_mc + pz_mc*pz_mc + (hasElectron ? SimParticles[0].mass*SimParticles[0].mass : MCParticles[0].mass*MCParticles[0].mass))")
                    .Define("theta_mc", "std::atan2(std::sqrt(px_mc*px_mc + py_mc*py_mc), pz_mc)")
                    .Define("phi_mc_rad", "std::atan2(py_mc, px_mc)")
                    .Define("phi_mc", "TMath::RadToDeg()*phi_mc_rad")
                    // Invariant kinematics (head-on ep collider):
                    //   Initial electron: k0 = (Ee, 0,0, -|pe|)
                    //   Initial proton:  P  = (Ep, 0,0, +|pp|)
                    //   Scattered electron: k' = (E, px, py, pz)
                    //   Four-momentum transfer: q = k0 - k'
                    // Definitions:
                    //   Q2 = -q^2  (positive)
                    //   Bjorken x = Q2 / (2 P·q)
                    //   W^2 = M_p^2 + 2 P·q - Q2  (invariant mass of hadronic final state)
                    // These formulas correctly incorporate the proton beam energy in the CM frame via P·q.
                    // Q2 using initial electron beam along -z, full 4-vector definition
                    .Define("Q2_mc", [electronEnergy](double px, double py, double pz, double E){
                        const double me = 0.00051099895; // GeV
                        const double Ee = electronEnergy;
                        const double pe = std::sqrt(std::max(0.0, Ee*Ee - me*me));
                        const double dE = Ee - E;
                        const double dqx = -px;
                        const double dqy = -py;
                        const double dqz = -pe - pz; // initial e momentum along -z
                        const double q2 = dE*dE - (dqx*dqx + dqy*dqy + dqz*dqz);
                        return -q2; // ensure positive Q2
                    }, {"px_mc","py_mc","pz_mc","E_mc"})
                    // W using invariants: W^2 = Mp^2 + 2 P·q - Q^2
                    .Define("W_mc", [electronEnergy, protonEnergy](double Q2, double px, double py, double pz, double E){
                        const double Mp = 0.9382720813; // GeV
                        const double me = 0.00051099895; // GeV
                        const double Ee = electronEnergy;
                        const double Ep = protonEnergy;
                        const double pe = std::sqrt(std::max(0.0, Ee*Ee - me*me));
                        const double pp = std::sqrt(std::max(0.0, Ep*Ep - Mp*Mp));
                        // q = k0 - k'
                        const double qE = Ee - E;
                        const double qx = -px;
                        const double qy = -py;
                        const double qz = -pe - pz; // initial e along -z
                        // P·q with P0 = (Ep,0,0,+pp)
                        const double Pdotq = Ep*qE - (0.0*qx + 0.0*qy + pp*qz);
                        const double W2 = Mp*Mp + 2.0*Pdotq - Q2;
                        return W2 > 0 ? std::sqrt(W2) : 0.0;
                    }, {"Q2_mc","px_mc","py_mc","pz_mc","E_mc"})
                    .Define("log10Q2_mc", "Q2_mc>0 ? std::log10(Q2_mc) : -10.0")
                    // Bjorken x = Q2 / (2 P·q)
                    .Define("x_bj", [electronEnergy, protonEnergy](double Q2, double px, double py, double pz, double E){
                        const double Mp = 0.9382720813; // GeV
                        const double me = 0.00051099895; // GeV
                        const double Ee = electronEnergy;
                        const double Ep = protonEnergy;
                        const double pe = std::sqrt(std::max(0.0, Ee*Ee - me*me));
                        const double pp = std::sqrt(std::max(0.0, Ep*Ep - Mp*Mp));
                        const double qE = Ee - E;
                        const double qx = -px;
                        const double qy = -py;
                        const double qz = -pe - pz;
                        const double Pdotq = Ep*qE - (0.0*qx + 0.0*qy + pp*qz);
                        const double denom = 2.0 * Pdotq;
                        return (Q2 > 0.0 && denom > 0.0) ? (Q2 / denom) : 1e-10;
                    }, {"Q2_mc","px_mc","py_mc","pz_mc","E_mc"})
                    .Define("log10x_mc", "x_bj>0 ? std::log10(x_bj) : -10.0");

    // Denominator: events with a valid scattered electron
    auto denomDF = baseDF.Filter("hasElectron");
    // Numerator: additionally require single reconstructed particle
    auto filterDF = denomDF.Filter(particleCollectionName+".size()==1");

    // Plot x,y,z momentum resolution as a function of the MCParticle momentum
    auto momentumDF = filterDF
        .Define("reco_momentum", particleCollectionName+"[0].momentum")
        .Define("px_rec", "reco_momentum.x")
        .Define("py_rec", "reco_momentum.y")
        .Define("pz_rec", "reco_momentum.z")
        // Calculate theta and phi for reco
        .Define("theta_rec", "std::atan2(std::sqrt(px_rec*px_rec + py_rec*py_rec), pz_rec)")
        .Define("phi_rec_rad", "std::atan2(py_rec, px_rec)")
        .Define("phi_rec", "TMath::RadToDeg()*phi_rec_rad")
        // Calculate resolutions
        .Define("px_diff", "(px_rec - px_mc)")
        .Define("py_diff", "(py_rec - py_mc)")
        .Define("pz_res", "(pz_rec - pz_mc)/pz_mc")
        .Define("E_rec", particleCollectionName+"[0].energy")
        .Define("E_res", "(E_rec - E_mc)/E_mc")
        .Define("theta_diff", "(theta_rec - theta_mc)")
        .Define("phi_diff", "TMath::RadToDeg()*ROOT::VecOps::DeltaPhi(phi_rec_rad, phi_mc_rad)")
        // Reconstructed Q2 and x using reconstructed kinematics and full collider invariants
        .Define("Q2_rec", [electronEnergy](float px, float py, float pz, float E){
            const double me = 0.00051099895; // GeV
            const double Ee = electronEnergy;
            const double pe = std::sqrt(std::max(0.0, Ee*Ee - me*me));
            const double dE = Ee - E;
            const double dqx = -px;
            const double dqy = -py;
            const double dqz = -pe - pz;
            const double q2 = (dE*dE - (dqx*dqx + dqy*dqy + dqz*dqz));
            return -q2;
        }, {"px_rec","py_rec","pz_rec","E_rec"})
        .Define("x_bj_rec", [electronEnergy, protonEnergy](double Q2, float px, float py, float pz, float E){
            const double Mp = 0.9382720813; // GeV
            const double me = 0.00051099895; // GeV
            const double Ee = electronEnergy;
            const double Ep = protonEnergy;
            const double pe = std::sqrt(std::max(0.0, Ee*Ee - me*me));
            const double pp = std::sqrt(std::max(0.0, Ep*Ep - Mp*Mp));
            const double qE = Ee - E;
            const double qx = -px;
            const double qy = -py;
            const double qz = -pe - pz;
            const double Pdotq = Ep*qE - (0.0*qx + 0.0*qy + pp*qz);
            const double denom = 2.0 * Pdotq;
            return (Q2 > 0.0 && denom > 0.0) ? (Q2 / denom) : 1e-10;
        }, {"Q2_rec","px_rec","py_rec","pz_rec","E_rec"})
        .Define("log10Q2_rec", "Q2_rec>0 ? std::log10(Q2_rec) : -10.0")
        .Define("log10x_rec", "x_bj_rec>0 ? std::log10(x_bj_rec) : -10.0")
        .Define("Q2_diff", "(Q2_rec - Q2_mc)")
        .Define("x_diff", "(x_bj_rec - x_bj)")
        .Define("log10Q2_diff", "(log10Q2_rec - log10Q2_mc)")
        .Define("log10x_diff", "(log10x_rec - log10x_mc)");

    //Print the size of the original DataFrame
    std::cout << "Original DataFrame size: " << d0.Count().GetValue() << std::endl;
    //Print the size of the filtered DataFrame
    std::cout << "Filtered DataFrame size: " << filterDF.Count().GetValue() << std::endl;

    int   momentumBins      = 100;
    float momentumXRange[2] = {-0.1, 0.1};
    float momentumYRange[2] = {-0.1, 0.1};
    float momentumZRange[2] = {-electronEnergy, 0.0f};

    int   momentumResolutionBins      = 100;
    float momentumDifferenceXRange[2] = {-0.05, 0.05};
    float momentumDifferenceYRange[2] = {-0.05, 0.05};
    float momentumResolutionZRange[2] = {-0.1, 0.1};

    int   energyBins      = 100;
    int   thetaBins       = 100;
    int   phiBins         = 100;
    float energyRange[2]  = {0.0f, electronEnergy};
    float thetaRange[2]   = {3.134, TMath::Pi()}; // radians from 3.1 to pi
    float phiRange[2]     = {-180, 180}; // degrees from -180 to 180

    int   resolutionBins           = 100;
    float energyResolutionRange[2] = {-0.1, 0.1};
    float thetaResolutionRange[2]  = {-0.003, 0.003};
    float phiResolutionRange[2]    = {-90, 90}; // degrees from -90 to 90

    // Q2 and x binning configuration
    int   Q2Bins           = 100;
    float Q2Range[2]       = {0.0, 1.0};
    int   xBins            = 100;
    float xRange[2]        = {0.0, 0.1};
    int   log10Q2Bins      = 100;
    float log10Q2Range[2]  = {-10.0, 0.0};
    int   log10xBins       = 100;
    float log10xRange[2]   = {-13.0, 0.0};
    int   Q2DiffBins       = 100;
    float Q2DiffRange[2]   = {-0.2, 0.2};
    int   xDiffBins        = 100;
    float xDiffRange[2]    = {-0.02, 0.02};
    int   log10DiffBins    = 100;
    float log10Q2DiffRange[2] = {-2.0, 2.0};
    float log10xDiffRange[2]  = {-2.0, 2.0};

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

    // Plot Reconstructed energy, theta and phi resolutions as a function of each reconstructed value of energy, thata and phi
    auto E_res_vs_E_Hist          = momentumDF.Histo2D({"E_res_vs_E", "E resolution vs E reconstructed; E reconstructed [GeV]; E resolution [GeV]", energyBins, energyRange[0], energyRange[1], resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "E_rec", "E_res");
    auto E_res_vs_theta_Hist      = momentumDF.Histo2D({"E_res_vs_theta", "E resolution vs theta reconstructed; theta reconstructed [rad]; E resolution [GeV]", thetaBins, thetaRange[0], thetaRange[1], resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "theta_rec", "E_res");
    auto E_res_vs_phi_Hist        = momentumDF.Histo2D({"E_res_vs_phi", "E resolution vs phi reconstructed; phi reconstructed [deg]; E resolution [GeV]", phiBins, phiRange[0], phiRange[1], resolutionBins, energyResolutionRange[0], energyResolutionRange[1]}, "phi_rec", "E_res");
    auto theta_diff_vs_E_Hist     = momentumDF.Histo2D({"theta_diff_vs_E", "theta difference vs E reconstructed; E reconstructed [GeV]; theta difference [rad]", energyBins, energyRange[0], energyRange[1], resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "E_rec", "theta_diff");
    auto theta_diff_vs_theta_Hist = momentumDF.Histo2D({"theta_diff_vs_theta", "theta difference vs theta reconstructed; theta reconstructed [rad]; theta difference [rad]", thetaBins, thetaRange[0], thetaRange[1], resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "theta_rec", "theta_diff");
    auto theta_diff_vs_phi_Hist   = momentumDF.Histo2D({"theta_diff_vs_phi", "theta difference vs phi reconstructed; phi reconstructed [deg]; theta difference [rad]", phiBins, phiRange[0], phiRange[1], resolutionBins, thetaResolutionRange[0], thetaResolutionRange[1]}, "phi_rec", "theta_diff");
    auto phi_diff_vs_E_Hist       = momentumDF.Histo2D({"phi_diff_vs_E", "phi difference vs E reconstructed; E reconstructed [GeV]; phi difference [rad]", energyBins, energyRange[0], energyRange[1], resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "E_rec", "phi_diff");
    auto phi_diff_vs_theta_Hist   = momentumDF.Histo2D({"phi_diff_vs_theta", "phi difference vs theta reconstructed; theta reconstructed [rad]; phi difference [deg]", thetaBins, thetaRange[0], thetaRange[1], resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "theta_rec", "phi_diff");
    auto phi_diff_vs_phi_Hist     = momentumDF.Histo2D({"phi_diff_vs_phi", "phi difference vs phi reconstructed; phi reconstructed [deg]; phi difference [deg]", phiBins, phiRange[0], phiRange[1], resolutionBins, phiResolutionRange[0], phiResolutionRange[1]}, "phi_rec", "phi_diff");

    // Q2 and x resolution plots
    auto Q2_Hist = momentumDF.Histo2D({"Q2_vs_Q2", "Reconstructed vs MC Q^{2}; Q^{2} reconstructed [GeV^{2}]; Q^{2} MC [GeV^{2}]", Q2Bins, Q2Range[0], Q2Range[1], Q2Bins, Q2Range[0], Q2Range[1]}, "Q2_rec", "Q2_mc");
    auto x_Hist  = momentumDF.Histo2D({"x_vs_x", "Reconstructed vs MC x; x reconstructed; x MC", xBins, xRange[0], xRange[1], xBins, xRange[0], xRange[1]}, "x_bj_rec", "x_bj");
    auto log10Q2_comp_Hist = momentumDF.Histo2D({"log10Q2_vs_log10Q2", "Reconstructed vs MC log10(Q^{2}); log10(Q^{2}) reconstructed; log10(Q^{2}) MC", log10Q2Bins, log10Q2Range[0], log10Q2Range[1], log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "log10Q2_rec", "log10Q2_mc");
    auto log10x_comp_Hist  = momentumDF.Histo2D({"log10x_vs_log10x", "Reconstructed vs MC log10(x); log10(x) reconstructed; log10(x) MC", log10xBins, log10xRange[0], log10xRange[1], log10xBins, log10xRange[0], log10xRange[1]}, "log10x_rec", "log10x_mc");

    auto Q2_diff_Hist = momentumDF.Histo1D({"Q2_diff", "Q^{2} difference; Q^{2} difference [GeV^{2}]; Entries", Q2DiffBins, Q2DiffRange[0], Q2DiffRange[1]}, "Q2_diff");
    auto x_diff_Hist  = momentumDF.Histo1D({"x_diff", "x difference; x difference; Entries", xDiffBins, xDiffRange[0], xDiffRange[1]}, "x_diff");
    auto log10Q2_diff_Hist = momentumDF.Histo1D({"log10Q2_diff", "log10(Q^{2}) difference; log10(Q^{2}) difference; Entries", log10DiffBins, log10Q2DiffRange[0], log10Q2DiffRange[1]}, "log10Q2_diff");
    auto log10x_diff_Hist  = momentumDF.Histo1D({"log10x_diff", "log10(x) difference; log10(x) difference; Entries", log10DiffBins, log10xDiffRange[0], log10xDiffRange[1]}, "log10x_diff");

    auto Q2_diff_vs_Q2_Hist = momentumDF.Histo2D({"Q2_diff_vs_Q2", "Q^{2} difference vs MC Q^{2}; Q^{2} MC [GeV^{2}]; Q^{2} difference [GeV^{2}]", Q2Bins, Q2Range[0], Q2Range[1], Q2DiffBins, Q2DiffRange[0], Q2DiffRange[1]}, "Q2_mc", "Q2_diff");
    auto x_diff_vs_x_Hist   = momentumDF.Histo2D({"x_diff_vs_x", "x difference vs MC x; x MC; x difference", xBins, xRange[0], xRange[1], xDiffBins, xDiffRange[0], xDiffRange[1]}, "x_bj", "x_diff");
    auto log10Q2_diff_vs_log10Q2_Hist = momentumDF.Histo2D({"log10Q2_diff_vs_log10Q2", "log10(Q^{2}) difference vs MC log10(Q^{2}); log10(Q^{2}) MC; log10(Q^{2}) difference", log10Q2Bins, log10Q2Range[0], log10Q2Range[1], log10DiffBins, log10Q2DiffRange[0], log10Q2DiffRange[1]}, "log10Q2_mc", "log10Q2_diff");
    auto log10x_diff_vs_log10x_Hist   = momentumDF.Histo2D({"log10x_diff_vs_log10x", "log10(x) difference vs MC log10(x); log10(x) MC; log10(x) difference", log10xBins, log10xRange[0], log10xRange[1], log10DiffBins, log10xDiffRange[0], log10xDiffRange[1]}, "log10x_mc", "log10x_diff");

    // Acceptance histograms (denominator vs numerator)
    auto E_all_Hist     = denomDF.Histo1D({"E_all",     "MC Electron Energy; E [GeV]; Entries", energyBins, energyRange[0], energyRange[1]}, "E_mc");
    auto theta_all_Hist = denomDF.Histo1D({"theta_all", "MC Electron Theta; theta [rad]; Entries", thetaBins, thetaRange[0], thetaRange[1]}, "theta_mc");
    auto phi_all_Hist   = denomDF.Histo1D({"phi_all",   "MC Electron Phi; phi [deg]; Entries", phiBins, phiRange[0], phiRange[1]}, "phi_mc");
    auto log10Q2_all_Hist = denomDF.Histo1D({"log10Q2_all", "MC log10(Q^{2}); log10(Q^{2}) [GeV^{2}]; Entries", log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "log10Q2_mc");
    auto W_all_Hist     = denomDF.Histo1D({"W_all",     "MC W; W [GeV]; Entries", energyBins, 0.0, 5.0}, "W_mc");

    auto E_acc_Hist     = filterDF.Histo1D({"E_acc",     "Accepted Electron Energy; E [GeV]; Entries", energyBins, energyRange[0], energyRange[1]}, "E_mc");
    auto theta_acc_Hist = filterDF.Histo1D({"theta_acc", "Accepted Electron Theta; theta [rad]; Entries", thetaBins, thetaRange[0], thetaRange[1]}, "theta_mc");
    auto phi_acc_Hist   = filterDF.Histo1D({"phi_acc",   "Accepted Electron Phi; phi [deg]; Entries", phiBins, phiRange[0], phiRange[1]}, "phi_mc");
    auto log10Q2_acc_Hist = filterDF.Histo1D({"log10Q2_acc", "Accepted log10(Q^{2}); log10(Q^{2}) [GeV^{2}]; Entries", log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "log10Q2_mc");
    auto W_acc_Hist     = filterDF.Histo1D({"W_acc",     "Accepted W; W [GeV]; Entries", energyBins, 0.0, 5.0}, "W_mc");

    // 2D acceptance: log10Q2 vs E
    auto log10Q2_vs_E_all_Hist = denomDF.Histo2D({"log10Q2_vs_E_all", "MC log10(Q^{2}) vs E; E [GeV]; log10(Q^{2})", energyBins, energyRange[0], energyRange[1], log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "E_mc", "log10Q2_mc");
    auto log10Q2_vs_E_acc_Hist = filterDF.Histo2D({"log10Q2_vs_E_acc", "Accepted log10(Q^{2}) vs E; E [GeV]; log10(Q^{2})", energyBins, energyRange[0], energyRange[1], log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "E_mc", "log10Q2_mc");

    // 2D acceptance: log10Q2 vs log10x
    auto log10Q2_vs_log10x_all_Hist = denomDF.Histo2D({"log10Q2_vs_log10x_all", "MC log10(Q^{2}) vs log10(x); log10(x); log10(Q^{2})", log10xBins, log10xRange[0], log10xRange[1], log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "log10x_mc", "log10Q2_mc");
    auto log10Q2_vs_log10x_acc_Hist = filterDF.Histo2D({"log10Q2_vs_log10x_acc", "Accepted log10(Q^{2}) vs log10(x); log10(x); log10(Q^{2})", log10xBins, log10xRange[0], log10xRange[1], log10Q2Bins, log10Q2Range[0], log10Q2Range[1]}, "log10x_mc", "log10Q2_mc");

    TH1D* hE_acceptance     = (TH1D*)E_acc_Hist->Clone("hE_acceptance");     hE_acceptance->Divide(E_all_Hist.GetPtr());
    TH1D* hTheta_acceptance = (TH1D*)theta_acc_Hist->Clone("hTheta_acceptance"); hTheta_acceptance->Divide(theta_all_Hist.GetPtr());
    TH1D* hPhi_acceptance   = (TH1D*)phi_acc_Hist->Clone("hPhi_acceptance");   hPhi_acceptance->Divide(phi_all_Hist.GetPtr());
    TH1D* hLog10Q2_acceptance = (TH1D*)log10Q2_acc_Hist->Clone("hLog10Q2_acceptance"); hLog10Q2_acceptance->Divide(log10Q2_all_Hist.GetPtr());
    TH1D* hW_acceptance     = (TH1D*)W_acc_Hist->Clone("hW_acceptance");     hW_acceptance->Divide(W_all_Hist.GetPtr());

    // 2D acceptance
    TH2D* hLog10Q2_vs_E_acceptance = (TH2D*)log10Q2_vs_E_acc_Hist->Clone("hLog10Q2_vs_E_acceptance");
    hLog10Q2_vs_E_acceptance->Divide(log10Q2_vs_E_all_Hist.GetPtr());

    TH2D* hLog10Q2_vs_log10x_acceptance = (TH2D*)log10Q2_vs_log10x_acc_Hist->Clone("hLog10Q2_vs_log10x_acceptance");
    hLog10Q2_vs_log10x_acceptance->Divide(log10Q2_vs_log10x_all_Hist.GetPtr());

    auto styleAcc = [](TH1* h,const char* title){
        h->SetTitle(title);
        h->GetYaxis()->SetTitle("Acceptance");
        h->SetMinimum(0.0);
        h->SetMaximum(1.05);
        h->SetMarkerStyle(20);
        h->SetMarkerColor(kBlue+1);
        h->SetLineColor(kBlue+1);
    };
    styleAcc(hE_acceptance,     "Electron Acceptance vs E; E [GeV]");
    styleAcc(hTheta_acceptance, "Electron Acceptance vs theta; theta [rad]");
    styleAcc(hPhi_acceptance,   "Electron Acceptance vs phi; phi [deg]");
    styleAcc(hLog10Q2_acceptance, "Electron Acceptance vs log10(Q^{2}); log10(Q^{2})");
    styleAcc(hW_acceptance,     "Electron Acceptance vs W; W [GeV]");

    // Style 2D acceptance
    hLog10Q2_vs_E_acceptance->SetTitle("Electron Acceptance: log10(Q^{2}) vs E; E [GeV]; log10(Q^{2})");
    hLog10Q2_vs_E_acceptance->GetZaxis()->SetTitle("Acceptance");
    hLog10Q2_vs_E_acceptance->SetMinimum(0.0);
    hLog10Q2_vs_E_acceptance->SetMaximum(1.0);

    hLog10Q2_vs_log10x_acceptance->SetTitle("Electron Acceptance: log10(Q^{2}) vs log10(x); log10(x); log10(Q^{2})");
    hLog10Q2_vs_log10x_acceptance->GetZaxis()->SetTitle("Acceptance");
    hLog10Q2_vs_log10x_acceptance->SetMinimum(0.0);
    hLog10Q2_vs_log10x_acceptance->SetMaximum(1.0);

    TCanvas *cAcceptance1D = new TCanvas("acceptance_1D_canvas", "Electron Acceptance (1D)", 2400, 1600);
    cAcceptance1D->Divide(3,2);
    cAcceptance1D->cd(1); hE_acceptance->Draw("E1");
    cAcceptance1D->cd(2); hTheta_acceptance->Draw("E1");
    cAcceptance1D->cd(3); hPhi_acceptance->Draw("E1");
    cAcceptance1D->cd(4); hLog10Q2_acceptance->Draw("E1");
    cAcceptance1D->cd(5); hW_acceptance->Draw("E1");
    cAcceptance1D->SetGrid();
    cAcceptance1D->Update();
    cAcceptance1D->SaveAs(acceptanceCanvasName);

    TCanvas *cAcceptance2D = new TCanvas("acceptance_2D_canvas", "Electron Acceptance: log10(Q^{2}) vs E", 1200, 1000);
    cAcceptance2D->cd();
    gStyle->SetOptStat(0);
    hLog10Q2_vs_E_acceptance->Draw("COLZ");
    cAcceptance2D->SetRightMargin(0.15);
    cAcceptance2D->Update();
    cAcceptance2D->SaveAs(energyQ2acceptanceCanvasName);

    TCanvas *cAcceptanceXQ2 = new TCanvas("acceptance_xQ2_canvas", "Electron Acceptance: log10(Q^{2}) vs log10(x)", 1200, 1000);
    cAcceptanceXQ2->cd();
    gStyle->SetOptStat(0);
    hLog10Q2_vs_log10x_acceptance->Draw("COLZ");
    cAcceptanceXQ2->SetRightMargin(0.15);
    cAcceptanceXQ2->Update();
    cAcceptanceXQ2->SaveAs(xQ2acceptanceCanvasName);

    // Create canvas for Q2 and x resolution
    TCanvas *cQ2xResolution = new TCanvas("Q2_x_resolution_canvas", "Q^{2} and x Resolution", 3000, 1600);
    cQ2xResolution->Divide(4, 2);
    cQ2xResolution->cd(1);
    log10Q2_comp_Hist->Draw("colz");
    gPad->SetLogz();
    cQ2xResolution->cd(2);
    log10Q2_diff_Hist->Draw();
    cQ2xResolution->cd(3);
    log10Q2_diff_vs_log10Q2_Hist->Draw("colz");
    gPad->SetLogz();
    cQ2xResolution->cd(4);
    Q2_Hist->Draw("colz");
    gPad->SetLogz();
    cQ2xResolution->cd(5);
    log10x_comp_Hist->Draw("colz");
    gPad->SetLogz();
    cQ2xResolution->cd(6);
    log10x_diff_Hist->Draw();
    cQ2xResolution->cd(7);
    log10x_diff_vs_log10x_Hist->Draw("colz");
    gPad->SetLogz();
    cQ2xResolution->cd(8);
    x_Hist->Draw("colz");
    gPad->SetLogz();
    cQ2xResolution->SetGrid();
    cQ2xResolution->Update();
    cQ2xResolution->SaveAs(Q2xResolutionCanvasName);

    // Create canvas for momentum component plots
    TCanvas *cMomentum = new TCanvas("momentum_canvas", "Momentum Resolution", 3000, 1600);
    cMomentum->Divide(3, 2);
    cMomentum->cd(1);
    px_Hist->Draw("colz");
    gPad->SetLogz();
    cMomentum->cd(2);
    py_Hist->Draw("colz");
    gPad->SetLogz();
    cMomentum->cd(3);
    pz_Hist->Draw("colz");
    gPad->SetLogz();
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
    gPad->SetLogz();
    cEnergyThetaPhi->cd(2);
    theta_Hist->Draw("colz");
    gPad->SetLogz();
    cEnergyThetaPhi->cd(3);
    phi_Hist->Draw("colz");
    gPad->SetLogz();
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
    gPad->SetLogz();
    cResolutionVsMC->cd(2);
    E_res_vs_theta_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(3);
    E_res_vs_phi_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(4);
    theta_diff_vs_E_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(5);
    theta_diff_vs_theta_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(6);
    theta_diff_vs_phi_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(7);
    phi_diff_vs_E_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(8);
    phi_diff_vs_theta_Hist->Draw("colz");
    gPad->SetLogz();
    cResolutionVsMC->cd(9);
    phi_diff_vs_phi_Hist->Draw("colz");
    gPad->SetLogz();
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
    hPhi_vs_theta_mean->SetTitle("Mean Phi Offset vs theta MC; theta MC [rad]; Mean Phi Offset [deg]");
    hPhi_vs_theta_mean->SetMarkerStyle(20);
    hPhi_vs_theta_mean->SetMarkerColor(kBlue);
    hPhi_vs_theta_mean->SetMaximum(5); // Adjust maximum for better visibility
    hPhi_vs_theta_mean->SetMinimum(-5); // Adjust minimum for better visibility
    hPhi_vs_theta_mean->Draw();
    cResolutionGraphs->cd(6);
    hPhi_vs_theta_stddev->SetTitle("Std Dev Phi Offset vs theta MC; theta MC [rad]; Std Dev Phi Offset [deg]");
    hPhi_vs_theta_stddev->SetMarkerStyle(20);
    hPhi_vs_theta_stddev->SetMarkerColor(kRed);
    hPhi_vs_theta_stddev->SetMaximum(60); // Adjust maximum for better visibility
    hPhi_vs_theta_stddev->SetMinimum(0); // Adjust minimum for better visibility
    hPhi_vs_theta_stddev->Draw();
    cResolutionGraphs->SetGrid();
    cResolutionGraphs->Update();
    // Save the canvas as a PNG file
    cResolutionGraphs->SaveAs(resolutionGraphsCanvasName);

    TFile *f = new TFile(outFile,"RECREATE");
    cMomentum->Write();
    cEnergyThetaPhi->Write();
    cResolutionVsMC->Write();
    cResolutionGraphs->Write();
    cAcceptance1D->Write();
    cAcceptance2D->Write();
    cAcceptanceXQ2->Write();
    cQ2xResolution->Write();
    px_Hist->Write();
    py_Hist->Write();
    pz_Hist->Write();
    px_diff_Hist->Write();
    py_diff_Hist->Write();
    pz_res_Hist->Write();
    E_Hist->Write();
    theta_Hist->Write();
    phi_Hist->Write();
    E_res_Hist->Write();
    theta_diff_Hist->Write();
    phi_diff_Hist->Write();
    E_res_vs_E_Hist->Write();
    E_res_vs_theta_Hist->Write();
    E_res_vs_phi_Hist->Write();
    theta_diff_vs_E_Hist->Write();
    theta_diff_vs_theta_Hist->Write();
    theta_diff_vs_phi_Hist->Write();
    phi_diff_vs_E_Hist->Write();
    phi_diff_vs_theta_Hist->Write();
    phi_diff_vs_phi_Hist->Write();
    // Write Q2 and x resolution histograms
    Q2_Hist->Write();
    x_Hist->Write();
    log10Q2_comp_Hist->Write();
    log10x_comp_Hist->Write();
    Q2_diff_Hist->Write();
    x_diff_Hist->Write();
    log10Q2_diff_Hist->Write();
    log10x_diff_Hist->Write();
    Q2_diff_vs_Q2_Hist->Write();
    x_diff_vs_x_Hist->Write();
    log10Q2_diff_vs_log10Q2_Hist->Write();
    log10x_diff_vs_log10x_Hist->Write();
    // Write acceptance related histograms
    E_all_Hist->Write(); E_acc_Hist->Write(); hE_acceptance->Write();
    theta_all_Hist->Write(); theta_acc_Hist->Write(); hTheta_acceptance->Write();
    phi_all_Hist->Write(); phi_acc_Hist->Write(); hPhi_acceptance->Write();
    log10Q2_all_Hist->Write(); log10Q2_acc_Hist->Write(); hLog10Q2_acceptance->Write();
    W_all_Hist->Write(); W_acc_Hist->Write(); hW_acceptance->Write();
    log10Q2_vs_E_all_Hist->Write(); log10Q2_vs_E_acc_Hist->Write(); hLog10Q2_vs_E_acceptance->Write();
    log10Q2_vs_log10x_all_Hist->Write(); log10Q2_vs_log10x_acc_Hist->Write(); hLog10Q2_vs_log10x_acceptance->Write();

    hE_vs_E_mean->Write();
    hE_vs_E_stddev->Write();
    hTheta_vs_E_mean->Write();
    hTheta_vs_E_stddev->Write();
    hPhi_vs_theta_mean->Write();
    hPhi_vs_theta_stddev->Write();

    f->Close();


    
    // // Get mean and error on the mean of E, theta and phi resolutions
    // double mean_E_res = E_res_Hist->GetMean();
    // double mean_theta_res = theta_diff_Hist->GetMean();
    // double mean_phi_res = phi_diff_Hist->GetMean();
    // double mean_E_res_error = E_res_Hist->GetMeanError();
    // double mean_theta_res_error = theta_diff_Hist->GetMeanError();
    // double mean_phi_res_error = phi_diff_Hist->GetMeanError();

    // // Get standard deviation of E, theta and phi resolutions
    // double stddev_E_res = E_res_Hist->GetStdDev();
    // double stddev_theta_res = theta_diff_Hist->GetStdDev();
    // double stddev_phi_res = phi_diff_Hist->GetStdDev();
    // double stddev_E_res_error = E_res_Hist->GetStdDevError();
    // double stddev_theta_res_error = theta_diff_Hist->GetStdDevError();
    // double stddev_phi_res_error = phi_diff_Hist->GetStdDevError();

    // // Print the resolutions
    // std::cout << "Mean E offset: " << mean_E_res << " +/- " << mean_E_res_error << std::endl;
    // std::cout << "Mean theta offset: " << mean_theta_res << " +/- " << mean_theta_res_error << std::endl;
    // std::cout << "Mean phi offset: " << mean_phi_res << " +/- " << mean_phi_res_error << std::endl;
    // std::cout << "Standard deviation of E resolution: " << stddev_E_res << " +/- " << stddev_E_res_error << std::endl;
    // std::cout << "Standard deviation of theta resolution: " << stddev_theta_res << " +/- " << stddev_theta_res_error << std::endl;
    // std::cout << "Standard deviation of phi resolution: " << stddev_phi_res << " +/- " << stddev_phi_res_error << std::endl;

    // int pass = 0;

    // // Fail if mean is more than 20% of the standard deviation away from zero
    // if(std::abs(mean_E_res) > 0.2 * stddev_E_res) {
    //     std::cout << "Mean E offset is more than 20\% (" << 0.2 * stddev_E_res << ") of the standard deviation away from zero!" << std::endl;
    //     pass = 1;
    // }
    // if(std::abs(mean_theta_res) > 0.2 * stddev_theta_res) {
    //     std::cout << "Mean theta offset is more than 20\% (" << 0.2 * stddev_theta_res << ") of the standard deviation away from zero!" << std::endl;
    //     pass = 1;
    // }
    // if(std::abs(mean_phi_res) > 0.2 * stddev_phi_res) {
    //     std::cout << "Mean phi offset is more than 20\% (" << 0.2 * stddev_phi_res << ") of the standard deviation away from zero!" << std::endl;
    //     pass = 1;
    // }

    // // Resolution limits
    // double E_res_limit = 0.05; // 5% resolution
    // double theta_res_limit = 0.0001; // 1 mrad resolution
    // double phi_res_limit = 30; // 30 degrees resolution

    // // Fail if standard deviation is more than the limit
    // if(std::abs(stddev_E_res) > E_res_limit) {
    //     std::cout << "Standard deviation of E resolution is more than the limit of " << E_res_limit << "!" << std::endl;
    //     pass = 1;
    // }
    // if(std::abs(stddev_theta_res) > theta_res_limit) {
    //     std::cout << "Standard deviation of theta resolution is more than the limit of " << theta_res_limit << " radians!" << std::endl;
    //     pass = 1;
    // }
    // if(std::abs(stddev_phi_res) > phi_res_limit) {
    //     std::cout << "Standard deviation of phi resolution is more than the limit of " << phi_res_limit << " degrees!" << std::endl;
    //     pass = 1;
    // }

}

