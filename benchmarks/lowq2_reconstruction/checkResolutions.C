#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>
#include <fstream>

int checkResolutions(const TString inputFile="/home/simong/EIC/detector_benchmarks_anl/sim_output/beamline/acceptanceTestcurrent.edm4hep.root", const TString outputFile="test.json") {
    
    int fail = 0;

    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return 1; // Return error code
    }

    TH1F* E_res_Hist = (TH1F*)file->Get("E_res");
    TH1F* theta_diff_Hist = (TH1F*)file->Get("theta_diff");
    TH1F* phi_diff_Hist = (TH1F*)file->Get("phi_diff");


    double mean_E_res = E_res_Hist->GetMean();
    double mean_theta_res = theta_diff_Hist->GetMean();
    double mean_phi_res = phi_diff_Hist->GetMean();
    double mean_E_res_error = E_res_Hist->GetMeanError();
    double mean_theta_res_error = theta_diff_Hist->GetMeanError();
    double mean_phi_res_error = phi_diff_Hist->GetMeanError();

    // Get standard deviation of E, theta and phi resolutions
    double stddev_E_res = E_res_Hist->GetStdDev();
    double stddev_theta_res = theta_diff_Hist->GetStdDev();
    double stddev_phi_res = phi_diff_Hist->GetStdDev();
    double stddev_E_res_error = E_res_Hist->GetStdDevError();
    double stddev_theta_res_error = theta_diff_Hist->GetStdDevError();
    double stddev_phi_res_error = phi_diff_Hist->GetStdDevError();

    // Print the resolutions
    std::cout << "Mean E offset: " << mean_E_res << " +/- " << mean_E_res_error << std::endl;
    std::cout << "Mean theta offset: " << mean_theta_res << " +/- " << mean_theta_res_error << std::endl;
    std::cout << "Mean phi offset: " << mean_phi_res << " +/- " << mean_phi_res_error << std::endl;
    std::cout << "Standard deviation of E resolution: " << stddev_E_res << " +/- " << stddev_E_res_error << std::endl;
    std::cout << "Standard deviation of theta resolution: " << stddev_theta_res << " +/- " << stddev_theta_res_error << std::endl;
    std::cout << "Standard deviation of phi resolution: " << stddev_phi_res << " +/- " << stddev_phi_res_error << std::endl;

    // Fail if mean is more than 20% of the standard deviation away from zero
    if(std::abs(mean_E_res) > 0.2 * stddev_E_res) {
        std::cout << "Mean E offset is more than 20\% (" << 0.2 * stddev_E_res << ") of the standard deviation away from zero!" << std::endl;
        fail = 1;
    }
    if(std::abs(mean_theta_res) > 0.2 * stddev_theta_res) {
        std::cout << "Mean theta offset is more than 20\% (" << 0.2 * stddev_theta_res << ") of the standard deviation away from zero!" << std::endl;
        fail = 1;
    }
    if(std::abs(mean_phi_res) > 0.2 * stddev_phi_res) {
        std::cout << "Mean phi offset is more than 20\% (" << 0.2 * stddev_phi_res << ") of the standard deviation away from zero!" << std::endl;
        fail = 1;
    }

    // Resolution limits
    double E_res_limit = 0.05; // 5% resolution
    double theta_res_limit = 0.001; // 1 mrad resolution
    double phi_res_limit = 0; // 30 degrees resolution

    // Fail if standard deviation is more than the limit
    if(std::abs(stddev_E_res) > E_res_limit) {
        std::cout << "E resolution is more than the limit of " << E_res_limit << "!" << std::endl;
        fail = 1;
    }
    if(std::abs(stddev_theta_res) > theta_res_limit) {
        std::cout << "Theta resolution is more than the limit of " << theta_res_limit << " radians!" << std::endl;
        fail = 1;
    }
    if(std::abs(stddev_phi_res) > phi_res_limit) {
        std::cout << "Phi resolution is more than the limit of " << phi_res_limit << " degrees!" << std::endl;
        fail = 1;
    }

    // Create json output file containing the resolutions, errors and offsets
    std::ofstream jsonFile(outputFile);
    if (!jsonFile.is_open()) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return 1; // Return error code
    }
    jsonFile << "{\n";
    jsonFile << "  \"mean_E_res\": " << mean_E_res << ",\n";
    jsonFile << "  \"mean_theta_res\": " << mean_theta_res << ",\n";
    jsonFile << "  \"mean_phi_res\": " << mean_phi_res << ",\n";
    jsonFile << "  \"mean_E_res_error\": " << mean_E_res_error << ",\n";
    jsonFile << "  \"mean_theta_res_error\": " << mean_theta_res_error << ",\n";
    jsonFile << "  \"mean_phi_res_error\": " << mean_phi_res_error << ",\n";
    jsonFile << "  \"stddev_E_res\": " << stddev_E_res << ",\n";
    jsonFile << "  \"stddev_theta_res\": " << stddev_theta_res << ",\n";
    jsonFile << "  \"stddev_phi_res\": " << stddev_phi_res << ",\n";
    jsonFile << "  \"stddev_E_res_error\": " << stddev_E_res_error << ",\n";
    jsonFile << "  \"stddev_theta_res_error\": " << stddev_theta_res_error << ",\n";
    jsonFile << "  \"stddev_phi_res_error\": " << stddev_phi_res_error << ",\n";
    jsonFile << "  \"fail\": " << fail << "\n";
    jsonFile << "}\n";
    jsonFile.close();

    return fail;

}