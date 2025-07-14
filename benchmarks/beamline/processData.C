#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHit.h"
#include <iostream>

void processData(const TString inputFile="/home/simong/EIC/detector_benchmarks_anl/sim_output/beamline/acceptanceTestcurrent.edm4hep.root", const TString outputFile="test.root", const int desired_cellID = 66757) {
    // Define the branches to read
    std::vector<std::string> branches = {
        "BackwardsBeamlineHits.cellID",
        "BackwardsBeamlineHits.position.x",
        "BackwardsBeamlineHits.position.y",
        "BackwardsBeamlineHits.position.z",
        "BackwardsBeamlineHits.momentum.x",
        "BackwardsBeamlineHits.momentum.y",
        "BackwardsBeamlineHits.momentum.z",
        "MCParticles.momentum.x",
        "MCParticles.momentum.y",
        "MCParticles.momentum.z"
    };

    float momentum_tolerance = 0.0001; // Define the momentum tolerance for filtering

    // Create a ROOT DataFrame to read the input files
    ROOT::RDataFrame df("events", inputFile);

    //Filter on events with only a single MCParticle
    auto filterDF = df.Filter("MCParticles.size()==1")
                    .Define("desiredHits", "BackwardsBeamlineHits[BackwardsBeamlineHits.cellID=="+std::to_string(desired_cellID)+"]")
                    .Filter("desiredHits.size()==1")
                    .Define("features", "std::array<float, 6>{static_cast<float>(desiredHits[0].position.x), static_cast<float>(desiredHits[0].position.y), static_cast<float>(desiredHits[0].position.z), static_cast<float>(desiredHits[0].momentum.x), static_cast<float>(desiredHits[0].momentum.y), static_cast<float>(desiredHits[0].momentum.z)}")
                    .Define("targets", "std::array<float, 3>{static_cast<float>(MCParticles[0].momentum.x), static_cast<float>(MCParticles[0].momentum.y), static_cast<float>(MCParticles[0].momentum.z)}")
                    .Define("features_momentum", [](const std::array<float, 6>& features) {
                        return std::sqrt(features[3] * features[3] + features[4] * features[4] + features[5] * features[5]);
                    }, {"features"})
                    .Define("targets_momentum", [](const std::array<float, 3>& targets) {
                        return std::sqrt(targets[0] * targets[0] + targets[1] * targets[1] + targets[2] * targets[2]);
                    }, {"targets"})
                    .Filter([momentum_tolerance](float features_momentum, float targets_momentum) {
                        float relative_difference = std::abs(features_momentum - targets_momentum) / targets_momentum;
                        return relative_difference <= momentum_tolerance;
                    }, {"features_momentum", "targets_momentum"});;

    // Save the filtered data to a new ROOT file
    filterDF.Snapshot("events", outputFile, {"features","targets"});

    std::cout << "Filtered data saved to " << outputFile << std::endl;
}