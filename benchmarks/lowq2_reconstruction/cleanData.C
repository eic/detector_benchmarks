#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHit.h"
#include <iostream>

void cleanData(const TString inputFile="/home/simong/EIC/detector_benchmarks_anl/sim_output/beamline/acceptanceTestcurrent.edm4hep.root", const TString outputFile="test.root", const double BeamEnergy=18.0, const int desired_cellID = 66757, const bool appendTruth = true) {
    
    float momentum_tolerance = 0.01; // Define the momentum tolerance for filtering

    // Create a ROOT DataFrame to read the input files
    ROOT::RDataFrame df("events", inputFile);

    //Filter on events with only a single MCParticle which hasn't scattered as it enters the drift volume
    auto filterDF = df.Define("SimParticles", "MCParticles[MCParticles.generatorStatus==1]")
                    .Filter("SimParticles.size()==1")
                    .Define("beamlineHit", "BackwardsBeamlineHits[BackwardsBeamlineHits.cellID=="+std::to_string(desired_cellID)+"]")
                    .Filter("beamlineHit.size()==1")
                    .Define("features", [BeamEnergy](const ROOT::VecOps::RVec<edm4hep::SimTrackerHitData>& hit) {
                        return std::array<float, 6>{
                            static_cast<float>(hit[0].position.x),
                            static_cast<float>(hit[0].position.y),
                            static_cast<float>(hit[0].position.z),
                            static_cast<float>(hit[0].momentum.x / BeamEnergy),
                            static_cast<float>(hit[0].momentum.y / BeamEnergy),
                            static_cast<float>(hit[0].momentum.z / BeamEnergy)
                        };
                    }, {"beamlineHit"})
                    .Define("targets", [BeamEnergy](const ROOT::VecOps::RVec<edm4hep::MCParticleData>& mcps) {
                            return std::array<float, 3>{
                                static_cast<float>(mcps[0].momentum.x / BeamEnergy),
                                static_cast<float>(mcps[0].momentum.y / BeamEnergy),
                                static_cast<float>(mcps[0].momentum.z / BeamEnergy)
                            };
                        }, {"SimParticles"})
                    .Define("features_momentum", [](const std::array<float, 6>& features) {
                        return std::sqrt(features[3] * features[3] + features[4] * features[4] + features[5] * features[5]);
                    }, {"features"})
                    .Define("targets_momentum", [](const std::array<float, 3>& targets) {
                        return std::sqrt(targets[0] * targets[0] + targets[1] * targets[1] + targets[2] * targets[2]);
                    }, {"targets"})
                    .Filter([momentum_tolerance](float features_momentum, float targets_momentum) {
                        float relative_difference = std::abs(features_momentum - targets_momentum) / targets_momentum;
                        return relative_difference <= momentum_tolerance;
                    }, {"features_momentum", "targets_momentum"});

    auto taggerDF = filterDF.Filter("_TaggerTrackerFeatureTensor_shape[0]==1");

    // Save the filtered data to a new ROOT file
    taggerDF.Snapshot("events", outputFile, {"_TaggerTrackerFeatureTensor_floatData","_TaggerTrackerFeatureTensor_shape","_TaggerTrackerTargetTensor_floatData","features_momentum","targets_momentum"});

    // Print the size of the original DataFrame
    // std::cout << "Original DataFrame size: " << df.Count().GetValue() << std::endl;
    // // Print the size of the filtered DataFrame
    // std::cout << "Filtered DataFrame size: " << filterDF.Count().GetValue() << std::endl;
    
    // std::cout << "Tagger filtered DataFrame size" << taggerDF.Count().GetValue() << std::endl;

    // std::cout << "Filtered data saved to " << outputFile << std::endl;

    // If appendTruth is true, add the truth information
    if (appendTruth) {
        // Open the output file in update mode
        ROOT::RDF::RSnapshotOptions opts;
        opts.fMode = "update";
        auto aliasDF = filterDF.Redefine("_TaggerTrackerFeatureTensor_floatData", "features")
                                .Redefine("_TaggerTrackerTargetTensor_floatData", "targets");
        // filterDF = Concatenate({aliasDF, filterDF});
        aliasDF.Snapshot("truthevents", outputFile, {"_TaggerTrackerFeatureTensor_floatData", "_TaggerTrackerFeatureTensor_shape", "_TaggerTrackerTargetTensor_floatData"}, opts);
        std::cout << "Truth information appended to " << outputFile << std::endl;

        // std::cout << "Total events after appending truth: " << aliasDF.Count().GetValue() << std::endl;
    }




}