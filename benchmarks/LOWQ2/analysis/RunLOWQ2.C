#include "LOWQ2Benchmarks.h"
#include <cstdlib>

void RunLOWQ2(std::string inputFileName="Brems_input.root", std::string outputFileName="plots/LOWQ2BremsRecon3.root",
                double eventCrossSection=0, std::string compactName="/opt/detector/epic-nightly/share/epic/epic.xml") {
    
    //Set implicit multi-threading
    ROOT::EnableImplicitMT();
    
    // Output script running conditions
    std::cout << "Running LOWQ2 benchmarks with the following parameters:" << std::endl;
    std::cout << "  - xml file: " << compactName << std::endl;

    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(compactName);

    double luminosity = 1e34; // [cm^-2 s^-1]
    double eBeamEnergy = 18.0; // [GeV]
    double pBeamEnergy = 275.0; // [GeV]
    double eventCrossSectionQR = 0.0551; // [mb]
    double eventCrossSectionBrems = 171.3; // [mb]
    double bunchSpacing = 10.15*1e-9; // [s]

    double eventRateQR    = luminosity * eventCrossSectionQR * 1e-27; // [Hz]
    double eventRateBrems = luminosity * eventCrossSectionBrems * 1e-27; // [Hz]
    double bunchRate      = 1.0 / bunchSpacing; // [Hz]

    double eventRate      = luminosity * eventCrossSection * 1e-27; // [Hz]

    //LOWQ2Benchmarks("/scratch/EIC/ReconOut/QR_new.root","plots/LOWQ2QRRecon2.root",detector,eventRateQR);
    LOWQ2Benchmarks("/scratch/EIC/ReconOut/Brems_new.root","plots/LOWQ2BremsRecon2.root",detector,eventRateBrems);
    LOWQ2Benchmarks("/scratch/EIC/ReconOut/Brems_new.root","plots/LOWQ2BremsRecon3.root",detector,bunchRate);
    //LOWQ2Benchmarks(inputFileName,outputFileName,detector,eventRate);

}