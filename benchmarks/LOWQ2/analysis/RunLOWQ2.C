#include "LOWQ2Benchmarks.h"
#include <cstdlib>

void RunLOWQ2(  std::string inputFileName  = "Brems_input.root", 
                std::string outputFileName = "plots/LOWQ2QRRecon3.root",
                std::string compactName    = "/opt/detector/epic-nightly/share/epic/epic.xml",
                bool   inputIsTimeBased    = false,      // true if the event sample is time-based, false if it is event-based
                double timeWindow          = 10.15*1e-9, //[s]
                double eventCrossSection   = 0.0551,     // [mb]
                double luminosity          = 1e34,       // [cm^-2 s^-1]
            ) {
    
    //Set implicit multi-threading
    ROOT::EnableImplicitMT();
    
    // Output script running conditions
    std::cout << "Running LOWQ2 benchmarks with the following parameters:" << std::endl;
    std::cout << "  - input file: " << inputFileName << std::endl;
    std::cout << "  - output file: " << outputFileName << std::endl;
    std::cout << "  - xml file: " << compactName << std::endl;
    std::cout << "  - input is time-based: " << inputIsTimeBased << std::endl;
    std::cout << "  - time window: " << timeWindow << " s" << std::endl;
    std::cout << "  - event cross section: " << eventCrossSection << " mb" << std::endl;
    std::cout << "  - luminosity: " << luminosity << " cm^-2 s^-1" << std::endl;

    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(compactName);

    // eventCrossSectionBrems = 171.3; [mb]
    // eventCrossSectionQR    = 0.0551; [mb]

    double eventRate      = luminosity * eventCrossSection * 1e-27; // [Hz]

    if(inputIsTimeBased){
        eventRate = 1.0 / timeWindow; // [Hz]
    }
   
    LOWQ2Benchmarks(inputFileName,outputFileName,detector,eventRate);

}