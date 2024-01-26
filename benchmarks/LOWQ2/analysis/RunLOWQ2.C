#include "LOWQ2Benchmarks.h"

void RunLOWQ2(){
    
    std::string compactName = "/home/simong/EIC/epic/epic_18x275.xml";
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(compactName);

    double luminosity = 1e34; // [cm^-2 s^-1]
    double eBeamEnergy = 18.0; // [GeV]
    double pBeamEnergy = 275.0; // [GeV]
    double eventCrossSectionQR = 0.0551; // [mb]
    double eventCrossSectionBrems = 1.0; // [mb]
    double bunchSpacing = 10.15*1e-9; // [s]

    double eventRateQR    = luminosity * eventCrossSectionQR * 1e-27; // [Hz]
    double eventRateBrems = luminosity * eventCrossSectionBrems * 1e-27; // [Hz]
    double bunchRate      = 1.0 / bunchSpacing; // [Hz]

    LOWQ2Benchmarks("/scratch/EIC/ReconOut/tempEventsQR2.root","LOWQ2QR.root",detector,eventRateQR);
    LOWQ2Benchmarks("/scratch/EIC/ReconOut/qr_18x275_ab/qr_18x275_ab0_recon.edm4hep.root","LOWQ2QROld.root",detector,eventRateQR);
    LOWQ2Benchmarks("/scratch/EIC/G4out/qr_18x275_new.edm4hep*.root","LOWQ2QRRates.root",detector,eventRateQR);
    LOWQ2Benchmarks("/scratch/EIC/G4out/brems_18x275_new.edm4hep*.root","LOWQ2BremsRates.root",detector,bunchRate);


}