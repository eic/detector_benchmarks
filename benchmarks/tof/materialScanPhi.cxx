#import <iostream>
#import <vector>

#include <DDRec/Material.h>
#include <DDRec/MaterialScan.h>

#include <TH1F.h>
#include <THStack.h>

Color_t color(const Material& m) {
  if      (m.name() == std::string("Silicon"))         return kBlue;
  else if (m.name() == std::string("Kapton"))          return kOrange;
  else if (m.name() == std::string("Aluminum"))        return kRed;
  else if (m.name() == std::string("NOVEC7200"))       return kRed;
  else if (m.name() == std::string("CarbonFoam"))      return kGray;
  else if (m.name() == std::string("CFRPMix"))         return kGray;
  else if (m.name() == std::string("CFRPMix2"))        return kGray;
  else if (m.name() == std::string("Air"))             return kWhite;
  else if (m.name() == std::string("Vacuum"))          return kWhite;
  else {
    std::cout << "Unknown material: " << m.name() << std::endl;
    return kRed;
  }
}

void materialScanPhi(
    double phimin = -3.15, // minimum eta
    double phimax = +3.15, // maximum eta
    double phistep = 0.1, // steps in eta
    double rmin = 0.00, // minium radius to scan from
    double rmax = 100.0, // maximum radius to scan to
    double eta = 0.0, // eta for material scan
    double rhomin = 0.0, // minimum distance from z axis
    double rhomax = 10000.0, // maximum distance from z axis
    double znmax = 10000.0, // maximum negative endcap z plane (positive number)
    double zpmax = 10000.0 // maximum positive endcap z plane (positive number)
) {
  // check inputs
  if (phimin > phimax || rmax <= 0.0) {
    std::cout << "Error: ordered phi range required" << std::endl;
    return -1;
  }

  // get material scans
  size_t total{0};
  std::vector<dd4hep::rec::MaterialVec> scan;
  double x0{0}, y0{0}, z0{0};
  for (double phi = phimin; phi <= phimax + 0.5*phistep; phi += phistep) {

    double theta = 2.0 * (atan(1) - atan(exp(-eta)));
    double r0 = max(rmin, rhomin / cos(theta));
    double x0 = r0 * cos(theta) * cos(phi);
    double y0 = r0 * cos(theta) * sin(phi);
    double z0 = r0 * sin(theta);

    double r = min((theta > 0? zpmax: -znmax) / sin(theta), min(rmax, rhomax / cos(theta)));
    double x = r * cos(theta) * cos(phi);
    double y = r * cos(theta) * sin(phi);
    double z = r * sin(theta);

    scan.emplace_back(gMaterialScan->scan(x0,y0,z0,x,y,z));
    total += scan.back().size();
  }

  // start creating histograms for stacking:
  // - start with first material layer at central eta bin
  // - pop front layers first positive
  size_t layer = 0;
  std::vector<TH1F> histograms;
  while (total > 0) {
    // find next layer, starting from center bin outwards
    size_t phi0 = static_cast<unsigned int>(std::rint(abs((phimax-phimin)/phistep/2)));
    for (size_t i = 0, j = phi0;
         i < scan.size();
         ++i, j += (2*phi0 < scan.size()? -1: +1) * (i <= 2*min(phi0,scan.size()-phi0-1)? (i%2==0? -i: +i): -1)) {
      if (scan.at(j).size() == 0) continue;
      // define layer
      auto layer_anchor = scan.at(j).at(0);
      histograms.emplace_back(Form("h%zu",layer++),layer_anchor.first.name(),scan.size(),phimin,phimax);
      //histograms.back().SetLineColor(color(layer_anchor.first));
      histograms.back().SetFillColor(color(layer_anchor.first));
      // add all bins to this layer
      for (size_t bin = 0; bin < scan.size(); bin++) {
        double X0{0};
        for (auto& mat_len: scan.at(bin)) {
          if (mat_len.first.name() == layer_anchor.first.name()) {
            X0 +=  mat_len.second / mat_len.first.radLength();
            scan.at(bin).erase(scan.at(bin).begin()); total--;
          } else {
            break;
          }
        }
        histograms.back().SetBinContent(bin+1, X0); // bins start at 1
      }
    }
  }
  std::cout << histograms.size() << " histograms created" << std::endl;

  // plot histograms as stack
  THStack hs("hs",Form("Material Scan (%.0f cm < rho < %.0f cm, -%.0f cm < z < %.0f cm)", rhomin, rhomax, znmax, zpmax));
  for (auto& h: histograms) {
    hs.Add(&h);
  }
  TCanvas cs("cs","Material Scan",1920,1080);
  auto pad = cs.cd();
  pad->SetLogy();
  hs.Draw();
  hs.GetXaxis()->SetTitle("phi");
  hs.GetYaxis()->SetTitle("Fraction X0");
  hs.SetMinimum(0.5e-3);
  hs.SetMaximum(4e-2);
  cs.SaveAs("materialScanPhi.png");
  cs.SaveAs("materialScanPhi.pdf");

}
