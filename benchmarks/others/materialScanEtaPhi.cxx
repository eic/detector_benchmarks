#import <iostream>
#import <vector>

#include <DDRec/Material.h>
#include <DDRec/MaterialScan.h>

#include <TH1F.h>
#include <THStack.h>

void materialScanEtaPhi(
    double etamin = -2.0, // minimum eta
    double etamax = +1.0, // maximum eta
    double etastep = 0.2, // steps in eta
    double phimin = 0.0, // minimum phi
    double phimax = 2.0 * M_PI, // maximum phi
    double phistep = M_PI / 2, // steps in phi
    double rhomax = 10000.0, // maximum distance from z axis
    double znmax = 10000.0, // maximum negative endcap z plane (positive number)
    double zpmax = 10000.0 // maximum positive endcap z plane (positive number)
) {
  // check inputs
  if (etamin > etamax) {
    std::cout << "Error: ordered eta range required" << std::endl;
    return;
  }
  if (phimin > phimax) {
    std::cout << "Error: ordered phi range required" << std::endl;
    return;
  }
  if (rhomax <= 0.0) {
    std::cout << "Error: positive rhomax required" << std::endl;
    return;
  }
  if (znmax <= 0.0) {
    std::cout << "Error: positive znmax required" << std::endl;
    return;
  }
  if (zpmax <= 0.0) {
    std::cout << "Error: positive zpmax required" << std::endl;
    return;
  }

  // get material scans
  size_t total{0};
  double x0{0}, y0{0}, z0{0};
  TH2F h2("h2","Material Scan Eta vs Phi",
    (etamax-etamin)/etastep+1,etamin,etamax,
    (phimax-phimin)/phistep+1,phimin,phimax
  );
  for (double eta = etamin; eta <= etamax + 0.5*etastep; eta += etastep) {
    for (double phi = phimin; phi <= phimax + 0.5*phistep; phi += phistep) {
      double theta = 2.0 * (atan(1) - atan(exp(-eta))); // |theta| < 90 deg, cos(theta) > 0, sin(theta) can be 0 
      double r = min((theta < 0? -znmax: zpmax) / sin(theta), rhomax / cos(theta)); // theta == 0 results in min(+inf, rhomax)
      double x = r * cos(theta) * cos(phi);
      double y = r * cos(theta) * sin(phi);
      double z = r * sin(theta); 
      auto scan = gMaterialScan->scan(x0,y0,z0,x,y,z);
      double X0{0};
      for (auto& mat_len: scan)
        X0 +=  mat_len.second / mat_len.first.radLength();
      h2.Fill(eta,phi,X0);
    }
  }

  // plot histograms as stack
  TCanvas cs("cs","Material Scan Eta vs Phi",1920,1080);
  auto pad = cs.cd();
  cs.SetLogz();
  h2.Draw("colz");
  h2.GetXaxis()->SetTitle("eta");
  h2.GetYaxis()->SetTitle("phi");
  cs.SaveAs("materialScanEtaPhi.png");
  cs.SaveAs("materialScanEtaPhi.pdf");
}
