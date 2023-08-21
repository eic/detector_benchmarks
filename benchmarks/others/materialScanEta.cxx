#import <iostream>
#import <vector>

#include <DDRec/Material.h>
#include <DDRec/MaterialScan.h>

#include <TH1F.h>
#include <THStack.h>

Color_t color(const Material& m) {
  if      (m.name() == std::string("Silicon"))         return kGray;
  else if (m.name() == std::string("Aluminum"))        return kAzure;
  else if (m.name() == std::string("CarbonFiber"))     return kGray;
  else if (m.name() == std::string("Beryllium"))       return kGreen;
  else if (m.name() == std::string("Gold"))            return kYellow;
  else if (m.name() == std::string("Mylar"))           return kGreen;
  else if (m.name() == std::string("Kapton"))          return kGreen;
  else if (m.name() == std::string("Copper"))          return kGreen;
  else if (m.name() == std::string("C2F6_DRICH"))      return kOrange;
  else if (m.name() == std::string("Ar10CO2"))         return kOrange;
  else if (m.name() == std::string("Aerogel"))         return kPink;
  else if (m.name() == std::string("AerogelOptical"))  return kPink;
  else if (m.name() == std::string("Aerogel_DRICH"))   return kPink;
  else if (m.name() == std::string("Lead"))            return kBlack;
  else if (m.name() == std::string("Steel235"))        return kGray+2;
  else if (m.name() == std::string("TungstenDens24"))  return kBlack;
  else if (m.name() == std::string("Polystyrene"))     return kGray;
  else if (m.name() == std::string("PolystyreneFoam")) return kGray;
  else if (m.name() == std::string("Epoxy"))           return kGray;
  else if (m.name() == std::string("PlasticScint"))    return kGray;
  else if (m.name() == std::string("AcrylicOptical"))  return kGray;
  else if (m.name() == std::string("Acrylic_DRICH"))   return kGray;
  else if (m.name() == std::string("Quartz"))          return kViolet;
  else if (m.name() == std::string("Air"))             return kBlue;
  else if (m.name() == std::string("AirOptical"))      return kBlue;
  else if (m.name() == std::string("Vacuum"))          return kWhite;
  else {
    std::cout << "Unknown material: " << m.name() << std::endl;
    return kRed;
  }
}

void materialScanEta(
    double etamin = -2.0, // minimum eta
    double etamax = +1.0, // maximum eta
    double etastep = 0.2, // steps in eta
    double phi = 0.0, // phi angle for material scan
    double rhomax = 10000.0, // maximum distance from z axis
    double znmax = 10000.0, // maximum negative endcap z plane (positive number)
    double zpmax = 10000.0 // maximum positive endcap z plane (positive number)
) {
  // check inputs
  if (etamin > etamax) {
    std::cout << "Error: ordered eta range required" << std::endl;
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
  std::vector<dd4hep::rec::MaterialVec> scan;
  double x0{0}, y0{0}, z0{0};
  for (double eta = etamin; eta <= etamax + 0.5*etastep; eta += etastep) {
    double theta = 2.0 * (atan(1) - atan(exp(-eta))); // |theta| < 90 deg, cos(theta) > 0, sin(theta) can be 0
    double r = min((theta < 0? -znmax: zpmax) / sin(theta), rhomax / cos(theta)); // theta == 0 results in min(+inf, rhomax)
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
    size_t eta0 = static_cast<unsigned int>(std::rint(abs((etamax-etamin)/etastep/2)));
    for (size_t i = 0, j = eta0;
         i < scan.size();
         ++i, j += (2*eta0 < scan.size()? -1: +1) * (i <= 2*min(eta0,scan.size()-eta0-1)? (i%2==0? -i: +i): -1)) {
      if (scan.at(j).size() == 0) continue;
      // define layer
      auto layer_anchor = scan.at(j).at(0);
      histograms.emplace_back(Form("h%zu",layer++),layer_anchor.first.name(),scan.size(),etamin,etamax);
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
  THStack hs("hs",Form("Material Scan (rho < %.0f cm, -%.0f cm < z < %.0f cm)", rhomax, znmax, zpmax));
  for (auto& h: histograms) {
    hs.Add(&h);
  }
  TCanvas cs("cs","Material Scan",1920,1080);
  auto pad = cs.cd();
  pad->SetLogy();
  hs.Draw();
  hs.GetXaxis()->SetTitle("eta");
  hs.GetYaxis()->SetTitle("Fraction X0");
  hs.SetMinimum(2.5e-3);
  cs.SaveAs("materialScanEta.png");
  cs.SaveAs("materialScanEta.pdf");
}
