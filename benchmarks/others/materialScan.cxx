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

void materialScan(
    double etamin = -2.0,
    double etamax = +1.0,
    double etastep = 0.2,
    double rmax = 100.0,
    double phi = 0.0)
{
  // check inputs
  if (etamin > etamax || rmax <= 0.0) {
    std::cout << "Error: ordered eta range required" << std::endl;
    return -1;
  }

  // get material scans
  size_t total{0};
  std::vector<dd4hep::rec::MaterialVec> scan;
  double x0{0}, y0{0}, z0{0};
  for (double eta = etamin; eta <= etamax + 0.5*etastep; eta += etastep) {
    double theta = 2.0 * atan(exp(-eta));
    double x = rmax * cos(phi) * sin(theta);
    double y = rmax * sin(phi) * sin(theta);
    double z = rmax * cos(theta);
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
  THStack hs("hs",Form("Material Scan (max distance %.0f cm)",rmax));
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
  hs.SetMaximum(100.);
  cs.SaveAs("materialScan.png");
  cs.SaveAs("materialScan.pdf");
}
