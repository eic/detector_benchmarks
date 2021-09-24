// material scan of a cylinder volume for athena
// the scan starts at 0 for now
// based on https://dd4hep.web.cern.ch/dd4hep/reference/MaterialScan_8cpp_source.html
// dd4hep also provides command line utilities materialBudget and materialScan
//        
// Shujie Li, Aug 2021

R__LOAD_LIBRARY(libDDCore.so)
R__LOAD_LIBRARY(libDDG4.so)
R__LOAD_LIBRARY(libDDG4IO.so)
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Detector.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/CellIDPositionConverter.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TGeoMedium.h"
#include "TGeoManager.h"
#include "DDRec/MaterialScan.h"
#include "DDRec/MaterialManager.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Printout.h"
#include "fmt/core.h"

#include <iostream>
#include <fstream>

using namespace dd4hep;
using namespace dd4hep::rec;

  // for silicon tracker, the outer barrel rmax ~45, outer endcap zmax~121
void test_matscan(const char* compact = "athena.xml",double zmax=130, double rmax = 61){

  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  detector.fromCompact(compact);
  MaterialScan matscan(detector); 
  TString detname = "all";
  fmt::print("\n");
  fmt::print("All detector subsystem names:\n");
  for(const auto&  d : detector.detectors() ) {
    fmt::print("  {}\n", d.first);
  }
  // to do: material scan of given dtectors.

  // TString det_list[14]={"TrackerBarrel_Inner","TrackerBarrel_Outer","TrackerEndcapN_Inner","TrackerEndcapN_Outer","TrackerEndcapP_Inner","TrackerEndcapP_Outer","TrackerSubAssembly_Inner","TrackerSubAssembly_Outer","VertexBarrel","VertexBarrelSubAssembly","VertexEndcapN","VertexEndcapP","VertexEndcapSubAssembly","cb_DIRC"};
  // for (int dd=0;dd<14;dd++){
  // TString detname = det_list[dd];
  // matscan.setDetector(detname.Data());      

  const char* fmt1 = "%8.3f %8.3f %8.3f %3d %-20s %3.0f  %8.4f %11.4f %11.4f %11.4f %11.4f %8.3f %8.3f %8.3f\n";

  // the beam vacuum is missing if starting from origin and ending at negative z. need to check later
  double x0=0,y0=0,z0=0.001,x1,y1,z1;  // cm
  double epsilon=1e-4; // (mm) default 1e-4: Materials with a thickness smaller than epsilon (default 1e-4=1mu

  const double DegToRad=3.141592653589793/180.0;
  double phi,dphi=0.5;    // Degree
  double dz=0.5, dr=0.5;    // cm

  vector<double> lx,ly,lz;
  // prepare coord. for barrel
  for(z1=-zmax;z1<=zmax;z1+=dz){
    for(phi=0.;phi<=360;phi+=dphi){
      x1 = cos(phi*DegToRad)*rmax;
      y1 = sin(phi*DegToRad)*rmax;
      lx.push_back(x1);
      ly.push_back(y1);
      lz.push_back(z1);
    }
  }
  // prepare coord. for endcaps
  for (double r=1; r<=rmax; r+=dr){
    for(phi=0.;phi<=360;phi+=dphi){
      x1 = cos(phi*DegToRad)*r;
      y1 = sin(phi*DegToRad)*r;
      lx.push_back(x1);
      ly.push_back(y1);
      lz.push_back(zmax);
      lx.push_back(x1);
      ly.push_back(y1);
      lz.push_back(-zmax);
    }
  }

  // loop over coord ponits for material scan
  long nn = lx.size();
  TString fname = Form("%s_z%g_r%g.dat",detname.Data(),zmax,rmax);
  FILE * pFile;
  pFile = fopen (fname,"w");

  for(int ii=0;ii<nn;ii++){
    x1=lx[ii]; y1=ly[ii]; z1=lz[ii];
    if(ii%5000==0)
      cout<<ii<<"/"<<nn<<":"<<x1<<" "<<y1<<" "<<z1<<endl;
    
    // if ((x1*x1+y1*y1)>(60*60))
      // continue;
    Vector3D p0(x0, y0, z0), p1(x1, y1, z1);
    Vector3D end, direction;
    direction = (p1-p0).unit();
  	const auto& placements = matscan.scan(x0,y0,z0,x1,y1,z1,epsilon); 
  	// matscan.print(x0,y0,z0,x1,y1,z1,epsilon);

    double sum_x0 = 0;
    double sum_lambda = 0;
    double path_length = 0, total_length = 0;


  	for (unsigned i=0;i<placements.size();i++){

      TGeoMaterial* mat=placements[i].first->GetMaterial();
      double length = placements[i].second;
      double nx0     = length / mat->GetRadLen();
      double nLambda = length / mat->GetIntLen();
      sum_x0        += nx0;
      sum_lambda    += nLambda;
      path_length   += length;
      total_length  += length;
      end = p0 + total_length * direction;


      fprintf(pFile, fmt1,x1, y1, z1, i+1, mat->GetName(), mat->GetZ(),
                mat->GetDensity(), mat->GetRadLen(), 
                length, path_length, sum_x0,end[0], end[1], end[2]);



  	} 
    // cout<<detname<<"  "<<x1<<","<<y1<<","<<z1<<endl;
    // cout<<x1<<","<<y1<<","<<z1<<": "<<placements.size()<<"  "<<sum_x0<<"  "<<total_length<<endl;
  }	
  fclose (pFile);
}

