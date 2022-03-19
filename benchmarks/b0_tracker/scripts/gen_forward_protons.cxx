#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
#include <random>
#include <cmath>
#include <math.h>

#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/RotationY.h"
#include "TRandom.h"

#include "common_bench/particles.h"

using namespace HepMC3;

/** Generate electrons in the central region.
 *  This is for testing detectors in the "barrel" region.
 */
void gen_forward_protons(int n_events = 100, 
                     const char* out_fname = "forward_protons.hepmc")
{
  double cos_theta_min = std::cos(0.5*(M_PI/180.0));
  double cos_theta_max = std::cos(0.0*(M_PI/180.0));

  const double M_p = common_bench::particleMap.at(2212).mass;

  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom *r1 = new TRandom();

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 111 - pi0
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 2212, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, M_p), 2212, 4);

    // Define momentum
    Double_t p     = r1->Uniform(200.0, 275.0);
    Double_t phi   = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t costh = r1->Uniform(cos_theta_min, cos_theta_max);
    Double_t th    = std::acos(costh);
    Double_t px    = p * std::cos(phi) * std::sin(th);
    Double_t py    = p * std::sin(phi) * std::sin(th);
    Double_t pz    = p * std::cos(th);


    ROOT::Math::XYZVector p0 = {px,py,pz};

    //ROOT::Math::Rotation3D r = (-0.025);
    ROOT::Math::RotationY r(-0.025);
    auto p_rot = r*p0;


    // Generates random vectors, uniformly distributed over the surface of a
    // sphere of given radius, in this case momentum.
    // r1->Sphere(px, py, pz, p);

    //std::cout << std::sqrt(px*px + py*py + pz*pz) - p << " is zero? \n";

    // type 1 is final state
    // pdgid 11 - electron 0.510 MeV/c^2
    GenParticlePtr p3 = std::make_shared<GenParticle>(
        FourVector(p_rot.x(), p_rot.y(), p_rot.z(), sqrt(p * p + (M_p * M_p))), 2212, 1);

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);

    v1->add_particle_out(p3);
    evt.add_vertex(v1);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 10000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
