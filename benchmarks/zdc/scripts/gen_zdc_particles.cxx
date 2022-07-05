//////////////////////////////////////////////////////////////
// ZDC detector
// Single Neutron dataset
// J.KIM 05/07/2021
//////////////////////////////////////////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include <random>

#include "TMath.h"
#include "TRandom.h"

#include <Math/Vector3D.h>
#include <Math/RotationY.h>

using namespace HepMC3;

void gen_zdc_particles(int n_events = 1e6, const std::string& particle = "neutron", double p_start = 125.0, double p_end = 145.0, const std::string& out_fname = "./data/zdc_neutrons.hepmc") {
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom* r1 = new TRandom();

  // Set crossing angle [rad]
  double crossing_angle = -0.025;

  // Constraining the solid angle, but larger than that subtended by the
  // detector
  // https://indico.bnl.gov/event/7449/contributions/35966/attachments/27177/41430/EIC-DWG-Calo-03192020.pdf
  // See a figure on slide 26
  double cos_theta_min = std::cos(0.0);
  double cos_theta_max = std::cos(0.005);

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 111 - pi0
    // pdgid 2212 - proton
    // pdgid 2112 - neutron
    GenParticlePtr p1 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum
    Double_t p_mag    = r1->Uniform(p_start, p_end);
    Double_t phi      = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t costheta = r1->Uniform(cos_theta_min, cos_theta_max);
    Double_t theta    = std::acos(costheta);
    Double_t p0_x     = p_mag * std::cos(phi) * std::sin(theta);
    Double_t p0_y     = p_mag * std::sin(phi) * std::sin(theta);
    Double_t p0_z     = p_mag * std::cos(theta);

    // Rotate the vector in Y by crossing angle when particles are being generated
    ROOT::Math::XYZVector p0{p0_x, p0_y, p0_z};
    ROOT::Math::RotationY r(crossing_angle);
    auto p = r * p0;
    auto px = p.X();
    auto py = p.Y();
    auto pz = p.Z();

    // type 1 is final state
    // pdgid 11 - electron 0.510 MeV/c^2
    // pdgid 22 - photon massless
    // pdgid 2112 - neutron 939.565 MeV/c^2
    GenParticlePtr p3;
    if (particle == "neutron") {
      p3 = std::make_shared<GenParticle>(FourVector(px, py, pz, std::hypot(p_mag, 0.939565)), 2112, 1);
    } else if (particle == "photon") {
      p3 = std::make_shared<GenParticle>(FourVector(px, py, pz, p_mag), 22, 1);
    } else {
      std::cout << "Particle type " << particle << " not recognized!" << std::endl;
      exit(-1);
    }

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
