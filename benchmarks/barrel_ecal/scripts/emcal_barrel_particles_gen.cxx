//////////////////////////////////////////////////////////////
// EMCAL Barrel detector
// Generate particles with particle gun
// J.Kim 04/02/21
// M.Zurek 05/05/21
//////////////////////////////////////////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

#include "TMath.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include <random>
#include <fmt/core.h>

#include "emcal_barrel_common_functions.h"

using namespace HepMC3;

void emcal_barrel_particles_gen(std::string out_fname, int n_events = 1e6, double e_start = 0.0, double e_end = 20.0, std::string particle_name = "electron") {
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom* r1 = new TRandom();

  // Constraining the solid angle, but larger than that subtended by the
  // detector
  // https://indico.bnl.gov/event/7449/contributions/35966/attachments/27177/41430/EIC-DWG-Calo-03192020.pdf
  // See a figure on slide 26
  double cos_theta_min = std::cos(M_PI * (45.0 / 180.0));
  double cos_theta_max = std::cos(M_PI * (135.0 / 180.0));

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 2212 - proton
    GenParticlePtr p1 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum
    Double_t p        = r1->Uniform(e_start, e_end);
    Double_t phi      = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t costheta = r1->Uniform(cos_theta_min, cos_theta_max);
    Double_t theta    = std::acos(costheta);
    Double_t px       = p * std::cos(phi) * std::sin(theta);
    Double_t py       = p * std::sin(phi) * std::sin(theta);
    Double_t pz       = p * std::cos(theta);
    
    // type 1 is final state
    auto [id, mass] = extract_particle_parameters(particle_name);
    GenParticlePtr p3 = std::make_shared<GenParticle>(FourVector(px, py, pz, sqrt(p * p + (mass * mass))), id, 1);

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

