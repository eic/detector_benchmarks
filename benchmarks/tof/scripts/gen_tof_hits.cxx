#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include <iostream>
#include<random>
#include<cmath>
#include <math.h>
#include <TMath.h>

using namespace HepMC3;

/** Generate multiple electrons/positron tracks in the central region.
 *  This is for testing detectors in the "barrel" region.
 */
void gen_tof_hits(int n_events = 100,
                         const char* out_fname = "tof_hits.hepmc",
                         int n_parts = 2)
{
  double cos_theta_min = std::cos( 1.0*(M_PI/180.0));
  double cos_theta_max = std::cos(189.0*(M_PI/180.0));

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
    Double_t p     = r1->Uniform(0.2, 4);
    Double_t costh = cos(9.385/180*3.1415926); //r1->Uniform(cos_theta_min, cos_theta_max);
    Double_t th    = std::acos(costh);
    Double_t phi   = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t px    = p * std::cos(phi) * std::sin(th);
    Double_t py    = p * std::sin(phi) * std::sin(th);
    Double_t pz    = (events_parsed % 2 == 0) ? p * std::cos(th) : -1 * p * std::cos(th);

    for (int ip = 0; ip < n_parts; ip++) {
      GenParticlePtr p1 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
      GenParticlePtr p2 = std::make_shared<GenParticle>(FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

      // Define momentum
        phi   = r1->Uniform(0.0, 2.0 * M_PI);
        px    = p * std::cos(phi) * std::sin(th);
        py    = p * std::sin(phi) * std::sin(th);
        pz    = p * std::cos(th);
        GenParticlePtr p3 = std::make_shared<GenParticle>(FourVector(px, py, pz, sqrt(p * p + (0.000511 * 0.000511))),
                                                            ((ip % 2 == 0) ? 11 : -11), 1);

        phi   = r1->Uniform(0.0, 2.0 * M_PI);
        px    = p * std::cos(phi) * std::sin(th);
        py    = p * std::sin(phi) * std::sin(th);
        pz    = p * std::cos(th);
        GenParticlePtr p4 = std::make_shared<GenParticle>(FourVector(px, py, pz, sqrt(p * p + (0.938272 * 0.938272))),
                                                          ((ip % 2 == 0) ? 2212 : -2212), 1);

        phi   = r1->Uniform(0.0, 2.0 * M_PI);
        px    = p * std::cos(phi) * std::sin(th);
        py    = p * std::sin(phi) * std::sin(th);
        pz    = p * std::cos(th);
        GenParticlePtr p5 = std::make_shared<GenParticle>(FourVector(px, py, pz, sqrt(p * p + (0.493677 * 0.493677))),
                                                          ((ip % 2 == 0) ? 321 : -321), 1);
        phi   = r1->Uniform(0.0, 2.0 * M_PI);
        px    = p * std::cos(phi) * std::sin(th);
        py    = p * std::sin(phi) * std::sin(th);
        pz    = p * std::cos(th);
        GenParticlePtr p6 = std::make_shared<GenParticle>(FourVector(px, py, pz, sqrt(p * p + (0.139570 * 0.139570))),
                                                          ((ip % 2 == 0) ? 211 : -211), 1);

      GenVertexPtr v1 = std::make_shared<GenVertex>();
      v1->add_particle_in(p1);
      v1->add_particle_in(p2);
      v1->add_particle_out(p3);
      v1->add_particle_out(p4);
      v1->add_particle_out(p5);
      v1->add_particle_out(p6);
      evt.add_vertex(v1);
    }

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
