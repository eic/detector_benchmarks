#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <iostream>
#include <random>
#include <cmath>
#include <math.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

using namespace HepMC3;

// Generate single particles for PFA validation
// Particles: e-, e+, pi-, pi+, gamma, etc.
// Angular ranges: forward (2-28°), barrel (45-135°), backward (130-177°)
// Fixed momentum per event (for resolution studies)
void gen_particles(
                    int n_events = 10000, 
                    const char* out_fname = "gen_particles.hepmc", 
                    TString particle_name = "e-",
                    double th_min = 3., // Minimum polar angle, in degrees
		    double th_max = 28., // Maximum polar angle, in degrees
		    double phi_min = 0., // Minimum azimuthal angle, in degrees
                    double phi_max = 360., // Maximum azimuthal angle, in degrees
                    double p = 10.,  // Momentum in GeV/c
		    int useCrossingAngle = 1  // 0 = no rotation, 1 = -25 mrad crossing angle
                  )
{ 
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom3 *r1 = new TRandom3(0); //Use time as random seed
  
  // Getting generated particle information
  TDatabasePDG *pdg = new TDatabasePDG();
  TParticlePDG *particle = pdg->GetParticle(particle_name);
  
  if (!particle) {
    std::cerr << "ERROR: Unknown particle " << particle_name << std::endl;
    return;
  }
  
  const double mass = particle->Mass();
  const int pdgID = particle->PdgCode();

  std::cout << "Generating " << n_events << " events of " << particle_name 
            << " (PDG ID: " << pdgID << ", mass: " << mass << " GeV)" << std::endl;
  std::cout << "Momentum: " << p << " GeV/c" << std::endl;
  std::cout << "Polar angle range: " << th_min << " - " << th_max << " degrees" << std::endl;
  std::cout << "Azimuthal angle range: " << phi_min << " - " << phi_max << " degrees" << std::endl;
  std::cout << "Crossing angle: " << (useCrossingAngle ? "-25 mrad" : "0") << std::endl;

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {

    //Set the event number
    evt.set_event_number(events_parsed);

    // Beam particles for event record
    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum with respect to proton direction (before crossing angle rotation)
    double phi   = r1->Uniform(phi_min * TMath::DegToRad(), phi_max * TMath::DegToRad());
    double th    = r1->Uniform(th_min * TMath::DegToRad(), th_max * TMath::DegToRad());

    // Fixed momentum for resolution studies
    double pevent = p;

    double px    = pevent * std::cos(phi) * std::sin(th);
    double py    = pevent * std::sin(phi) * std::sin(th);
    double pz    = pevent * std::cos(th);
    TVector3 pvec(px, py, pz); 

    // Rotate to lab coordinate system with crossing angle
    double cross_angle = -25. / 1000. * useCrossingAngle; //in Rad
    TVector3 pbeam_dir(sin(cross_angle), 0, cos(cross_angle)); //proton beam direction
    pvec.RotateY(-pbeam_dir.Theta()); // Theta is returned positive, beam in negative X
    
    // type 1 is final state
    double energy = sqrt(pevent * pevent + (mass * mass));
    GenParticlePtr p3 = std::make_shared<GenParticle>(
        FourVector(pvec.X(), pvec.Y(), pvec.Z(), energy),
        pdgID, 1);

    // Vertex at origin
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    double vt = 0.;

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    evt.shift_position_by(FourVector(vx, vy, vz, vt));
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);
    v1->add_particle_out(p3);
    evt.add_vertex(v1);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 1000 == 0 && events_parsed > 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
