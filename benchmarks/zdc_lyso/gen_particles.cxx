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

// Generate single electron as input to the Insert simulation.
// --
// We generate events with a constant polar angle with respect to
// the proton direction and then rotate so that the events are given
// in normal lab coordinate system
// --
void gen_particles(
                    int n_events = 10, 
                    const char* out_fname = "gen_particles.hepmc", 
                    TString particle_name = "e-",
                    double th_min = 3., // Minimum polar angle, in degrees
		    double th_max = 3., // Maximum polar angle, in degrees
		    double phi_min = 0., // Minimum azimuthal angle, in degrees
                    double phi_max = 360., // Maximum azimuthal angle, in degrees
                    double p = 10.,  // Momentum in GeV/c
		    int dist = 0  //Momentum distribution: 0=fixed, 1=uniform, 2=Gaussian
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
  const double mass = particle->Mass();
  const int pdgID = particle->PdgCode();

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {

    //Set the event number
    evt.set_event_number(events_parsed);

    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 111 - pi0
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum with respect to proton direction
    double phi   = r1->Uniform(phi_min*TMath::DegToRad(),phi_max*TMath::DegToRad());
    double th    = r1->Uniform(th_min*TMath::DegToRad(),th_max*TMath::DegToRad());

    //Total momentum distribution
    double pevent = -1;
    if(dist==0){ //fixed
	pevent = p;
    }
    else if(dist==1){ //Uniform: +-50% variation
	pevent = p*(1. + r1->Uniform(-0.5,0.5) );
    }
    else if(dist==2){  //Gaussian: Sigma = 0.1*mean
	while(pevent<0) //Avoid negative values
		pevent = r1->Gaus(p,0.1*p);
    }

    double px    = pevent * std::cos(phi) * std::sin(th);
    double py    = pevent * std::sin(phi) * std::sin(th);
    double pz    = pevent * std::cos(th);
    TVector3 pvec(px,py,pz); 

    //Rotate to lab coordinate system
    double cross_angle = -25./1000.; //in Rad
    TVector3 pbeam_dir(sin(cross_angle),0,cos(cross_angle)); //proton beam direction
    pvec.RotateY(-pbeam_dir.Theta()); // Theta is returned positive, beam in negative X
    // type 1 is final state
    // pdgid 11 - electron 0.510 MeV/c^2
    GenParticlePtr p3 = std::make_shared<GenParticle>(
        FourVector(
            pvec.X(), pvec.Y(), pvec.Z(),
            sqrt(pevent*pevent + (mass * mass))),
        pdgID, 1);

    //If wanted, set non-zero vertex
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
    if (events_parsed % 100 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
