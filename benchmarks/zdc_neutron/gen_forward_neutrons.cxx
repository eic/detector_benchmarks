#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TRandom3.h"
#include "TVector3.h"

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <iostream>
#include <random>
#include <cmath>
#include <math.h>
#include <TMath.h>

using namespace HepMC3;

//Generate single neutrons in the forward detector region.
void gen_forward_neutrons(int n_events = 10000, UInt_t seed = 0, const char* out_fname = "fwd_neutrons.hepmc")
{

  double theta_min = 0.0; //in mRad
  double theta_max = 6.5; //in mRad
  double cost_min = std::cos(theta_max/1000.) ; //Minimum value of cos(theta)
  double cost_max = std::cos(theta_min/1000.) ; //Maximum value of cos(theta)

  double p_min = 100.; //in GeV/c
  double p_max = 100.; //in GeV/c

  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom3 *r1 = new TRandom3(seed); //Default = 0, which uses clock to set seed
  cout<<"Random number seed is "<<r1->GetSeed()<<"!"<<endl;

  // Getting generated particle information
  TDatabasePDG *pdg = new TDatabasePDG();
  TParticlePDG *particle = pdg->GetParticle("neutron");
  const double mass = particle->Mass();
  const int pdgID = particle->PdgCode();

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {

    //Set the event number
    evt.set_event_number(events_parsed);

    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum with respect to EIC proton beam direction
    Double_t p     = r1->Uniform(p_min,p_max); //Uniform momentum generation
    Double_t phi   = r1->Uniform(0.0, 2.0 * M_PI); //Uniform phi generation

    Double_t cost = r1->Uniform(cost_min,cost_max); //Uniform cos(theta) generation
    Double_t th    = std::acos(cost);
    
    Double_t px    = p * std::cos(phi) * std::sin(th);
    Double_t py    = p * std::sin(phi) * std::sin(th);
    Double_t pz    = p * std::cos(th);
    Double_t E     = sqrt(p*p + mass*mass);

    //Rotate to lab coordinate system
    TVector3 pvec(px,py,pz); 
    double cross_angle = -25./1000.; //in Rad
    TVector3 pbeam_dir(sin(cross_angle),0,cos(cross_angle)); //proton beam direction
    pvec.RotateY(-pbeam_dir.Theta()); // Theta is returned positive, beam in negative X

    // type 1 is final state
    GenParticlePtr p3 = std::make_shared<GenParticle>(
        FourVector(pvec.X(), pvec.Y(), pvec.Z(), E), pdgID, 1 );

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
