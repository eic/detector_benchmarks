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
#include <TMath.h>

using namespace HepMC3;

std::tuple<double, int, double> GetParticleInfo(TDatabasePDG* pdg, TString particle_name)
{
  TParticlePDG *particle = pdg->GetParticle(particle_name);
  const double mass = particle->Mass();
  const int pdgID = particle->PdgCode();
  const double lifetime = particle->Lifetime();
  return std::make_tuple(mass, pdgID, lifetime);
}
// Calculates the decay length of a particle. Samples from an exponential decay.
double GetDecayLength(TRandom3* r1, double lifetime, double mass, double momentum_magnitude)
{ 
  double c_speed = TMath::C() * 1000.; // speed of light im mm/sec
  double average_decay_length = (momentum_magnitude/mass) * lifetime * c_speed;
  return r1->Exp(average_decay_length);
}

// Generate single pi0 mesons and decay them to 2 photons
void gen_pi0_decay(int n_events = 100000, UInt_t seed = 0, char* out_fname = "pi0_decay.hepmc",
		      double p_min = 60., // in GeV/c
		      double p_max = 275.) // in GeV/c
{

  const double theta_min = 0.0; // in mRad
  const double theta_max = 4; // in mRad
  //const double p_min = 100.; // in GeV/c
  //const double p_max = 275.; // in GeV/c

  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom3 *r1 = new TRandom3(seed); //Default = 0, which uses clock to set seed
  cout<<"Random number seed is "<<r1->GetSeed()<<"!"<<endl;

  // Getting generated particle information
  TDatabasePDG *pdg = new TDatabasePDG();
    
  auto pi0_info = GetParticleInfo(pdg, "pi0");
  double pi0_mass = std::get<0>(pi0_info);
  int pi0_pdgID = std::get<1>(pi0_info);
  double pi0_lifetime = std::get<2>(pi0_info);

  auto photon_info = GetParticleInfo(pdg, "gamma");
  double photon_mass = std::get<0>(photon_info);
  int photon_pdgID = std::get<1>(photon_info);

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
    Double_t pi0_p     = r1->Uniform(p_min, p_max);
    Double_t pi0_phi   = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t pi0_th    = r1->Uniform(theta_min/1000., theta_max/1000.); // Divide by 1000 for radians
    Double_t pi0_px    = pi0_p * TMath::Cos(pi0_phi) * TMath::Sin(pi0_th);
    Double_t pi0_py    = pi0_p * TMath::Sin(pi0_phi) * TMath::Sin(pi0_th);
    Double_t pi0_pz    = pi0_p * TMath::Cos(pi0_th);
    Double_t pi0_E     = TMath::Sqrt(pi0_p*pi0_p + pi0_mass*pi0_mass);

    // Rotate to lab coordinate system
    TVector3 pi0_pvec(pi0_px, pi0_py, pi0_pz); 
    double cross_angle = -25./1000.; // in Rad
    TVector3 pbeam_dir(TMath::Sin(cross_angle), 0, TMath::Cos(cross_angle)); //proton beam direction
    pi0_pvec.RotateY(cross_angle); // Theta is returned positive, beam in negative X

    // type 2 is state that will decay
    GenParticlePtr p_pi0 = std::make_shared<GenParticle>(
        FourVector(pi0_pvec.X(), pi0_pvec.Y(), pi0_pvec.Z(), pi0_E), pi0_pdgID, 2 );
    
    // Generating pi0 particle, will be generated at origin
    // Must have input electron + proton for vertex
    GenVertexPtr pi0_initial_vertex = std::make_shared<GenVertex>();
    pi0_initial_vertex->add_particle_in(p1);
    pi0_initial_vertex->add_particle_in(p2);
    pi0_initial_vertex->add_particle_out(p_pi0);
    evt.add_vertex(pi0_initial_vertex);

    // Generate gamma1 + gamma1 in pi0 rest frame
    TLorentzVector gamma1_rest, gamma2_rest;

    // Generating uniformly along a sphere
    double cost_gamma1_rest = r1->Uniform(-1,1);
    double th_gamma1_rest = TMath::ACos(cost_gamma1_rest);
    double sint_gamma1_rest = TMath::Sin(th_gamma1_rest);

    double phi_gamma1_rest = r1->Uniform(-1.*TMath::Pi(),1.*TMath::Pi());
    double cosp_gamma1_rest = TMath::Cos(phi_gamma1_rest);
    double sinp_gamma1_rest = TMath::Sin(phi_gamma1_rest);

    // Calculate energy of each particle in the pi0 rest frame
    // This is half the pi0 mass
    double E_gamma1_rest = pi0_mass/2;
    double E_gamma2_rest = pi0_mass/2;

    // Both particles will have the same momentum, and are massless
    double momentum_rest = pi0_mass/2;

    gamma1_rest.SetE(E_gamma1_rest);
    gamma1_rest.SetPx( momentum_rest * sint_gamma1_rest * cosp_gamma1_rest );
    gamma1_rest.SetPy( momentum_rest * sint_gamma1_rest * sinp_gamma1_rest );
    gamma1_rest.SetPz( momentum_rest * cost_gamma1_rest );

    gamma2_rest.SetE(E_gamma2_rest);
    gamma2_rest.SetPx( -gamma1_rest.Px() );
    gamma2_rest.SetPy( -gamma1_rest.Py() );
    gamma2_rest.SetPz( -gamma1_rest.Pz() );

    // Boost both gammas to lab frame
    TLorentzVector pi0_lab(pi0_pvec.X(), pi0_pvec.Y(), pi0_pvec.Z(), pi0_E);
    TVector3 pi0_boost = pi0_lab.BoostVector();
    TLorentzVector gamma1_lab, gamma2_lab;  
    gamma1_lab = gamma1_rest; 
    gamma1_lab.Boost(pi0_boost);
    gamma2_lab = gamma2_rest;
    gamma2_lab.Boost(pi0_boost);

    // Calculating position for pi0 decay
    TVector3 pi0_unit = pi0_lab.Vect().Unit();
    double pi0_decay_length = GetDecayLength(r1, pi0_lifetime, pi0_mass, pi0_lab.P());
    TVector3 pi0_decay_position = pi0_unit * pi0_decay_length;
    double pi0_decay_time = pi0_decay_length / pi0_lab.Beta() ; // Decay time in lab frame in length units (mm)
  
    // Generating vertex for pi0 decay
    GenParticlePtr p_gamma1 = std::make_shared<GenParticle>(
      FourVector(gamma1_lab.Px(), gamma1_lab.Py(), gamma1_lab.Pz(), gamma1_lab.E()), photon_pdgID, 1 );

    GenParticlePtr p_gamma2 = std::make_shared<GenParticle>(
      FourVector(gamma2_lab.Px(), gamma2_lab.Py(), gamma2_lab.Pz(), gamma2_lab.E()), photon_pdgID, 1 );

    GenVertexPtr v_pi0_decay = std::make_shared<GenVertex>(FourVector(pi0_decay_position.X(), pi0_decay_position.Y(), pi0_decay_position.Z(), pi0_decay_time));
    v_pi0_decay->add_particle_in(p_pi0);
    v_pi0_decay->add_particle_out(p_gamma1);
    v_pi0_decay->add_particle_out(p_gamma2);

    evt.add_vertex(v_pi0_decay);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }
    double zdc_z=35800;
    TVector3 extrap_gamma1=pi0_decay_position+gamma1_lab.Vect()*((zdc_z-pbeam_dir.Dot(pi0_decay_position))/(pbeam_dir.Dot(gamma1_lab.Vect())));
    TVector3 extrap_gamma2=pi0_decay_position+gamma2_lab.Vect()*((zdc_z-pbeam_dir.Dot(pi0_decay_position))/(pbeam_dir.Dot(gamma2_lab.Vect())));
    if (extrap_gamma1.Angle(pbeam_dir)<0.004 && extrap_gamma2.Angle(pbeam_dir)<0.004 && pi0_decay_position.Dot(pbeam_dir)<zdc_z)
      hepmc_output.write_event(evt);
    if (events_parsed % 1000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();

  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
