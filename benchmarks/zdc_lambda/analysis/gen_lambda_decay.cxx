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

// Generate single lambda mesons and decay them to a neutron + 2 photons
void gen_lambda_decay(int n_events = 100000, UInt_t seed = 0, char* out_fname = "lambda_decay.hepmc",
		      double p_min = 100., // in GeV/c
		      double p_max = 275.) // in GeV/c
{

  const double theta_min = 0.0; // in mRad
  const double theta_max = 3.0; // in mRad
  //const double p_min = 100.; // in GeV/c
  //const double p_max = 275.; // in GeV/c

  WriterAscii hepmc_output(out_fname);
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom3 *r1 = new TRandom3(seed); //Default = 0, which uses clock to set seed
  cout<<"Random number seed is "<<r1->GetSeed()<<"!"<<endl;

  // Getting generated particle information
  TDatabasePDG *pdg = new TDatabasePDG();
  
  auto lambda_info = GetParticleInfo(pdg, "Lambda0");
  double lambda_mass = std::get<0>(lambda_info);
  int lambda_pdgID = std::get<1>(lambda_info);
  double lambda_lifetime = std::get<2>(lambda_info);

  auto neutron_info = GetParticleInfo(pdg, "neutron");
  double neutron_mass = std::get<0>(neutron_info);
  int neutron_pdgID = std::get<1>(neutron_info);

  auto pi0_info = GetParticleInfo(pdg, "pi0");
  double pi0_mass = std::get<0>(pi0_info);
  int pi0_pdgID = std::get<1>(pi0_info);
  double pi0_lifetime = std::get<2>(pi0_info);

  auto photon_info = GetParticleInfo(pdg, "gamma");
  double photon_mass = std::get<0>(photon_info);
  int photon_pdgID = std::get<1>(photon_info);

  int accepted_events = 0;
  while (accepted_events < n_events) {

    //Set the event number
    evt.set_event_number(accepted_events);

    // FourVector(px,py,pz,e,pdgid,status)
    // type 4 is beam
    // pdgid 11 - electron
    // pdgid 2212 - proton
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    // Define momentum with respect to EIC proton beam direction
    Double_t lambda_p     = r1->Uniform(p_min, p_max);
    Double_t lambda_phi   = r1->Uniform(0.0, 2.0 * M_PI);
    Double_t lambda_th    = r1->Uniform(theta_min/1000., theta_max/1000.); // Divide by 1000 for radians
    Double_t lambda_px    = lambda_p * TMath::Cos(lambda_phi) * TMath::Sin(lambda_th);
    Double_t lambda_py    = lambda_p * TMath::Sin(lambda_phi) * TMath::Sin(lambda_th);
    Double_t lambda_pz    = lambda_p * TMath::Cos(lambda_th);
    Double_t lambda_E     = TMath::Sqrt(lambda_p*lambda_p + lambda_mass*lambda_mass);

    // Rotate to lab coordinate system
    TVector3 lambda_pvec(lambda_px, lambda_py, lambda_pz); 
    double cross_angle = -25./1000.; // in Rad
    TVector3 pbeam_dir(TMath::Sin(cross_angle), 0, TMath::Cos(cross_angle)); //proton beam direction
    lambda_pvec.RotateY(cross_angle); // Theta is returned positive, beam in negative X

    // type 2 is state that will decay
    GenParticlePtr p_lambda = std::make_shared<GenParticle>(
        FourVector(lambda_pvec.X(), lambda_pvec.Y(), lambda_pvec.Z(), lambda_E), lambda_pdgID, 2 );
    
    // Generating lambda particle, will be generated at origin
    // Must have input electron + proton for vertex
    GenVertexPtr lambda_initial_vertex = std::make_shared<GenVertex>();
    lambda_initial_vertex->add_particle_in(p1);
    lambda_initial_vertex->add_particle_in(p2);
    lambda_initial_vertex->add_particle_out(p_lambda);
    evt.add_vertex(lambda_initial_vertex);

    // Generate neutron + pi0 in lambda rest frame
    TLorentzVector neutron_rest, pi0_rest;

    // Generating uniformly along a sphere
    double cost_neutron_rest = r1->Uniform(-1,1);
    double th_neutron_rest = TMath::ACos(cost_neutron_rest);
    double sint_neutron_rest = TMath::Sin(th_neutron_rest);

    double phi_neutron_rest = r1->Uniform(-1.*TMath::Pi(),1.*TMath::Pi());
    double cosp_neutron_rest = TMath::Cos(phi_neutron_rest);
    double sinp_neutron_rest = TMath::Sin(phi_neutron_rest);

    // Calculate energy of each particle in the lambda rest frame
    // See problem 3.19 in Introduction to Elementary Particles, 2nd edition by D. Griffiths
    double E_neutron_rest = (-TMath::Power(pi0_mass, 2.) + TMath::Power(lambda_mass, 2.) + TMath::Power(neutron_mass, 2.) ) / (2. * lambda_mass) ;
    double E_pi0_rest = (-TMath::Power(neutron_mass, 2.) + TMath::Power(lambda_mass, 2.) + TMath::Power(pi0_mass, 2.) ) / (2. * lambda_mass) ;

    // Both particles will have the same momentum, so just use neutron variables
    double momentum_rest = TMath::Sqrt( E_neutron_rest*E_neutron_rest - neutron_mass*neutron_mass );

    neutron_rest.SetE(E_neutron_rest);
    neutron_rest.SetPx( momentum_rest * sint_neutron_rest * cosp_neutron_rest );
    neutron_rest.SetPy( momentum_rest * sint_neutron_rest * sinp_neutron_rest );
    neutron_rest.SetPz( momentum_rest * cost_neutron_rest );

    pi0_rest.SetE(E_pi0_rest);
    pi0_rest.SetPx( -neutron_rest.Px() );
    pi0_rest.SetPy( -neutron_rest.Py() );
    pi0_rest.SetPz( -neutron_rest.Pz() );

    // Boost neutron & pion to lab frame
    TLorentzVector lambda_lab(lambda_pvec.X(), lambda_pvec.Y(), lambda_pvec.Z(), lambda_E);
    TVector3 lambda_boost = lambda_lab.BoostVector();
    TLorentzVector neutron_lab, pi0_lab;  
    neutron_lab = neutron_rest; 
    neutron_lab.Boost(lambda_boost);
    pi0_lab = pi0_rest;
    pi0_lab.Boost(lambda_boost);

    // Calculating position for lambda decay
    TVector3 lambda_unit = lambda_lab.Vect().Unit();
    double lambda_decay_length = GetDecayLength(r1, lambda_lifetime, lambda_mass, lambda_lab.P());
    TVector3 lambda_decay_position = lambda_unit * lambda_decay_length;
    double lambda_decay_time = lambda_decay_length / lambda_lab.Beta() ; // Decay time in lab frame in length units (mm)
  
    // Generating vertex for lambda decay
    GenParticlePtr p_neutron = std::make_shared<GenParticle>(
      FourVector(neutron_lab.Px(), neutron_lab.Py(), neutron_lab.Pz(), neutron_lab.E()), neutron_pdgID, 1 );

    GenParticlePtr p_pi0 = std::make_shared<GenParticle>(
      FourVector(pi0_lab.Px(), pi0_lab.Py(), pi0_lab.Pz(), pi0_lab.E()), pi0_pdgID, 2 );

    GenVertexPtr v_lambda_decay = std::make_shared<GenVertex>(FourVector(lambda_decay_position.X(), lambda_decay_position.Y(), lambda_decay_position.Z(), lambda_decay_time));
    v_lambda_decay->add_particle_in(p_lambda);
    v_lambda_decay->add_particle_out(p_neutron);
    v_lambda_decay->add_particle_out(p_pi0);

    evt.add_vertex(v_lambda_decay);

    // Generate two photons from pi0 decay
    TLorentzVector gamma1_rest, gamma2_rest;

    // Generating uniformly along a sphere
    double cost_gamma1_rest = r1->Uniform(-1,1);
    double th_gamma1_rest = TMath::ACos(cost_gamma1_rest);
    double sint_gamma1_rest = TMath::Sin(th_gamma1_rest);

    double phi_gamma1_rest = r1->Uniform(-1.*TMath::Pi(),1.*TMath::Pi());
    double cosp_gamma1_rest = TMath::Cos(phi_gamma1_rest);
    double sinp_gamma1_rest = TMath::Sin(phi_gamma1_rest);

    // Photons are massless so they each get equal energies
    gamma1_rest.SetE(pi0_mass/2.);
    gamma1_rest.SetPx( (pi0_mass/2.)*sint_gamma1_rest*cosp_gamma1_rest );
    gamma1_rest.SetPy( (pi0_mass/2.)*sint_gamma1_rest*sinp_gamma1_rest );
    gamma1_rest.SetPz( (pi0_mass/2.)*cost_gamma1_rest );

    gamma2_rest.SetE(pi0_mass/2.);
    gamma2_rest.SetPx( -gamma1_rest.Px() );
    gamma2_rest.SetPy( -gamma1_rest.Py() );
    gamma2_rest.SetPz( -gamma1_rest.Pz() );

    // Boost neutron & pion to lab frame
    TVector3 pi0_boost = pi0_lab.BoostVector();
    TLorentzVector gamma1_lab, gamma2_lab;
    gamma1_lab = gamma1_rest; 
    gamma1_lab.Boost(pi0_boost);
    gamma2_lab = gamma2_rest; 
    gamma2_lab.Boost(pi0_boost);
  
    GenParticlePtr p_gamma1 = std::make_shared<GenParticle>(
      FourVector(gamma1_lab.Px(), gamma1_lab.Py(), gamma1_lab.Pz(), gamma1_lab.E()), photon_pdgID, 1 );

    GenParticlePtr p_gamma2 = std::make_shared<GenParticle>(
      FourVector(gamma2_lab.Px(), gamma2_lab.Py(), gamma2_lab.Pz(), gamma2_lab.E()), photon_pdgID, 1 );

    // Generate pi0 at same position as the lambda. Approximating pi0 decay as instantaneous
    GenVertexPtr v_pi0_decay = std::make_shared<GenVertex>(FourVector(lambda_decay_position.X(), lambda_decay_position.Y(), lambda_decay_position.Z(), lambda_decay_time));
    v_pi0_decay->add_particle_in(p_pi0);
    v_pi0_decay->add_particle_out(p_gamma1);
    v_pi0_decay->add_particle_out(p_gamma2);

    //std::cout<<  lambda_pvec.Angle(pbeam_dir) << " " << neutron_lab.Angle(pbeam_dir) << " " << gamma1_lab.Angle(pbeam_dir) << " " << gamma2_lab.Angle(pbeam_dir) << std::endl;
    
    evt.add_vertex(v_pi0_decay);

    if (accepted_events == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }
    double zdc_z=35800;
    TVector3 extrap_neutron=lambda_decay_position+neutron_lab.Vect()*((zdc_z-pbeam_dir.Dot(lambda_decay_position))/(pbeam_dir.Dot(neutron_lab.Vect())));
    TVector3 extrap_gamma1=lambda_decay_position+gamma1_lab.Vect()*((zdc_z-pbeam_dir.Dot(lambda_decay_position))/(pbeam_dir.Dot(gamma1_lab.Vect())));
    TVector3 extrap_gamma2=lambda_decay_position+gamma2_lab.Vect()*((zdc_z-pbeam_dir.Dot(lambda_decay_position))/(pbeam_dir.Dot(gamma2_lab.Vect())));
    if (extrap_neutron.Angle(pbeam_dir)<0.004 && extrap_gamma1.Angle(pbeam_dir)<0.004 && extrap_gamma2.Angle(pbeam_dir)<0.004 && lambda_decay_position.Dot(pbeam_dir)<zdc_z)
      accepted_events++;
      hepmc_output.write_event(evt);
    if (accepted_events % 1000 == 0) {
      std::cout << "Event: " << accepted_events << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();

  std::cout << "Events parsed and written: " << accepted_events << std::endl;
}
