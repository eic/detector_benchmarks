#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>
#include <random>
#include <cmath>
#include <math.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <vector>

using namespace HepMC3;

// Generate multi-particle events for PFA validation
// Event types:
//   0 = di-jet (back-to-back particles)
//   1 = multi-hadron (multiple pions/kaons/protons)
//   2 = EM shower (e+e-gamma combinations)
void gen_multiparticle(
                    int n_events = 1000, 
                    const char* out_fname = "gen_multiparticle.hepmc", 
                    int event_type = 0,  // 0=di-jet, 1=multi-hadron, 2=EM shower
                    double total_energy = 20.,  // Total energy in GeV
                    int useCrossingAngle = 1  // 0 = no rotation, 1 = -25 mrad crossing angle
                  )
{ 
  WriterAscii hepmc_output(out_fname);
  int events_parsed = 0;
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom3 *r1 = new TRandom3(0); //Use time as random seed
  
  // Getting particle database
  TDatabasePDG *pdg = new TDatabasePDG();

  std::cout << "Generating " << n_events << " multi-particle events" << std::endl;
  std::cout << "Event type: " << event_type << std::endl;
  std::cout << "Total energy: " << total_energy << " GeV" << std::endl;
  std::cout << "Crossing angle: " << (useCrossingAngle ? "-25 mrad" : "0") << std::endl;

  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {

    //Set the event number
    evt.set_event_number(events_parsed);

    // Beam particles for event record
    GenParticlePtr p1 =
        std::make_shared<GenParticle>(FourVector(0.0, 0.0, 10.0, 10.0), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, 0.0, 0.938), 2212, 4);

    GenVertexPtr v1 = std::make_shared<GenVertex>();
    evt.shift_position_by(FourVector(0., 0., 0., 0.));
    v1->add_particle_in(p1);
    v1->add_particle_in(p2);

    // Crossing angle rotation
    double cross_angle = -25. / 1000. * useCrossingAngle; //in Rad
    TVector3 pbeam_dir(sin(cross_angle), 0, cos(cross_angle));

    if (event_type == 0) {
      // Di-jet: back-to-back particles
      // Split energy between two jets
      double e1 = total_energy * r1->Uniform(0.3, 0.7);
      double e2 = total_energy - e1;
      
      // Random direction for first jet
      double phi1 = r1->Uniform(0, 2 * TMath::Pi());
      double th1 = r1->Uniform(10 * TMath::DegToRad(), 170 * TMath::DegToRad());
      
      // Second jet back-to-back
      double phi2 = phi1 + TMath::Pi();
      double th2 = TMath::Pi() - th1;
      
      // Particle types: mostly pions with some kaons
      std::vector<TString> hadrons = {"pi+", "pi-", "pi0", "K+", "K-"};
      
      // First jet (2-4 particles)
      int n1 = r1->Integer(2, 5);
      for (int i = 0; i < n1; i++) {
        TString pname = hadrons[r1->Integer(hadrons.size())];
        TParticlePDG *part = pdg->GetParticle(pname);
        double mass = part->Mass();
        int pdgID = part->PdgCode();
        
        // Energy split with some spread
        double e_part = (e1 / n1) * r1->Gaus(1.0, 0.2);
        e_part = std::max(e_part, mass + 0.1); // Ensure kinematically allowed
        double p_part = sqrt(e_part * e_part - mass * mass);
        
        // Angular spread around jet axis
        double phi = phi1 + r1->Gaus(0, 0.2);
        double th = th1 + r1->Gaus(0, 0.1);
        
        double px = p_part * std::cos(phi) * std::sin(th);
        double py = p_part * std::sin(phi) * std::sin(th);
        double pz = p_part * std::cos(th);
        TVector3 pvec(px, py, pz);
        pvec.RotateY(-pbeam_dir.Theta());
        
        GenParticlePtr p_out = std::make_shared<GenParticle>(
            FourVector(pvec.X(), pvec.Y(), pvec.Z(), e_part), pdgID, 1);
        v1->add_particle_out(p_out);
      }
      
      // Second jet (2-4 particles)
      int n2 = r1->Integer(2, 5);
      for (int i = 0; i < n2; i++) {
        TString pname = hadrons[r1->Integer(hadrons.size())];
        TParticlePDG *part = pdg->GetParticle(pname);
        double mass = part->Mass();
        int pdgID = part->PdgCode();
        
        double e_part = (e2 / n2) * r1->Gaus(1.0, 0.2);
        e_part = std::max(e_part, mass + 0.1);
        double p_part = sqrt(e_part * e_part - mass * mass);
        
        double phi = phi2 + r1->Gaus(0, 0.2);
        double th = th2 + r1->Gaus(0, 0.1);
        
        double px = p_part * std::cos(phi) * std::sin(th);
        double py = p_part * std::sin(phi) * std::sin(th);
        double pz = p_part * std::cos(th);
        TVector3 pvec(px, py, pz);
        pvec.RotateY(-pbeam_dir.Theta());
        
        GenParticlePtr p_out = std::make_shared<GenParticle>(
            FourVector(pvec.X(), pvec.Y(), pvec.Z(), e_part), pdgID, 1);
        v1->add_particle_out(p_out);
      }
      
    } else if (event_type == 1) {
      // Multi-hadron: isotropic distribution of hadrons
      std::vector<TString> hadrons = {"pi+", "pi-", "pi0", "K+", "K-", "proton", "neutron"};
      
      int n_particles = r1->Integer(3, 8);
      for (int i = 0; i < n_particles; i++) {
        TString pname = hadrons[r1->Integer(hadrons.size())];
        TParticlePDG *part = pdg->GetParticle(pname);
        double mass = part->Mass();
        int pdgID = part->PdgCode();
        
        double e_part = (total_energy / n_particles) * r1->Gaus(1.0, 0.3);
        e_part = std::max(e_part, mass + 0.1);
        double p_part = sqrt(e_part * e_part - mass * mass);
        
        double phi = r1->Uniform(0, 2 * TMath::Pi());
        double costh = r1->Uniform(-0.9, 0.9); // Avoid very forward/backward
        double th = acos(costh);
        
        double px = p_part * std::cos(phi) * std::sin(th);
        double py = p_part * std::sin(phi) * std::sin(th);
        double pz = p_part * std::cos(th);
        TVector3 pvec(px, py, pz);
        pvec.RotateY(-pbeam_dir.Theta());
        
        GenParticlePtr p_out = std::make_shared<GenParticle>(
            FourVector(pvec.X(), pvec.Y(), pvec.Z(), e_part), pdgID, 1);
        v1->add_particle_out(p_out);
      }
      
    } else if (event_type == 2) {
      // EM shower: e+, e-, gamma combinations
      int n_particles = r1->Integer(2, 6);
      
      for (int i = 0; i < n_particles; i++) {
        TString pname;
        double rand = r1->Uniform(0, 1);
        if (rand < 0.33) pname = "e-";
        else if (rand < 0.66) pname = "e+";
        else pname = "gamma";
        
        TParticlePDG *part = pdg->GetParticle(pname);
        double mass = part->Mass();
        int pdgID = part->PdgCode();
        
        double e_part = (total_energy / n_particles) * r1->Gaus(1.0, 0.3);
        e_part = std::max(e_part, mass + 0.1);
        double p_part = sqrt(e_part * e_part - mass * mass);
        
        // EM showers tend to be more collimated
        double phi = r1->Uniform(0, 2 * TMath::Pi());
        double th = r1->Gaus(90 * TMath::DegToRad(), 20 * TMath::DegToRad());
        th = std::max(10 * TMath::DegToRad(), std::min(170 * TMath::DegToRad(), th));
        
        double px = p_part * std::cos(phi) * std::sin(th);
        double py = p_part * std::sin(phi) * std::sin(th);
        double pz = p_part * std::cos(th);
        TVector3 pvec(px, py, pz);
        pvec.RotateY(-pbeam_dir.Theta());
        
        GenParticlePtr p_out = std::make_shared<GenParticle>(
            FourVector(pvec.X(), pvec.Y(), pvec.Z(), e_part), pdgID, 1);
        v1->add_particle_out(p_out);
      }
    }

    evt.add_vertex(v1);

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    if (events_parsed % 100 == 0 && events_parsed > 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;
}
