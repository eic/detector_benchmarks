// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#pragma once

// general
#include <map>
#include <spdlog/spdlog.h>

// ROOT
#include <TVector3.h>

// IRT
#include <IRT/ParametricSurface.h>

// DD4hep
#include <Evaluator/DD4hepUnits.h>

// data model
#include <edm4hep/EDM4hepVersion.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/utils/kinematics.h>

namespace benchmarks {

  // Tools class, filled with miscellanous helper functions
  class Tools {
    public:

      // -------------------------------------------------------------------------------------
      // physics constants

      // h*c constant, for wavelength <=> energy conversion [GeV*nm]
      // obtained from `dd4hep::h_Planck * dd4hep::c_light / (dd4hep::GeV * dd4hep::nm)`
      static constexpr double HC = 1.239841875e-06;

      // opticalphoton PDG
      static constexpr int PHOTON_PDG = -22;

      // -------------------------------------------------------------------------------------
      // PDG mass lookup

      // local PDG mass database
      // FIXME: cannot use `TDatabasePDG` since it is not thread safe; until we
      // have a proper PDG database service, we hard-code the masses we need;
      // use Tools::GetPDGMass for access
      static std::unordered_map<int,double> GetPDGMasses() {
        return std::unordered_map<int,double>{
          { 11,   0.000510999 },
          { 211,  0.13957     },
          { 321,  0.493677    },
          { 2212, 0.938272    }
        };
      }

      static double GetPDGMass(int pdg) {
        double mass;
        try { mass = GetPDGMasses().at(std::abs(pdg)); }
        catch(const std::out_of_range& e) {
          throw std::runtime_error(fmt::format("RUNTIME ERROR: unknown PDG={} in algorithms/pid/Tools::GetPDGMass",pdg));
        }
        return mass;
      }

      static int GetNumPDGs() { return GetPDGMasses().size(); };


      // -------------------------------------------------------------------------------------
      // get photon wavelength
      static double GetPhotonWavelength(const edm4hep::SimTrackerHit& hit, std::shared_ptr<spdlog::logger> log) {
#if EDM4HEP_BUILD_VERSION >= EDM4HEP_VERSION(0, 99, 0)
        auto phot = hit.getParticle();
#else
        auto phot = hit.getMCParticle();
#endif
        if(!phot.isAvailable()) {
          log->error("no MCParticle in hit");
          return -1.0;
        }
        if(phot.getPDG() != PHOTON_PDG) {
          log->warn("hit MCParticle is not an opticalphoton; PDG is {}; ignoring", phot.getPDG());
          return -1.0;
        }
        auto momentum   = edm4hep::utils::p(phot);
        auto wavelength = momentum > 0 ? HC / momentum : 0.0;
        return wavelength;
      }

      // -------------------------------------------------------------------------------------
      // binning
      static constexpr int    n_bins         = 100;
      static constexpr int    momentum_bins  = 500;
      static constexpr double momentum_max   = 70;
      static constexpr int    eta_bins       = 50;
      static constexpr double eta_min        = 0;
      static constexpr double eta_max        = 5;
      static constexpr int    npe_max        = 50;
      static constexpr int    nphot_max      = 3000;
      static constexpr int    theta_bins     = 1500;
      static constexpr double theta_max      = 300;
      static constexpr double thetaResid_max = 100;
      static constexpr int    phi_bins       = 100;
      static constexpr int    wl_bins        = 120;
      static constexpr double wl_min         = 0;
      static constexpr double wl_max         = 1200;

  }; // class Tools
} // namespace benchmarks
