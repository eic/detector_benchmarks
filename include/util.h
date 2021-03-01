#ifndef UTIL_H
#define UTIL_H

// TODO: should probably be moved to a global benchmark utility library

#include <algorithm>
#include <cmath>
#include <exception>
#include <fmt/core.h>
#include <limits>
#include <string>
#include <vector>

#include <Math/Vector4D.h>

#include "dd4pod/Geant4ParticleCollection.h"
#include "eicd/TrackParametersCollection.h"

namespace util {

  // Exception definition for unknown particle errors
  // FIXME: A utility exception base class should be included in the analysis
  //        utility library, so we can skip most of this boilerplate
  class unknown_particle_error : public std::exception {
  public:
    unknown_particle_error(std::string_view particle) : m_particle{particle} {}
    virtual const char* what() const throw()
    {
      return fmt::format("Unknown particle type: {}", m_particle).c_str();
    }
    virtual const char* type() const throw() { return "unknown_particle_error"; }

  private:
    const std::string m_particle;
  };

  // Simple function to return the appropriate PDG mass for the particles
  // we care about for this process.
  // FIXME: consider something more robust (maybe based on hepPDT) to the
  //        analysis utility library
  inline double get_pdg_mass(std::string_view part)
  {
    if (part == "electron") {
      return 0.0005109989461;
    } else if (part == "muon") {
      return .1056583745;
    } else if (part == "jpsi") {
      return 3.0969;
    } else if (part == "upsilon") {
      return 9.49630;
    } else if (part == "proton"){
      return 0.938272;
    } else {
      throw unknown_particle_error{part};
    }
  }

  // Get a vector of 4-momenta from raw tracking info, using an externally
  // provided particle mass assumption.
  inline auto momenta_from_tracking(const std::vector<eic::TrackParametersData>& tracks,
                                    const double                                 mass)
  {
    std::vector<ROOT::Math::PxPyPzMVector> momenta{tracks.size()};
    // transform our raw tracker info into proper 4-momenta
    std::transform(tracks.begin(), tracks.end(), momenta.begin(), [mass](const auto& track) {
      // make sure we don't divide by zero
      if (fabs(track.qOverP) < 1e-9) {
        return ROOT::Math::PxPyPzMVector{};
      }
      const double p  = fabs(1. / track.qOverP);
      const double px = p * cos(track.phi) * sin(track.theta);
      const double py = p * sin(track.phi) * sin(track.theta);
      const double pz = p * cos(track.theta);
      return ROOT::Math::PxPyPzMVector{px, py, pz, mass};
    });
    return momenta;
  }

  // Get a vector of 4-momenta from the simulation data.
  // TODO: Add PID selector (maybe using ranges?)
  inline auto momenta_from_simulation(const std::vector<dd4pod::Geant4ParticleData>& parts)
  {
    std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
    // transform our simulation particle data into 4-momenta
    std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
      return ROOT::Math::PxPyPzMVector{part.psx, part.psy, part.psz, part.mass};
    });
    return momenta;
  }

  // Find the decay pair candidates from a vector of particles (parts),
  // with invariant mass closest to a desired value (pdg_mass)
  inline std::pair<ROOT::Math::PxPyPzMVector, ROOT::Math::PxPyPzMVector>
  find_decay_pair(const std::vector<ROOT::Math::PxPyPzMVector>& parts, const double pdg_mass)
  {
    int    first     = -1;
    int    second    = -1;
    double best_mass = -1;

    // go through all particle combinatorics, calculate the invariant mass
    // for each combination, and remember which combination is the closest
    // to the desired pdg_mass
    for (size_t i = 0; i < parts.size(); ++i) {
      for (size_t j = i + 1; j < parts.size(); ++j) {
        const double new_mass{(parts[i] + parts[j]).mass()};
        if (fabs(new_mass - pdg_mass) < fabs(best_mass - pdg_mass)) {
          first     = i;
          second    = j;
          best_mass = new_mass;
        }
      }
    }
    if (first < 0) {
      return {{}, {}};
    }
    return {parts[first], parts[second]};
  }

  // Calculate the magnitude of the momentum of a vector of 4-vectors
  inline auto mom(const std::vector<ROOT::Math::PxPyPzMVector>& momenta)
  {
    std::vector<double> P(momenta.size());
    // transform our raw tracker info into proper 4-momenta
    std::transform(momenta.begin(), momenta.end(), P.begin(),
                   [](const auto& mom) { return mom.P(); });
    return P;
  }
  // Calculate the transverse momentum of a vector of 4-vectors
  inline auto pt(const std::vector<ROOT::Math::PxPyPzMVector>& momenta)
  {
    std::vector<double> pt(momenta.size());
    // transform our raw tracker info into proper 4-momenta
    std::transform(momenta.begin(), momenta.end(), pt.begin(),
                   [](const auto& mom) { return mom.pt(); });
    return pt;
  }

  // Calculate the azimuthal angle phi of a vector of 4-vectors
  inline auto phi(const std::vector<ROOT::Math::PxPyPzMVector>& momenta)
  {
    std::vector<double> phi(momenta.size());
    // transform our raw tracker info into proper 4-momenta
    std::transform(momenta.begin(), momenta.end(), phi.begin(),
                   [](const auto& mom) { return mom.phi(); });
    return phi;
  }
  // Calculate the pseudo-rapidity of a vector of particles
  inline auto eta(const std::vector<ROOT::Math::PxPyPzMVector>& momenta)
  {
    std::vector<double> eta(momenta.size());
    // transform our raw tracker info into proper 4-momenta
    std::transform(momenta.begin(), momenta.end(), eta.begin(),
                   [](const auto& mom) { return mom.eta(); });
    return eta;
  }

  //=========================================================================================================

} // namespace util

#endif
