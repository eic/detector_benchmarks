// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#pragma once

#include <spdlog/spdlog.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/utils/vector_utils.h>

namespace benchmarks {

  class ChargedParticle {

    public:

      ChargedParticle(std::shared_ptr<spdlog::logger> log) : m_log(log) {};
      ~ChargedParticle() {}

      // accessors
      double GetMomentum() { return m_momentum; };
      double GetEta()      { return m_eta;      };
      int    GetPDG()      { return m_pdg;      };

      // set info from a single MCParticle
      void SetSingleParticle(const edm4hep::MCParticle& mc_part);
      // set info from a (thrown) MCParticle of a collection
      void SetSingleParticle(const edm4hep::MCParticleCollection& mc_parts);

    private:

      // members
      double m_momentum;
      double m_eta;
      int    m_pdg = 0;

      // logger
      std::shared_ptr<spdlog::logger> m_log;

  };

}
