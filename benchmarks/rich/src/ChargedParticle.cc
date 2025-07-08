// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "ChargedParticle.h"

namespace benchmarks {

  // set info from a single MCParticle
  void ChargedParticle::SetSingleParticle(const edm4hep::MCParticle& mc_part) {
    m_momentum = edm4hep::utils::magnitude(mc_part.getMomentum());
    m_eta      = edm4hep::utils::eta(mc_part.getMomentum());
    m_pdg      = mc_part.getPDG();
  }

  // set info from a (thrown) MCParticle of a collection
  void ChargedParticle::SetSingleParticle(const edm4hep::MCParticleCollection& mc_parts) {
    int n_thrown = 0;
    for(const auto& mc_part : mc_parts) {
      if(mc_part.getGeneratorStatus() == 1) {
        this->SetSingleParticle(mc_part);
        n_thrown++;
      }
    }
    if(n_thrown == 0) m_log->warn("ChargedParticle::SetSingleParticle: no particle found");
    if(n_thrown >  1) m_log->warn("ChargedParticle::SetSingleParticle: found {} particles (more than one)", n_thrown);
  }

}
