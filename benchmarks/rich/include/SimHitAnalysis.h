// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#pragma once

#include <spdlog/spdlog.h>

#include <TH1D.h>
#include <TH2D.h>

#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/MCParticleCollection.h>

#include "Tools.h"
#include "ChargedParticle.h"

namespace benchmarks {

  class SimHitAnalysis {

    public:
      SimHitAnalysis() = default;
      ~SimHitAnalysis() {}

      // algorithm methods
      void AlgorithmInit(std::shared_ptr<spdlog::logger>& logger);
      void AlgorithmProcess(
          const edm4hep::MCParticleCollection&    mc_parts,
          const edm4hep::SimTrackerHitCollection& sim_hits
          );
      void AlgorithmFinish();

    private:

      // histograms
      TH1D *m_nphot_dist;
      TH2D *m_nphot_vs_p;
      TH1D *m_nphot_vs_p__transient; // transient (not written)
      TH1D *m_phot_spectrum;

      // logger
      std::shared_ptr<spdlog::logger> m_log;

  };

}
