// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#pragma once

#include <spdlog/spdlog.h>

#include <TH2D.h>
#include <TGraphErrors.h>

#include <edm4eic/MCRecoParticleAssociationCollection.h>
#include <edm4hep/utils/kinematics.h>

#include "Tools.h"

namespace benchmarks {

  class ReconstructedParticleAnalysis {

    public:
      ReconstructedParticleAnalysis() = default;
      ~ReconstructedParticleAnalysis() {}

      // algorithm methods
      void AlgorithmInit(std::shared_ptr<spdlog::logger>& logger);
      void AlgorithmProcess(
          const edm4eic::MCRecoParticleAssociationCollection& in_assocs
          );
      void AlgorithmFinish();

    private:

      // plots
      TGraphErrors *m_correct_vs_p;
      TH2D         *m_correct_vs_p_transient; // (not written)

      // additional objects
      std::shared_ptr<spdlog::logger> m_log;

  };

}
