// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#pragma once

#include <spdlog/spdlog.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <edm4eic/RawTrackerHitCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>

#include "Tools.h"

namespace benchmarks {

  class RawHitAnalysis {

    public:

      RawHitAnalysis() = default;
      ~RawHitAnalysis() {}

      // algorithm methods
      void AlgorithmInit(std::shared_ptr<spdlog::logger>& logger);
      void AlgorithmProcess(
          const edm4eic::RawTrackerHitCollection& raw_hits,
          const edm4eic::MCRecoTrackerHitAssociationCollection& assocs
          );
      void AlgorithmFinish();

    private:

      // binning
      const int adc_max = std::pow(2,10);
      const int tdc_max = std::pow(2,10);

      // histograms
      TH1D *m_adc_dist;
      TH1D *m_tdc_dist;
      TH2D *m_tdc_vs_adc;
      TH1D *m_phot_spectrum;
      TH1D *m_nhits_dist;

      // logging
      std::shared_ptr<spdlog::logger> m_log;

  };

}
