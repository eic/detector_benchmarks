// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "RawHitAnalysis.h"
#include <edm4eic/EDM4eicVersion.h>

namespace benchmarks {

  // AlgorithmInit
  //---------------------------------------------------------------------------
  void RawHitAnalysis::AlgorithmInit(
      std::shared_ptr<spdlog::logger>& logger
      )
  {
    m_log = logger;

    // define histograms
    m_adc_dist = new TH1D("adc_dist", "ADC Distribution;counts", adc_max, 0, adc_max);
    m_tdc_dist = new TH1D("tdc_dist", "TDC Distribution;counts", tdc_max, 0, tdc_max);
    m_tdc_vs_adc = new TH2D("tdc_vs_adc", "TDC vs. ADC;ADC counts;TDC counts",
        adc_max/2, 0, adc_max/2,
        tdc_max/2, 0, tdc_max/2
        );
    m_phot_spectrum = new TH1D(
        "phot_spectrum_rec",
        "Photon wavelength for digitized hits;#lambda [nm], at emission point",
        Tools::wl_bins, Tools::wl_min, Tools::wl_max
        );
    m_nhits_dist = new TH1D("nhits_dist", "Number of digitized hits;N_{hits}",
        Tools::npe_max, 0, Tools::npe_max
        );

    // format histograms
    auto format1D = [] (auto h) {
      h->SetLineColor(kBlack);
      h->SetFillColor(kBlack);
    };
    format1D(m_adc_dist);
    format1D(m_tdc_dist);
    m_phot_spectrum->SetLineColor(kViolet+2);
    m_phot_spectrum->SetFillColor(kViolet+2);
    m_nhits_dist->SetLineColor(kCyan+2);
    m_nhits_dist->SetFillColor(kCyan+2);
  }


  // AlgorithmProcess
  //---------------------------------------------------------------------------
  void RawHitAnalysis::AlgorithmProcess(
      const edm4eic::RawTrackerHitCollection& raw_hits,
      const edm4eic::MCRecoTrackerHitAssociationCollection& assocs
      )
  {

    // fill nhits
    m_nhits_dist->Fill(raw_hits.size());

    // loop over all raw hits (including noise)
    for(const auto& raw_hit : raw_hits) {
      auto adc = raw_hit.getCharge();
      auto tdc = raw_hit.getTimeStamp();
      m_adc_dist->Fill(adc);
      m_tdc_dist->Fill(tdc);
      m_tdc_vs_adc->Fill(adc,tdc);
    }

    // loop over hits with associations (no noise)
    for(const auto& assoc : assocs) {
#if EDM4EIC_VERSION_MAJOR >= 6
      auto hit = assoc.getSimHit();
#else
    for(const auto& hit : assoc.getSimHits()) {
#endif
      auto wavelength = Tools::GetPhotonWavelength(hit, m_log);
      if(wavelength>=0)
        m_phot_spectrum->Fill(wavelength);
#if EDM4EIC_VERSION_MAJOR >= 6
#else
    }
#endif
    }
  }


  // AlgorithmFinish
  //---------------------------------------------------------------------------
  void RawHitAnalysis::AlgorithmFinish() {
  }

}
