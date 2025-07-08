// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "SimHitAnalysis.h"

namespace benchmarks {

  // AlgorithmInit
  //---------------------------------------------------------------------------
  void SimHitAnalysis::AlgorithmInit(
      std::shared_ptr<spdlog::logger>& logger
      )
  {
    m_log = logger;

    // initialize histograms
    m_nphot_dist = new TH1D("nphot_dist", "Number of incident photons;N_{photons}",
        Tools::nphot_max, 0, Tools::nphot_max
        );
    m_nphot_vs_p = new TH2D("nphot_vs_p", "Number of incident photons vs. Thrown Momentum;p [GeV];N_{photons}",
        Tools::momentum_bins, 0, Tools::momentum_max,
        Tools::nphot_max,     0, Tools::nphot_max
        );
    m_nphot_vs_p__transient = new TH1D("nphot_vs_p__transient", "",
        m_nphot_vs_p->GetNbinsX(),
        m_nphot_vs_p->GetXaxis()->GetXmin(),
        m_nphot_vs_p->GetXaxis()->GetXmax()
        );
    m_phot_spectrum = new TH1D(
        "phot_spectrum_sim",
        "Incident photon wavelength;#lambda [nm], at emission point",
        Tools::wl_bins, Tools::wl_min, Tools::wl_max
        );
    m_phot_spectrum->SetLineColor(kViolet+2);
    m_phot_spectrum->SetFillColor(kViolet+2);
    m_nphot_dist->SetLineColor(kCyan+2);
    m_nphot_dist->SetFillColor(kCyan+2);
  }


  // AlgorithmProcess
  //---------------------------------------------------------------------------
  void SimHitAnalysis::AlgorithmProcess(
      const edm4hep::MCParticleCollection&    mc_parts,
      const edm4hep::SimTrackerHitCollection& sim_hits
      )
  {

    // get thrown momentum from a single MCParticle
    auto part = std::make_unique<ChargedParticle>(m_log);
    part->SetSingleParticle(mc_parts);
    auto thrown_momentum = part->GetMomentum();

    // fill nphot distribution
    m_nphot_dist->Fill(sim_hits.size());

    // loop over hits
    for(const auto& sim_hit : sim_hits) {

      // - fill 1D thrown momentum histogram `m_nphot_vs_p__transient` for each (pre-digitized) sensor hit
      // FIXME: get thrown momentum from multi-track events; but in this attempt
      // below the parent momentum is not consistent; instead for now we use
      // `thrown_momentum` from truth above
      /*
      auto photon = sim_hit.getMCParticle();
      if(photon.parents_size()>0) thrown_momentum = edm4hep::utils::p(photon.getParents(0));
      */
      m_nphot_vs_p__transient->Fill(thrown_momentum);

      // fill photon spectrum
      auto wavelength = Tools::GetPhotonWavelength(sim_hit, m_log);
      if(wavelength>=0)
        m_phot_spectrum->Fill(wavelength);
    }

    // use `m_nphot_vs_p__transient` results to fill 2D hist `m_nphot_vs_p` for this event
    for(int b=1; b<=m_nphot_vs_p__transient->GetNbinsX(); b++) {
      auto nphot = m_nphot_vs_p__transient->GetBinContent(b);
      if(nphot>0) {
        auto momentum = m_nphot_vs_p__transient->GetBinCenter(b);
        m_nphot_vs_p->Fill(momentum,nphot);
      }
    }

    // clear `m_nphot_vs_p__transient` to be ready for the next event
    m_nphot_vs_p__transient->Reset();

  }


  // AlgorithmFinish
  //---------------------------------------------------------------------------
  void SimHitAnalysis::AlgorithmFinish() {
    // delete transient histograms, so they don't get written
    delete m_nphot_vs_p__transient;
  }

}
