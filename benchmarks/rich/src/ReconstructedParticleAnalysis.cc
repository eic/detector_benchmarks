// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#include "ReconstructedParticleAnalysis.h"

namespace benchmarks {

  // AlgorithmInit
  //---------------------------------------------------------------------------
  void ReconstructedParticleAnalysis::AlgorithmInit(std::shared_ptr<spdlog::logger>& logger) {
    m_log = logger;

    // initialize plots
    m_correct_vs_p_transient = new TH2D("correct_vs_p_transient", "Fraction with Correct PDG vs. Thrown Momentum;p [GeV];fraction",
        Tools::momentum_bins, 0, Tools::momentum_max,
        2,                    0, 2
        );
    m_correct_vs_p = new TGraphErrors();
    m_correct_vs_p->SetName("correct_vs_p");
    m_correct_vs_p->SetTitle(m_correct_vs_p_transient->GetTitle());
    m_correct_vs_p->GetXaxis()->SetTitle(m_correct_vs_p_transient->GetXaxis()->GetTitle());
    m_correct_vs_p->GetYaxis()->SetTitle(m_correct_vs_p_transient->GetYaxis()->GetTitle());
    m_correct_vs_p->SetMarkerStyle(kFullCircle);
  }


  // AlgorithmProcess
  //---------------------------------------------------------------------------
  void ReconstructedParticleAnalysis::AlgorithmProcess(
      const edm4eic::MCRecoParticleAssociationCollection& in_assocs
      )
  {
    if(!in_assocs.isValid()) { m_log->error("invalid input collection 'in_assocs'"); return; }

    // loop over input associations
    for(const auto& assoc : in_assocs) {

      // get particles
      m_log->trace("{:-^50}"," Particle ");
      auto& simpart = assoc.getSim();
      auto& recpart = assoc.getRec();
      if(!recpart.isAvailable()) { m_log->warn("reconstructed particle not available for this association"); continue; }
      if(!simpart.isAvailable()) { m_log->warn("simulated particle not available for this association");     continue; }

      // get PDG values
      auto simpart_pdg      = simpart.getPDG(); // from MC truth
      auto recpart_pdg      = recpart.getPDG(); // PDG member of ReconstructedParticle
      auto recpart_pid_used = recpart.getParticleIDUsed();
      auto pid_pdg          = recpart_pid_used.isAvailable() ? recpart_pid_used.getPDG() : 0;
      if(!recpart_pid_used.isAvailable())
        m_log->error("reconstructed particle does not have particleIDUsed relation");
      m_log->trace("PDGs: sim.PDG | rec.PDG | PDG from PID  =  {:^6} | {:^6} | {:^6}", simpart_pdg, recpart_pdg, pid_pdg);

      // check goodness of PID: if near zero, assume PID was not used for this particle and skip
      if(recpart.getGoodnessOfPID()<0.01) {
        m_log->warn("skip particle with low goodnessOfPID (likely has no PID)");
        continue;
      }

      // get momenta
      auto simpart_p = edm4hep::utils::p(simpart);
      // auto recpart_p = edm4hep::utils::p(recpart); // unused

      // fill plots
      m_correct_vs_p_transient->Fill(
          simpart_p,
          simpart_pdg == pid_pdg ? 1 : 0
          );

    } // end loop over input associations

  }


  // AlgorithmFinish
  //---------------------------------------------------------------------------
  void ReconstructedParticleAnalysis::AlgorithmFinish() {

    // compute fraction of correct PDG for each x-axis bin
    for(int bx=1; bx<=m_correct_vs_p_transient->GetNbinsX(); bx++) {
      auto n_wrong  = m_correct_vs_p_transient->GetBinContent(bx,1);
      auto n_right  = m_correct_vs_p_transient->GetBinContent(bx,2);
      auto n_total  = n_wrong + n_right;
      auto momentum = m_correct_vs_p_transient->GetXaxis()->GetBinCenter(bx);
      if(n_total>0) {
        auto frac = n_right / n_total;
        m_correct_vs_p->AddPoint(momentum, frac);
      }
    }

    // delete transient histograms, so they don't get written
    delete m_correct_vs_p_transient;

    // write objects which are not automatically saved
    m_correct_vs_p->Write();
  }

}
