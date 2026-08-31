// NeutronThresholdsUtil.h
#ifndef NEUTRON_THRESHOLDS_UTIL_H
#define NEUTRON_THRESHOLDS_UTIL_H

#include <vector>
#include <iostream>
#include <TH2D.h>


// Constants for the neutron threshold analysis
namespace NeutronThresholds {
    const double E_MIP = 0.00075;

    // Time thresholds in ns
    const std::vector<double> T_MAX = {25, 50, 100, 200, 500};

    // Energy thresholds based on MIP
    const std::vector<double> E_TH = {0.5 * E_MIP, 0.25 * E_MIP, 0.1 * E_MIP, 0.05 * E_MIP, 0};

    // Initialize the hits passed counters
    inline std::vector<std::vector<double>> createHitsPassedMatrix() {
        return std::vector<std::vector<double>>(5, std::vector<double>(5, 0));
    }

    // Process the contributions for a hit
    template<typename ContribCollection>
    void processContributions(const ContribCollection& contrib, std::vector<std::vector<double>>& hits_passed, std::vector<std::vector<double>>& hits_passed_telap, 
        TH1* h_contrib_time = nullptr, TH1* h_contrib_energy = nullptr, double E_TH_correction = 1.0, bool debug = false
    ) {
        std::vector<double> E_sum(5, 0);
        std::vector<double> E_sum_telap(5, 0);
        double t_min = 1100;

        for (unsigned c = 0; c < contrib.size(); ++c) {
            if (contrib.at(c).getTime() < t_min) { t_min = contrib.at(c).getTime(); }
        }

        if(debug) std::cout << contrib.at(0).getTime() << std::endl;

        for (unsigned c = 0; c < contrib.size(); ++c) {
            if(debug) std::cout << "hit time = " << contrib.at(c).getTime() << std::endl;

            // double t_min_per_contrib = std::min(1100, contrib.at(c).getTime());

            if (h_contrib_time) h_contrib_time->Fill(contrib.at(c).getTime());
            if (h_contrib_energy) h_contrib_energy->Fill(contrib.at(c).getEnergy());

            for (int itm = 0; itm < T_MAX.size(); itm++) {
                if (contrib.at(c).getTime() < T_MAX[itm]) {
                    E_sum[itm] += contrib.at(c).getEnergy();
                }
                if ((contrib.at(c).getTime() - t_min) < T_MAX[itm]) {
                    E_sum_telap[itm] += contrib.at(c).getEnergy();
                }
            }
        }

        for (int iet = 0; iet < E_TH.size(); iet++) {
            for (int ies = 0; ies < E_sum.size(); ies++) {
                if (E_sum[ies] > E_TH[iet]*E_TH_correction) {
                    hits_passed[iet][ies] += 1;
                }
                if (E_sum_telap[ies] > E_TH[iet]*E_TH_correction) {
                    hits_passed_telap[iet][ies] += 1;
                }
            }
        }
    }

    // Fill histograms with the hits_passed data
    void fillThresholdHistograms(
        const std::vector<std::vector<double>>& hits_passed,
        const std::vector<std::vector<double>>& hits_passed_telap,
        TH2* h_energy_vs_time,
        TH2* h_energy_vs_telap,
        TH2* h_energy_vs_time_total,
        double E_TH_correction = 1.0
    ) {
        for (int iet = 0; iet < E_TH.size(); iet++) {
            for (int itm = 0; itm < T_MAX.size(); itm++) {
                if (hits_passed[iet][itm] > 0) {
                    h_energy_vs_time->Fill(E_TH[iet]*E_TH_correction, T_MAX[itm]);
                }
                if (hits_passed_telap[iet][itm] > 0) {
                    h_energy_vs_telap->Fill(E_TH[iet]*E_TH_correction, T_MAX[itm]);
                }
                h_energy_vs_time_total->Fill(E_TH[iet]*E_TH_correction, T_MAX[itm]);
            }
        }
    }
}

#endif // NEUTRON_THRESHOLDS_UTIL_H
