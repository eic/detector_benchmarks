/*
 * HistogramsSim.h
 *
 *  Created on: 21 apr 2025
 *      Author: Sam Corey
 */

#ifndef HISTOGRAMSNEUTRONTHRESHOLDS_H_
#define HISTOGRAMSNEUTRONTHRESHOLDS_H_

#include <TH1.h>
#include <TH1D.h>

#include <TH2.h>
#include <TH2D.h>

#include <TMath.h>

using namespace TMath;

void CreateHistogamsNeutronThresholds();

// HcalEndcapNHit

TH1D *h_nHCal_hit_contrib_energy;
TH1D *h_nHCal_hit_contrib_time;

TH2D *h_nHCal_hit_contrib_energy_vs_time;
TH2D *h_nHCal_hit_contrib_energy_vs_telap;
TH2D *h_nHCal_hit_contrib_energy_vs_time_total;


void CreateHistogamsNeutronThresholds()
{

	// HcalEndcapNHit
	h_nHCal_hit_contrib_energy = new TH1D("h_nHCal_hit_contrib_time", "Backwards HCal Hit Contribution Time; t; counts", 50, 0, 1000);
	h_nHCal_hit_contrib_time = new TH1D("h_nHCal_hit_contrib_energy", "Backwards HCal Hit Contribution Energy; Energy [GeV]; counts", 50, 0, 0.00000002);

	h_nHCal_hit_contrib_energy_vs_time = new TH2D("h_nHCal_hit_contrib_2D_E_vs_t", "Backwards HCal Hit Contribution Energy vs. time; Energy [GeV]; time [ns]; counts", 50, 0, 0.00000002, 50, 0, 1000);
	h_nHCal_hit_contrib_energy_vs_telap = new TH2D("h_nHCal_hit_contrib_energy_vs_telap", "Backwards HCal Hit Contribution Energy vs. Time; Energy threshold [GeV]; t max; counts", 50, 0, 0.0005, 50, 0, 1100);
	h_nHCal_hit_contrib_energy_vs_time_total = new TH2D("h_nHCal_hit_contrib_energy_vs_time_total", "Backwards HCal Hit Contribution Energy vs. Time; Energy threshold [GeV]; t max; counts", 50, 0, 0.0005, 50, 0, 1100);

}

void DeleteHistogamsNeutronThresholds()
{

	delete h_nHCal_hit_contrib_energy;
	delete h_nHCal_hit_contrib_time;

	delete h_nHCal_hit_contrib_energy_vs_time;
	delete h_nHCal_hit_contrib_energy_vs_telap;
	delete h_nHCal_hit_contrib_energy_vs_time_total;

}


#endif /* HISTOGRAMSDEPTHCHECK_H_ */
