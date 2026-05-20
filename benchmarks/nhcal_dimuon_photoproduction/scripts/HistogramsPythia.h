/*
 * HistogramsPythia.h
 *
 *  Created on: 30 sie 2024
 *      Author: Khaless
 */

#ifndef HISTOGRAMSPYTHIA_H_
#define HISTOGRAMSPYTHIA_H_

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"

int CreateHistograms();


void MakeLogBins(double *array, int nbins, double binLo, double binHi)
{
	//double exp = (nbinEdg-1)/nBinPerDecimal;

    // Calculate the logarithmic bin edges
    double logMin = log10(binLo);
    double logMax = log10(binHi);
    double binWidth = (logMax - logMin) / nbins;

    for (int i = 0; i <= nbins; ++i) {
    	array[i] = pow(10, logMin + i * binWidth);
    }
}

// Event

TH1F *h_Events;
TH1F *h_Events_types;
TH1F *h_Events_cuts;

TH1F *h_Events_Diffractive;

TH1F *h_Events_nPartonsOut;

TH1F *h_XsecGen;
TH1F *h_XsecSel;

TH2F *h_XsecGen_err;
TH2F *h_XsecSel_err;

TH1F *h_Event_nPart_final;
TH1F *h_Event_nJets;
TH1F *h_Event_nJets_meas;
TH1F *h_Event_nJets_meas_no_nHCal;

TH2F *h_Event_xQ2;
TH2F *h_Event_yQ2;
TH2F *h_Event_xy;

TH2F *h_Event_nHCal_xQ2;
TH2F *h_Event_nHCal_yQ2;
TH2F *h_Event_nHCal_xy;

TH1F *h_Event_nPion_p;
TH1F *h_Event_nPion_n;
TH1F *h_Event_nKaon_p;
TH1F *h_Event_nKaon_n;
TH1F *h_Event_nProton_p;
TH1F *h_Event_nProton_n;
TH1F *h_Event_nElectron_p;
TH1F *h_Event_nElectron_n;

TH1F *h_Event_nNeutron;
TH1F *h_Event_nGamma;

// special


TH1F *h_Event_Q2;
TH1F *h_Event_x;
TH1F *h_Event_y;

TH1F *h_Event_jets_Q2;
TH1F *h_Event_jets_x;
TH1F *h_Event_jets_y;

TH1F *h_Event_nHCal_0_Q2;
TH1F *h_Event_nHCal_0_x;
TH1F *h_Event_nHCal_0_y;

TH1F *h_Event_nHCal_1_Q2;
TH1F *h_Event_nHCal_1_x;
TH1F *h_Event_nHCal_1_y;

TH1F *h_Event_nHCal_2_Q2;
TH1F *h_Event_nHCal_2_x;
TH1F *h_Event_nHCal_2_y;

TH1F *h_Event_AllHCal_Q2;
TH1F *h_Event_AllHCal_x;
TH1F *h_Event_AllHCal_y;


TH1F *h_Event_JetMeas_nHCal_0_Q2;
TH1F *h_Event_JetMeas_nHCal_0_x;
TH1F *h_Event_JetMeas_nHCal_0_y;

TH1F *h_Event_JetMeas_nHCal_1_Q2;
TH1F *h_Event_JetMeas_nHCal_1_x;
TH1F *h_Event_JetMeas_nHCal_1_y;

TH1F *h_Event_JetMeas_nHCal_2_Q2;
TH1F *h_Event_JetMeas_nHCal_2_x;
TH1F *h_Event_JetMeas_nHCal_2_y;

TH1F *h_Event_JetMeas_AllHCal_Q2;
TH1F *h_Event_JetMeas_AllHCal_x;
TH1F *h_Event_JetMeas_AllHCal_y;

TH1F *h_Event_JetMeas_no_nHCal_Q2;
TH1F *h_Event_JetMeas_no_nHCal_x;
TH1F *h_Event_JetMeas_no_nHCal_y;

TH2F *h_Event_HCal_jets;
TH2F *h_Event_HCal_jets_meas;
TH2F *h_Event_HCal_jets_meas_no_nHCal;

//full
TH2F *h_Event_HCal_jets_meas_full;


// Outgoing partons

TH1F *h_Partons_status;
TH1F *h_Partons_types;
TH1F *h_Partons_types_anti;

TH2F *h_Partons_eta;
TH2F *h_Partons_phi;
TH2F *h_Partons_p;
TH2F *h_Partons_pT;


TH2F *h_Parton_eta_pT;
TH2F *h_Parton_eta_p;
TH2F *h_Parton_eta_E;

TH2F *h_Parton_x_eta;
TH2F *h_Parton_y_eta;

TH2F *h_Parton_x_eta1;
TH2F *h_Parton_x_eta2;

TH2F *h_Parton_y_eta1;
TH2F *h_Parton_y_eta2;

// Particles

TH1F *h_Particle_eta;

TH2F *h_Particle_eta_p;
TH2F *h_Particle_eta_pT;
TH2F *h_Particle_eta_E;


// weighted
TH1F *h_Particle_eta_wE;

// eta, momentum
TH2F *h_Particle_pion_p_eta_p;
TH2F *h_Particle_pion_n_eta_p;
TH2F *h_Particle_Kaon_p_eta_p;
TH2F *h_Particle_Kaon_n_eta_p;
TH2F *h_Particle_proton_p_eta_p;
TH2F *h_Particle_proton_n_eta_p;
TH2F *h_Particle_Electron_p_eta_p;
TH2F *h_Particle_Electron_n_eta_p;

TH2F *h_Particle_Neutron_eta_p;
TH2F *h_Particle_Gamma_eta_p;

// eta, transverse momentum pT
TH2F *h_Particle_pion_p_eta_pT;
TH2F *h_Particle_pion_n_eta_pT;
TH2F *h_Particle_Kaon_p_eta_pT;
TH2F *h_Particle_Kaon_n_eta_pT;
TH2F *h_Particle_proton_p_eta_pT;
TH2F *h_Particle_proton_n_eta_pT;
TH2F *h_Particle_Electron_p_eta_pT;
TH2F *h_Particle_Electron_n_eta_pT;

TH2F *h_Particle_Neutron_eta_pT;
TH2F *h_Particle_Gamma_eta_pT;

// eta, energy
TH2F *h_Particle_Pion_p_eta_E;
TH2F *h_Particle_Pion_n_eta_E;
TH2F *h_Particle_Kaon_p_eta_E;
TH2F *h_Particle_Kaon_n_eta_E;
TH2F *h_Particle_Proton_p_eta_E;
TH2F *h_Particle_Proton_n_eta_E;
TH2F *h_Particle_Electron_p_eta_E;
TH2F *h_Particle_Electron_n_eta_E;

TH2F *h_Particle_Neutron_eta_E;
TH2F *h_Particle_Gamma_eta_E;


// Jets
TH1D *h_Jet_nPart;
TH1D *h_Jet_mass;
TH1D *h_Jet_charge;
TH1D *h_Jet_E;
TH1D *h_Jet_p;
TH1D *h_Jet_pT;
TH1D *h_Jet_eta;
TH1D *h_Jet_deta;

TH2F *h_Jets_eta;
TH2F *h_Jets_phi;
TH2F *h_Jets_p;
TH2F *h_Jets_pT;
TH2F *h_Jets_E;

// Jets measured
TH1D *h_Jet_meas_nPart;
TH1D *h_Jet_meas_mass;
TH1D *h_Jet_meas_charge;
TH1D *h_Jet_meas_E;
TH1D *h_Jet_meas_p;
TH1D *h_Jet_meas_pT;
TH1D *h_Jet_meas_eta;
TH1D *h_Jet_meas_deta;

TH2F *h_Jets_meas_eta;
TH2F *h_Jets_meas_phi;
TH2F *h_Jets_meas_p;
TH2F *h_Jets_meas_pT;
TH2F *h_Jets_meas_E;

TH1D *h_Jet_bHCal_part_eta;
TH1D *h_Jet_meas_bHCal_part_eta;

TH2D *h_Jet_HCal_part_eta;
TH2D *h_Jet_meas_HCal_part_eta;


// Jets measured without nHCal
TH1D *h_Jet_meas_no_nHCal_nPart;
TH1D *h_Jet_meas_no_nHCal_mass;
TH1D *h_Jet_meas_no_nHCal_charge;
TH1D *h_Jet_meas_no_nHCal_E;
TH1D *h_Jet_meas_no_nHCal_p;
TH1D *h_Jet_meas_no_nHCal_pT;
TH1D *h_Jet_meas_no_nHCal_eta;
TH1D *h_Jet_meas_no_nHCal_deta;

TH2F *h_Jets_meas_no_nHCal_eta;
TH2F *h_Jets_meas_no_nHCal_phi;
TH2F *h_Jets_meas_no_nHCal_p;
TH2F *h_Jets_meas_no_nHCal_pT;
TH2F *h_Jets_meas_no_nHCal_E;

// measured jets vs. partons
TH2F *h_Jets_meas_Partons_eta;
TH2F *h_Jets_meas_Partons_phi;
TH2F *h_Jets_meas_Partons_E;

TH2F *h_Jet_meas_Parton_eta1;
TH2F *h_Jet_meas_Parton_phi1;
TH2F *h_Jet_meas_Parton_E1;
TH2F *h_Jet_meas_Parton_eta2;
TH2F *h_Jet_meas_Parton_phi2;
TH2F *h_Jet_meas_Parton_E2;


// temp

TH1F *hist_eta_energy_tmp;
TH1F *hist_eta_energy_denom_tmp;




int CreateHistograms()
{


	const int nbins_x = 200;

	const int nbinEdg_x = nbins_x+1;

	double *logBinsArray_x = new double[nbins_x];

	MakeLogBins(logBinsArray_x, nbins_x, 10e-7, 1.0);


	// Event

	h_Events = new TH1F("h_Events", "Number of events; events; counts", 10, 0.0, 10.0);
	h_Events_types = new TH1F("h_Events_types", "Number of events; type; counts", 10, 0.0, 10.0);

	h_Events_types->GetXaxis()->SetBinLabel(1, "All events");
	h_Events_types->GetXaxis()->SetBinLabel(2, "1 parton in ePIC acc.");
	h_Events_types->GetXaxis()->SetBinLabel(3, "2 partons in ePIC acc.");
	h_Events_types->GetXaxis()->SetBinLabel(4, "1 jet in ePIC acc.");
	h_Events_types->GetXaxis()->SetBinLabel(5, "2 jets in ePIC acc.");
	h_Events_types->GetXaxis()->SetBinLabel(6, "2 or more jets in ePIC acc.");

	h_Events_cuts = new TH1F("h_Events_cuts", "Number of events; selection; counts", 10, 0.0, 10.0);

	h_Events_cuts->GetXaxis()->SetBinLabel(1, "Generated");
	h_Events_cuts->GetXaxis()->SetBinLabel(2, "Diffractive A or B");
	h_Events_cuts->GetXaxis()->SetBinLabel(3, "Hard diffractive A or B");
	h_Events_cuts->GetXaxis()->SetBinLabel(4, "2 partons in ePIC acc.");

	h_Events_Diffractive = new TH1F("h_Events_Diffractive", "Number of diffractive events; events; counts", 10, 0.0, 10.0);

	h_Events_Diffractive->GetXaxis()->SetBinLabel(1, "Diffractive A");
	h_Events_Diffractive->GetXaxis()->SetBinLabel(2, "Diffractive B");
	h_Events_Diffractive->GetXaxis()->SetBinLabel(3, "Hard diffractive A");
	h_Events_Diffractive->GetXaxis()->SetBinLabel(4, "Hard diffractive B");

	h_XsecGen = new TH1F("h_XsecGen", "Generated event cross-section; #sigma [mb]; counts", 10000, 0.0, 0.001);
	h_XsecSel = new TH1F("h_XsecSel", "Selected event cross-section; #sigma [mb]; counts", 10000, 0.0, 0.001);

	h_XsecGen_err = new TH2F("h_XsecGen_err", "Generated event cross-section vs. uncertainty; #sigma [mb]; #sigma_{#sigma} [mb]; counts", 10000, 0.0, 0.001, 1000, 0.0, 0.00001);
	h_XsecSel_err = new TH2F("h_XsecSel_err", "Selected event cross-section vs. uncertainty; #sigma [mb]; #sigma_{#sigma} [mb]; counts", 10000, 0.0, 0.001, 1000, 0.0, 0.00001);


	h_Events_nPartonsOut = new TH1F("h_Events_nPartonsOut", "Number of outgoing partons; N_{out} [1]; counts", 21, -0.5, 20.5);

	h_Event_nPart_final = new TH1F("h_Event_nPart_final", "Number of final MC particles; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nJets = new TH1F("h_Event_nJets", "Number of jets; N_{jet} [1]; counts", 21, -0.5, 20.5);
	h_Event_nJets_meas = new TH1F("h_Event_nJets_meas", "Number of measured jets; N_{jet} [1]; counts", 21, -0.5, 20.5);
	h_Event_nJets_meas_no_nHCal = new TH1F("h_Event_nJets_meas_no_nHCal", "Number of measured jets wihout nHCal; N_{jet} [1]; counts", 21, -0.5, 20.5);


	h_Event_xQ2 = new TH2F("h_Event_xQ2", "Event Q^{2} vs. x; x [1]; Q^{2} [GeV^{2}/c^{2}]; counts", nbins_x, logBinsArray_x, 500, 0.0, 5.0);
	h_Event_yQ2 = new TH2F("h_Event_yQ2", "Event Q^{2} vs. inelasticity y; y [1]; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 1.0, 500, 0.0, 5.0);
	h_Event_xy = new TH2F("h_Event_xy", "Event inelasticity y vs. x; x [1]; y [1]; counts", nbins_x, logBinsArray_x, 1000, 0.0, 1.0);


	h_Event_nHCal_xQ2 = new TH2F("h_Event_nHCal_xQ2", "Event with nHCal activity Q^{2} vs. x; x [1]; Q^{2} [GeV^{2}/c^{2}]; counts", nbins_x, logBinsArray_x, 500, 0.0, 5.0);
	h_Event_nHCal_yQ2 = new TH2F("h_Event_nHCal_yQ2", "Event with nHCal activity Q^{2} vs. inelasticity y; y [1]; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 1.0, 500, 0.0, 5.0);
	h_Event_nHCal_xy = new TH2F("h_Event_nHCal_xy", "Event with nHCal activity inelasticity y vs. x; x [1]; y [1]; counts", nbins_x, logBinsArray_x, 1000, 0.0, 1.0);


	h_Event_nPion_p = new TH1F("h_Event_nPion_p", "Number of MC particles #pi^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nPion_n = new TH1F("h_Event_nPion_n", "Number of MC particles #pi^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nKaon_p = new TH1F("h_Event_nKaon_p", "Number of MC particles K^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nKaon_n = new TH1F("h_Event_nKaon_n", "Number of MC particles K^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nProton_p = new TH1F("h_Event_nProton_p", "Number of MC particles p^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nProton_n = new TH1F("h_Event_nProton_n", "Number of MC particles p^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nElectron_p = new TH1F("h_Event_nElectron_p", "Number of MC particles e^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nElectron_n = new TH1F("h_Event_nElectron_n", "Number of MC particles e^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);

	h_Event_nNeutron = new TH1F("h_Event_nNeutron", "Number of MC particles n; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_Event_nGamma = new TH1F("h_Event_nGamma", "Number of MC particles #gamma; N_{MC} [1]; counts", 2001, -0.5, 2000.5);


	// special

	h_Event_Q2 = new TH1F("h_Event_Q2", "Event Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_x = new TH1F("h_Event_x", "Event x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_y = new TH1F("h_Event_y", "Event inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_jets_Q2 = new TH1F("h_Event_jets_Q2", "Event with jets>=2 Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_jets_x = new TH1F("h_Event_jets_x", "Event with jets>=2 x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_jets_y = new TH1F("h_Event_jets_y", "Event with jets>=2 inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_nHCal_0_Q2 = new TH1F("h_Event_nHCal_0_Q2", "Event with 0 jets in nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_nHCal_0_x = new TH1F("h_Event_nHCal_0_x", "Event with 0 jets in nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_nHCal_0_y = new TH1F("h_Event_nHCal_0_y", "Event with 0 jets in nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_nHCal_1_Q2 = new TH1F("h_Event_nHCal_1_Q2", "Event with 1 jet in nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_nHCal_1_x = new TH1F("h_Event_nHCal_1_x", "Event with 1 jet in nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_nHCal_1_y = new TH1F("h_Event_nHCal_1_y", "Event with 1 jet in nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_nHCal_2_Q2 = new TH1F("h_Event_nHCal_2_Q2", "Event with 2 jets in nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_nHCal_2_x = new TH1F("h_Event_nHCal_2_x", "Event with 2 jets in nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_nHCal_2_y = new TH1F("h_Event_nHCal_2_y", "Event with 2 jets in nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_AllHCal_Q2 = new TH1F("h_Event_AllHCal_Q2", "Event with jets in any HCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_AllHCal_x = new TH1F("h_Event_AllHCal_x", "Event with jets in any HCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_AllHCal_y = new TH1F("h_Event_AllHCal_y", "Event with jets in any HCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);


	h_Event_JetMeas_nHCal_0_Q2 = new TH1F("h_Event_JetMeas_nHCal_0_Q2", "Event with 0 jets measured in nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_JetMeas_nHCal_0_x = new TH1F("h_Event_JetMeas_nHCal_0_x", "Event with 0 jets measured in nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_JetMeas_nHCal_0_y = new TH1F("h_Event_JetMeas_nHCal_0_y", "Event with 0 jets measured in nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_JetMeas_nHCal_1_Q2 = new TH1F("h_Event_JetMeas_nHCal_1_Q2", "Event with 1 jet measured in nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_JetMeas_nHCal_1_x = new TH1F("h_Event_JetMeas_nHCal_1_x", "Event with 1 jet measured in nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_JetMeas_nHCal_1_y = new TH1F("h_Event_JetMeas_nHCal_1_y", "Event with 1 jet measured in nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_JetMeas_nHCal_2_Q2 = new TH1F("h_Event_JetMeas_nHCal_2_Q2", "Event with 2 jets measured in nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_JetMeas_nHCal_2_x = new TH1F("h_Event_JetMeas_nHCal_2_x", "Event with 2 jets measured in nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_JetMeas_nHCal_2_y = new TH1F("h_Event_JetMeas_nHCal_2_y", "Event with 2 jets measured in nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);

	h_Event_JetMeas_AllHCal_Q2 = new TH1F("h_Event_JetMeas_AllHCal_Q2", "Event with jets measured in any HCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_JetMeas_AllHCal_x = new TH1F("h_Event_JetMeas_AllHCal_x", "Event with jets measured in any HCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_JetMeas_AllHCal_y = new TH1F("h_Event_JetMeas_AllHCal_y", "Event with jets measured in any HCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);


	h_Event_JetMeas_no_nHCal_Q2 = new TH1F("h_Event_JetMeas_no_nHCal_Q2", "Event with jets measured without nHCal Q^{2}; Q^{2} [GeV^{2}/c^{2}]; counts", 1000, 0.0, 10.0);
	h_Event_JetMeas_no_nHCal_x = new TH1F("h_Event_JetMeas_no_nHCal_x", "Event with jets measured without nHCal x; x [1]; counts", nbins_x, logBinsArray_x);
	h_Event_JetMeas_no_nHCal_y = new TH1F("h_Event_JetMeas_no_nHCal_y", "Event with jets measured without nHCal inelasticity y; y [1]; counts", 1000, 0.0, 1.0);



	h_Event_HCal_jets = new TH2F("h_Event_HCal_jets", "Event with dijets in HCals; Jet #1; Jet #2; counts", 4, 0.0, 4.0, 4, 0.0, 4.0);
	h_Event_HCal_jets_meas = new TH2F("h_Event_HCal_jets_meas", "Event with measured dijets in HCals; Jet #1; Jet #2; counts", 4, 0.0, 4.0, 4, 0.0, 4.0);
	h_Event_HCal_jets_meas_no_nHCal = new TH2F("h_Event_HCal_jets_meas_no_nHCal", "Event with measured dijets in HCals without nHCal; Jet #1; Jet #2; counts", 4, 0.0, 4.0, 4, 0.0, 4.0);


	h_Event_HCal_jets->GetXaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets->GetXaxis()->SetBinLabel(2, "Barrel HCal");
	h_Event_HCal_jets->GetXaxis()->SetBinLabel(3, "LFHCAL");
	h_Event_HCal_jets->GetXaxis()->SetBinLabel(4, "Any");

	h_Event_HCal_jets->GetYaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets->GetYaxis()->SetBinLabel(2, "Barrel HCal");
	h_Event_HCal_jets->GetYaxis()->SetBinLabel(3, "LFHCAL");
	h_Event_HCal_jets->GetYaxis()->SetBinLabel(4, "Any");


	h_Event_HCal_jets_meas->GetXaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets_meas->GetXaxis()->SetBinLabel(2, "Barrel HCal");
	h_Event_HCal_jets_meas->GetXaxis()->SetBinLabel(3, "LFHCAL");
	h_Event_HCal_jets_meas->GetXaxis()->SetBinLabel(4, "Any");

	h_Event_HCal_jets_meas->GetYaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets_meas->GetYaxis()->SetBinLabel(2, "Barrel HCal");
	h_Event_HCal_jets_meas->GetYaxis()->SetBinLabel(3, "LFHCAL");
	h_Event_HCal_jets_meas->GetYaxis()->SetBinLabel(4, "Any");


	h_Event_HCal_jets_meas_no_nHCal->GetXaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets_meas_no_nHCal->GetXaxis()->SetBinLabel(2, "Barrel HCal");
	h_Event_HCal_jets_meas_no_nHCal->GetXaxis()->SetBinLabel(3, "LFHCAL");
	h_Event_HCal_jets_meas_no_nHCal->GetXaxis()->SetBinLabel(4, "Any");

	h_Event_HCal_jets_meas_no_nHCal->GetYaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets_meas_no_nHCal->GetYaxis()->SetBinLabel(2, "Barrel HCal");
	h_Event_HCal_jets_meas_no_nHCal->GetYaxis()->SetBinLabel(3, "LFHCAL");
	h_Event_HCal_jets_meas_no_nHCal->GetYaxis()->SetBinLabel(4, "Any");


	// full
	h_Event_HCal_jets_meas_full = new TH2F("h_Event_HCal_jets_meas_full", "Event with measured dijets in HCals; Jet #1; Jet #2; counts", 6, 0.0, 6.0, 6, 0.0, 6.0);

	h_Event_HCal_jets_meas_full->GetXaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets_meas_full->GetXaxis()->SetBinLabel(2, "nHCal+Barrel HCal");
	h_Event_HCal_jets_meas_full->GetXaxis()->SetBinLabel(3, "Barrel HCal");
	h_Event_HCal_jets_meas_full->GetXaxis()->SetBinLabel(4, "Barrel HCal+LFHCAL");
	h_Event_HCal_jets_meas_full->GetXaxis()->SetBinLabel(5, "LFHCAL");
	h_Event_HCal_jets_meas_full->GetXaxis()->SetBinLabel(6, "Any");

	h_Event_HCal_jets_meas_full->GetYaxis()->SetBinLabel(1, "nHCal");
	h_Event_HCal_jets_meas_full->GetYaxis()->SetBinLabel(2, "nHCal+Barrel HCal");
	h_Event_HCal_jets_meas_full->GetYaxis()->SetBinLabel(3, "Barrel HCal");
	h_Event_HCal_jets_meas_full->GetYaxis()->SetBinLabel(4, "Barrel HCal+LFHCAL");
	h_Event_HCal_jets_meas_full->GetYaxis()->SetBinLabel(5, "LFHCAL");
	h_Event_HCal_jets_meas_full->GetYaxis()->SetBinLabel(6, "Any");



	// Outgoing partons
	h_Partons_status = new TH1F("h_Partons_status", "Outgoing parton status; status; counts", 401, -200.5, 200.5);
	h_Partons_types = new TH1F("h_Partons_types", "Outgoing parton type; type; counts", 41, -0.5, 40.5);
	h_Partons_types_anti = new TH1F("h_Partons_types_anti", "Outgoing parton type; type; counts", 41, -0.5, 40.5);

	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(1), "d");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(2), "u");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(3), "s");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(4), "c");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(5), "b");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(6), "t");
	//h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(7), "b^{#prime}");
	//h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(8), "t^{#prime}");

	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(11), "e^{-}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(12), "#nu_{e}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(13), "\mu^{-}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(14), "#nu_{#mu}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(15), "\tau^{-}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(16), "#nu_{#tau}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(17), "#tau^{#prime-}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(18), "#nu_{#tau^{#prime}}");

	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(21), "g");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(22), "#gamma");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(23), "Z^{0}");
	h_Partons_types->GetXaxis()->SetBinLabel(h_Partons_types->GetXaxis()->FindBin(24), "W^{+}");


	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(1), "d");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(2), "u");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(3), "s");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(4), "c");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(5), "b");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(6), "t");
	//h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(7), "b^{#prime}");
	//h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(8), "t^{#prime}");

	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(11), "e^{-}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(12), "#nu_{e}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(13), "#mu^{-}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(14), "#nu_{#mu}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(15), "#tau^{-}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(16), "#nu_{#tau}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(17), "#tau^{#prime-}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(18), "#nu_{#tau^{#prime}}");

	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(21), "g");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(22), "#gamma");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(23), "Z^{0}");
	h_Partons_types_anti->GetXaxis()->SetBinLabel(h_Partons_types_anti->GetXaxis()->FindBin(24), "W^{+}");


	h_Partons_eta = new TH2F("h_Partons_eta", "Outgoing partons #eta_{1} vs #eta_{2}; #eta_{1} [1]; #eta_{2} [1]; counts", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_Partons_phi = new TH2F("h_Partons_phi", "Outgoing partons #phi_{1} vs #phi_{2}; #phi_{1} [1]; #phi_{2} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());

	h_Partons_p = new TH2F("h_Partons_p", "Outgoing partons p_{1} vs p_{2}; p_{1} [GeV/c]; p_{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Partons_pT = new TH2F("h_Partons_pT", "Outgoing partons p_{T}^{1} vs p_{T}^{2}; p_{T}^{1} [GeV/c]; p_{T}^{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);

	h_Parton_y_eta1 = new TH2F("h_Parton_y_eta1", "Outgoing partons #eta_{1} vs. y; y [1]; #eta_{1} [1]; counts", 100, 0.0, 1.0, 100, -5.0, 5.0);

	h_Parton_eta_p = new TH2F("h_Parton_eta_p", "Outgoing parton p vs #eta; #eta [1]; p [GeV/c]; counts", 100, -5.0, 5.0, 200, 0.0, 20.0);
	h_Parton_eta_pT = new TH2F("h_Parton_eta_pT", "Outgoing parton p_{T} vs #eta; #eta [1]; p_{T} [GeV/c]; counts", 100, -5.0, 5.0, 200, 0.0, 20.0);
	h_Parton_eta_E = new TH2F("h_Parton_eta_E", "Outgoing parton E vs #eta; #eta [1]; E [GeV]; counts", 100, -5.0, 5.0, 200, 0.0, 20.0);

	h_Parton_x_eta = new TH2F("h_Parton_x_eta", "Outgoing parton #eta vs. x; x [1]; #eta [1]; counts", nbins_x, logBinsArray_x, 100, -5.0, 5.0);
	h_Parton_y_eta = new TH2F("h_Parton_y_eta", "Outgoing partons #eta vs. y; y [1]; #eta [1]; counts", 100, 0.0, 1.0, 100, -5.0, 5.0);

	h_Parton_x_eta1 = new TH2F("h_Parton_x_eta1", "Outgoing partons #eta_{1} vs. x; x [1]; #eta_{1} [1]; counts", nbins_x, logBinsArray_x, 100, -5.0, 5.0);
	h_Parton_x_eta2 = new TH2F("h_Parton_x_eta2", "Outgoing partons #eta_{2} vs. x; x [1]; #eta_{2} [1]; counts", nbins_x, logBinsArray_x, 100, -5.0, 5.0);

	h_Parton_y_eta1 = new TH2F("h_Parton_y_eta1", "Outgoing partons #eta_{1} vs. y; y [1]; #eta_{1} [1]; counts", 100, 0.0, 1.0, 100, -5.0, 5.0);
	h_Parton_y_eta2 = new TH2F("h_Parton_y_eta2", "Outgoing partons #eta_{2} vs. y; y [1]; #eta_{2} [1]; counts", 100, 0.0, 1.0, 100, -5.0, 5.0);


	// Particles
	h_Particle_eta = new TH1F("h_Particle_eta", "MC particle #eta; #eta; counts", 200, -10.0, 10.0);

	h_Particle_eta_wE = new TH1F("h_Particle_eta_wE", "MC particle #eta, E weighed; #eta; counts", 200, -10.0, 10.0);

	h_Particle_eta_p = new TH2F("h_Particle_eta_p", "MC particles #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_eta_pT = new TH2F("h_Particle_eta_pT", "MC particles #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_eta_E = new TH2F("h_Particle_eta_E", "MC particles #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);


    // eta, momentum
	h_Particle_pion_p_eta_p = new TH2F("h_Particle_Pion_p_eta_p", "MC particles #pi^{+} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_pion_n_eta_p = new TH2F("h_Particle_Pion_n_eta_p", "MC particles #pi^{-} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Kaon_p_eta_p = new TH2F("h_Particle_Kaon_p_eta_p", "MC particles K^{+} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Kaon_n_eta_p = new TH2F("h_Particle_Kaon_n_eta_p", "MC particles K^{-} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_proton_p_eta_p = new TH2F("h_Particle_Proton_p_eta_p", "MC particles p^{+} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_proton_n_eta_p = new TH2F("h_Particle_Proton_n_eta_p", "MC particles p^{-} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Electron_p_eta_p = new TH2F("h_Particle_Electron_p_eta_p", "MC particles e^{+} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Electron_n_eta_p = new TH2F("h_Particle_Electron_n_eta_p", "MC particles e^{-} #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);

	h_Particle_Neutron_eta_p = new TH2F("h_Particle_Neutron_eta_p", "MC particles n #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Gamma_eta_p = new TH2F("h_Particle_Gamma_eta_p", "MC particles #gamman #eta vs. momentum; #eta [1]; p_{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);


    // eta, transverse momentum pT
	h_Particle_pion_p_eta_pT = new TH2F("h_Particle_Pion_p_eta_pT", "MC particles #pi^{+} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_pion_n_eta_pT = new TH2F("h_Particle_Pion_n_eta_pT", "MC particles #pi^{-} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Kaon_p_eta_pT = new TH2F("h_Particle_Kaon_p_eta_pT", "MC particles K^{+} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Kaon_n_eta_pT = new TH2F("h_Particle_Kaon_n_eta_pT", "MC particles K^{-} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_proton_p_eta_pT = new TH2F("h_Particle_Proton_p_eta_pT", "MC particles p^{+} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_proton_n_eta_pT = new TH2F("h_Particle_Proton_n_eta_pT", "MC particles p^{-} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Electron_p_eta_pT = new TH2F("h_Particle_Electron_p_eta_pT", "MC particles e^{+} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Electron_n_eta_pT = new TH2F("h_Particle_Electron_n_eta_pT", "MC particles e^{-} #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);

	h_Particle_Neutron_eta_pT = new TH2F("h_Particle_Neutron_eta_pT", "MC particles n #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Gamma_eta_pT = new TH2F("h_Particle_Gamma_eta_pT", "MC particles #gamman #eta vs. momentum; #eta [1]; p_{T}^{MC} [GeV/c]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);


	// eta, energy
	h_Particle_Pion_p_eta_E = new TH2F("h_Particle_Pion_p_eta_E", "MC particles #pi^{+} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Pion_n_eta_E = new TH2F("h_Particle_Pion_n_eta_E", "MC particles #pi^{-} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Kaon_p_eta_E = new TH2F("h_Particle_Kaon_p_eta_E", "MC particles K^{+} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Kaon_n_eta_E = new TH2F("h_Particle_Kaon_n_eta_E", "MC particles K^{-} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Proton_p_eta_E = new TH2F("h_Particle_Proton_p_eta_E", "MC particles p^{+} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Proton_n_eta_E = new TH2F("h_Particle_Proton_n_eta_E", "MC particles p^{-} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Electron_p_eta_E = new TH2F("h_Particle_Electron_p_eta_E", "MC particles e^{+} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Electron_n_eta_E = new TH2F("h_Particle_Electron_n_eta_E", "MC particles e^{-} #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);

	h_Particle_Neutron_eta_E = new TH2F("h_Particle_Neutron_eta_E", "MC particles n #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);
	h_Particle_Gamma_eta_E = new TH2F("h_Particle_Gamma_eta_E", "MC particles #gamma #eta vs. energy; #eta [1]; E_{MC} [GeV]; counts", 200, -10.0, 10.0, 500, 0.0, 50.0);

	// Jets
	h_Jet_nPart = new TH1D("h_Jet_nPart", "Jet number of particles; N_{part} [1]; counts", 201, -0.5, 200.5);
	h_Jet_mass = new TH1D("h_Jet_mass", "Jet mass; m [GeV/c^{2}]; counts", 2000, 0.0, 20.0);
    h_Jet_charge = new TH1D("h_Jet_charge", "Jet charge; q [1]; counts", 101, -50.5, 50.5);
	h_Jet_E = new TH1D("h_Jet_E", "Jet energy; E [GeV]; counts", 500, 0.0, 50.0);
	h_Jet_p = new TH1D("h_Jet_p", "Jet momentum; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_Jet_pT = new TH1D("h_Jet_pT", "Jet transverse momentum; p_{T} [GeV/c]; counts", 500, 0.0, 50.0);
	h_Jet_eta = new TH1D("h_Jet_eta", "Jet #eta; #eta [1]; counts", 200, -5.0, 5.0);
	h_Jet_deta = new TH1D("h_Jet_deta", "Jet #Delta#eta; #Delta#eta [1]; counts", 400, -10.0, 10.0);

	h_Jets_eta = new TH2F("h_Jets_eta", "Jets #eta_{1} vs. #eta_{2}; #eta_{1} [1]; #eta_{2} [1]; counts", 200, -5.0, 5.0, 200, -5.0, 5.0);
	h_Jets_phi = new TH2F("h_Jets_phi", "Jets #phi_{1} vs. #phi_{2}; #phi_{1} [1]; #phi_{2} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());
	h_Jets_p = new TH2F("h_Jets_p", "Jets p_{1} vs. p_{2}; p_{1} [GeV/c]; p_{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Jets_pT = new TH2F("h_Jets_pT", "Jets p_{T}^{1} vs. p_{T}^{2}; p_{T}^{1} [GeV/c]; p_{T}^{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Jets_E = new TH2F("h_Jets_E", "Jets E_{1} vs. E_{2}; E_{1} [GeV/c]; E_{2} [GeV]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);


	// Jets measured
	h_Jet_meas_nPart = new TH1D("h_Jet_meas_nPart", "Jet measured number of particles; N_{part} [1]; counts", 201, -0.5, 200.5);
	h_Jet_meas_mass = new TH1D("h_Jet_meas_mass", "Jet measured mass; m [GeV/c^{2}]; counts", 2000, 0.0, 20.0);
    h_Jet_meas_charge = new TH1D("h_Jet_meas_charge", "Jet measured charge; q [1]; counts", 101, -50.5, 50.5);
	h_Jet_meas_E = new TH1D("h_Jet_meas_E", "Jet measured energy; E [GeV]; counts", 500, 0.0, 50.0);
	h_Jet_meas_p = new TH1D("h_Jet_meas_p", "Jet measured momentum; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_Jet_meas_pT = new TH1D("h_Jet_meas_pT", "Jet measured transverse momentum; p_{T} [GeV/c]; counts", 500, 0.0, 50.0);
	h_Jet_meas_eta = new TH1D("h_Jet_meas_eta", "Jet measured #eta; #eta [1]; counts", 200, -5.0, 5.0);
	h_Jet_meas_deta = new TH1D("h_Jet_meas_deta", "Jet measured #Delta#eta; #Delta#eta [1]; counts", 400, -10.0, 10.0);

	h_Jets_meas_eta = new TH2F("h_Jets_meas_eta", "Jets measured #eta_{1} vs. #eta_{2}; #eta_{1} [1]; #eta_{2} [1]; counts", 200, -5.0, 5.0, 200, -5.0, 5.0);
	h_Jets_meas_phi = new TH2F("h_Jets_meas_phi", "Jets measured #phi_{1} vs. #phi_{2}; #phi_{1} [1]; #phi_{2} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());
	h_Jets_meas_p = new TH2F("h_Jets_meas_p", "Jets measured p_{1} vs. p_{2}; p_{1} [GeV/c]; p_{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Jets_meas_pT = new TH2F("h_Jets_meas_pT", "Jets measured p_{T}^{1} vs. p_{T}^{2}; p_{T}^{1} [GeV/c]; p_{T}^{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Jets_meas_E = new TH2F("h_Jets_meas_E", "Jets measured E_{1} vs. E_{2}; E_{1} [GeV/c]; E_{2} [GeV]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);


	h_Jet_bHCal_part_eta = new TH1D("h_Jet_bHCal_part_eta", "Jet in bHCal particle #eta; #eta [1]; counts", 200, -5.0, 5.0);
	h_Jet_meas_bHCal_part_eta = new TH1D("h_Jet_meas_bHCal_part_eta", "Jet measured in bHCal particle #eta; #eta [1]; counts", 200, -5.0, 5.0);

	h_Jet_HCal_part_eta = new TH2D("h_Jet_HCal_part_eta", "Jet in HCals particle #eta; #eta [1]; HCal; counts", 200, -5.0, 5.0, 4, 0.0, 4.0);
	h_Jet_meas_HCal_part_eta = new TH2D("h_Jet_meas_HCal_part_eta", "Jet measured in HCals particle #eta; #eta [1]; HCal; counts", 200, -5.0, 5.0, 4, 0.0, 4.0);

	h_Jet_HCal_part_eta->GetYaxis()->SetBinLabel(1, "no HCal");
	h_Jet_HCal_part_eta->GetYaxis()->SetBinLabel(2, "nHCal");
	h_Jet_HCal_part_eta->GetYaxis()->SetBinLabel(3, "bHCal");
	h_Jet_HCal_part_eta->GetYaxis()->SetBinLabel(4, "LFHCAL");

	h_Jet_meas_HCal_part_eta->GetYaxis()->SetBinLabel(1, "no HCal");
	h_Jet_meas_HCal_part_eta->GetYaxis()->SetBinLabel(2, "nHCal");
	h_Jet_meas_HCal_part_eta->GetYaxis()->SetBinLabel(3, "bHCal");
	h_Jet_meas_HCal_part_eta->GetYaxis()->SetBinLabel(4, "LFHCAL");



	// Jets measured without nHCal
	h_Jet_meas_no_nHCal_nPart = new TH1D("h_Jet_meas_no_nHCal_nPart", "Jet measured without nHCal number of particles; N_{part} [1]; counts", 201, -0.5, 200.5);
	h_Jet_meas_no_nHCal_mass = new TH1D("h_Jet_meas_no_nHCal_mass", "Jet measured without nHCal mass; m [GeV/c^{2}]; counts", 2000, 0.0, 20.0);
    h_Jet_meas_no_nHCal_charge = new TH1D("h_Jet_meas_no_nHCal_charge", "Jet measured without nHCal charge; q [1]; counts", 101, -50.5, 50.5);
	h_Jet_meas_no_nHCal_E = new TH1D("h_Jet_meas_no_nHCal_E", "Jet measured without nHCal energy; E [GeV]; counts", 500, 0.0, 50.0);
	h_Jet_meas_no_nHCal_p = new TH1D("h_Jet_meas_no_nHCal_p", "Jet measured without nHCal momentum; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_Jet_meas_no_nHCal_pT = new TH1D("h_Jet_meas_no_nHCal_pT", "Jet measured without nHCal transverse momentum; p_{T} [GeV/c]; counts", 500, 0.0, 50.0);
	h_Jet_meas_no_nHCal_eta = new TH1D("h_Jet_meas_no_nHCal_eta", "Jet measured without nHCal #eta; #eta [1]; counts", 200, -5.0, 5.0);
	h_Jet_meas_no_nHCal_deta = new TH1D("h_Jet_meas_no_nHCal_deta", "Jet measured without nHCal #Delta#eta; #Delta#eta [1]; counts", 400, -10.0, 10.0);

	h_Jets_meas_no_nHCal_eta = new TH2F("h_Jets_meas_no_nHCal_eta", "Jets measured without nHCal #eta_{1} vs. #eta_{2}; #eta_{1} [1]; #eta_{2} [1]; counts", 200, -5.0, 5.0, 200, -5.0, 5.0);
	h_Jets_meas_no_nHCal_phi = new TH2F("h_Jets_meas_no_nHCal_phi", "Jets measured without nHCal #phi_{1} vs. #phi_{2}; #phi_{1} [1]; #phi_{2} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());
	h_Jets_meas_no_nHCal_p = new TH2F("h_Jets_meas_no_nHCal_p", "Jets measured without nHCal p_{1} vs. p_{2}; p_{1} [GeV/c]; p_{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Jets_meas_no_nHCal_pT = new TH2F("h_Jets_meas_no_nHCal_pT", "Jets measured without nHCal p_{T}^{1} vs. p_{T}^{2}; p_{T}^{1} [GeV/c]; p_{T}^{2} [GeV/c]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);
	h_Jets_meas_no_nHCal_E = new TH2F("h_Jets_meas_no_nHCal_E", "Jets measured without nHCal E_{1} vs. E_{2}; E_{1} [GeV/c]; E_{2} [GeV]; counts", 200, 0.0, 20.0, 200, 0.0, 20.0);


	// measured jets vs. partons
	h_Jets_meas_Partons_eta = new TH2F("h_Jets_meas_Partons_eta", "Jet #eta vs parton #eta; #eta_{parton} [1]; #eta_{jet} [1]; counts", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_Jets_meas_Partons_phi = new TH2F("h_Jets_meas_Partons_phi", "Jet #phi vs parton #phi; #phi_{parton} [1]; #phi_{jet} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());
	h_Jets_meas_Partons_E = new TH2F("h_Jets_meas_Partons_E", "Jet E vs parton E; E_{parton} [GeV]; E_{jet} [GeV]; counts", 100, 0.0, 50.0, 100, 0.0, 50.0);

	h_Jet_meas_Parton_eta1 = new TH2F("h_Jet_meas_Parton_eta1", "Jet #eta vs parton 1 #eta; #eta_{parton} [1]; #eta_{jet} [1]; counts", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_Jet_meas_Parton_phi1 = new TH2F("h_Jet_meas_Parton_phi1", "Jet #phi vs parton 1 #phi; #phi_{parton} [1]; #phi_{jet} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());
	h_Jet_meas_Parton_E1 = new TH2F("h_Jet_meas_Parton_E1", "Jet E vs parton 1 E; E_{parton} [GeV]; E_{jet} [GeV]; counts", 100, 0.0, 50.0, 100, 0.0, 50.0);

	h_Jet_meas_Parton_eta2 = new TH2F("h_Jet_meas_Parton_eta2", "Jet #eta vs parton 2 #eta; #eta_{parton} [1]; #eta_{jet} [1]; counts", 100, -5.0, 5.0, 100, -5.0, 5.0);
	h_Jet_meas_Parton_phi2 = new TH2F("h_Jet_meas_Parton_phi2", "Jet #phi vs parton 2 #phi; #phi_{parton} [1]; #phi_{jet} [1]; counts", 200, -2.0*TMath::Pi(), 2.0*TMath::Pi(), 200, -2.0*TMath::Pi(), 2.0*TMath::Pi());
	h_Jet_meas_Parton_E2 = new TH2F("h_Jet_meas_Parton_E2", "Jet E vs parton 2 E; E_{parton} [GeV]; E_{jet} [GeV]; counts", 100, 0.0, 50.0, 100, 0.0, 50.0);


	// temp
	const int nEtaBins = 7;
	double EtaBins[nEtaBins+1] = {-10, -4.14, -1.18, -1.1, 1.1, 1.2, 4.2, 10};

	hist_eta_energy_tmp = new TH1F("hist_eta_energy_tmp", "hist_eta_energy_tmp; #eta, E [GeV]", nEtaBins, EtaBins);
	hist_eta_energy_denom_tmp = new TH1F("hist_eta_energy_denom_tmp", "hist_eta_energy_denom_tmp; #eta, E [GeV]", nEtaBins, EtaBins);


	//hist_eta_energy_tmp->Sumw2();
	//hist_eta_energy_denom_tmp->Sumw2();


	return 1;
}


#endif /* HISTOGRAMSPYTHIA_H_ */
