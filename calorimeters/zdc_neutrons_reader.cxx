#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TH1F.h"
#include <iostream>
using namespace HepMC3;

void zdc_neutrons_reader(){

  //-------------------------------------
  ReaderAscii hepmc_input("data/neutrons_zdc.hepmc");
  int        events_parsed = 0;
  GenEvent   evt(Units::GEV, Units::MM);

  TH1F* h_neutron_energy = new TH1F("n energy","; E [GeV]",100,0,200);

  while(!hepmc_input.failed()) {
    // Read event from input file
    hepmc_input.read_event(evt);

    // If reading failed - exit loop
    if( hepmc_input.failed() ) break;

    for(const auto& v : evt.vertices() ) {
      for(const auto& p : v->particles_out() ) {
        if(p->pid() == 2112) {
          h_neutron_energy->Fill(p->momentum().e());
        }
      }
    }
    evt.clear();
    events_parsed++;
  }
  std::cout << "Events parsed and written: " << events_parsed << std::endl;

  TCanvas* c = new TCanvas();
  h_neutron_energy->Draw();
  c->SaveAs("results/zdc_neutrons_reader.png");
  c->SaveAs("results/zdc_neutrons_reader.pdf");
}


