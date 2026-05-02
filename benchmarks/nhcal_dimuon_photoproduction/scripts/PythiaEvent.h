/*
 * PythiaEvent.h
 *
 *  Created on: 24 lis 2020
 *      Author: Khaless
 */

#ifndef PythiaEvent_H
#define PythiaEvent_H

#include <vector>

// ROOT
#include "TObject.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

using namespace Pythia8;

class PythiaEvent : public TObject {
public:
	PythiaEvent() {}
  //PythiaEvent(const PythiaEvent& ev);
  ~PythiaEvent()	{}

  // event
  long eventId;
  int nParticlesFinal;
  float Q2; // momentum transfer
  float x; // momentum fraction
  float y; // inelasticity

  Particle scatteredEle;

  vector<Particle> particles;


/*
  // Upsilon candidate
  float p;
  float pt;
  float y;
  float eta;
  float phi;
  float m;
  short charge;

  // electron 1
  float p1;
  float pt1;
  float eta1;
  float phi1;
  float nSigmaElectron1;
  float nSigmaPion1;

  // electron 2
  float p2;
  float pt2;
  float eta2;
  float phi2;
  float nSigmaElectron2;
  float nSigmaPion2;*/

private:
  ClassDef(PythiaEvent,1);
};



#endif /* PythiaEvent_H */


