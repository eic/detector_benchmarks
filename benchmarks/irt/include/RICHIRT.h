#pragma once 
#include <podio/Reader.h>
#include <podio/Frame.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4eic/IrtRadiatorInfoCollection.h>
#include <edm4eic/IrtParticleCollection.h>
#include <edm4eic/TrackSegmentCollection.h>
#include <edm4eic/TrackCollection.h>
#include <TVector3.h>

#include "Setup.h"
#include "Confusion.h"
class RICHIRT : public Setup, public Confusion{
 public:
  RICHIRT();
  ~RICHIRT();
  void SetCollections();
  void Process();
  void BookHistos();
 private:
  std::string m_IrtRadiatorInfo;
  std::string m_IrtParticleCollection;
  std::string m_Tracks;
  std::string m_CentralCKFTracks;
  std::string m_MCParticles;
  
};
