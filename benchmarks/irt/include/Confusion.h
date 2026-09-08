#pragma once
#include <memory>
#include <TH2D.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

class Confusion {
 public:
  Confusion();
  ~Confusion();
 protected:
  std::unique_ptr<TH2D> m_ConfusionMatrix;
  enum class PID {
    Electron,
      Muon,
      Pion,
      Kaon,
      Proton,
      Count
      };
  void BookConfusionMatrix();
  int PDGToBin(int pdg);
  void FillConfusion(int,int);
};
