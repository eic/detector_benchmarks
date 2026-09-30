#include "Confusion.h"
Confusion::Confusion(){
  constexpr int N = 5;
  m_ConfusionMatrix = std::make_unique<TH2D>("Confusion",
					     "RICH PID Confusion Matrix;Reco PID;True PID",
					     N, 0, N,
					     N, 0, N
					     );

}
Confusion::~Confusion(){
  std::cout<<"Confusion Destructed"<<std::endl;
}
void Confusion::BookConfusionMatrix(){ 
  constexpr int N = 5;
  const char* labels[] = {"e", "#pi", "K", "p", "#mu"};

  for (int i = 0; i < N; ++i) {
    m_ConfusionMatrix->GetXaxis()->SetBinLabel(i + 1, labels[i]);
    m_ConfusionMatrix->GetYaxis()->SetBinLabel(i + 1, labels[i]);
  }
  m_ConfusionMatrix->Write();
}

void Confusion::FillConfusion(int reco, int truth) {
  m_ConfusionMatrix->Fill(reco,truth);
}
int Confusion::PDGToBin(int pdg) {
  switch (std::abs(pdg)) { 
  case 11:   return 0;
  case 211:  return 1;
  case 321:  return 2;
  case 2212: return 3;
  case 13:   return 4;
  default:   return -1;
  }
}
