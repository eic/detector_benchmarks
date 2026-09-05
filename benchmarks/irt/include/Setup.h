#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
class Setup{
public:
  Setup();
  ~Setup();
  void HistoInit();
  //void Setup(const char* irtRICH, const char* inFile, const char* outFile);
  // Setters
  void set_irtRICH(const char* det)   { m_IrtRICH =  det;    }
  void set_inFile (const char* input) { m_InFile  =  input;  }
  void set_outFile(const char* output){ m_OutFile =  output; }
  // Getters
  const std::string& get_irtRICH() const { return m_IrtRICH; }
  const std::string& get_inFile()  const { return m_InFile;  }
  const std::string& get_outFile() const { return m_OutFile; }
  
private:
  std::string m_IrtRICH;
  std::string m_InFile;
  std::string m_OutFile;
protected:
  //////Histogramming
  //Vars
  std::string m_histName;
  std::string m_histTitle;
  ///PFRICH (mom & eta)
  double m_PLowerPFRICH;
  double m_PUpperPFRICH;
  double m_ELowerPFRICH;
  double m_EUpperPFRICH;
  int m_nBinsMomPFRICH;
  int m_nBinsEtaPFRICH;
  // (Angle and Hits)
  int m_nBinsThetaPFRICH;
  double m_ThetaLowerPFRICH;
  double m_ThetaUpperPFRICH;
  double m_nBinsNpePFRICH;
  double m_NpeLowerPFRICH;
  double m_NpeUpperPFRICH;
  ///DRICH (mom & eta)
  int m_nBinsMomDRICH;
  double m_PLowerDRICH;
  double m_PUpperDRICH;
  double m_nBinsEtaDRICH;
  double m_ELowerDRICH;
  double m_EUpperDRICH;
  // (Angle and Hits)
  int m_nBinsThetaDRICH;
  double m_ThetaLowerGasDRICH;
  double m_ThetaUpperGasDRICH;
  double m_ThetaLowerAerogelDRICH;
  double m_ThetaUpperAerogelDRICH;
  int m_nBinsNpeDRICH;
  double m_NpeLowerDRICH;
  double m_NpeUpperDRICH;
  //TH1Ds IRT-Parts
  std::unique_ptr<TH1D> m_IrtPDG;
  std::unique_ptr<TH1D> m_IrtNpe;
  std::unique_ptr<TH1D> m_IrtHits;

  //Tracks
  std::unique_ptr<TH1D> m_RICH_Mom;
  std::unique_ptr<TH1D> m_RICH_Eta;
  
  //Radiator (Vector of NRadaitors)
  std::vector<std::unique_ptr<TH1D>> m_Rad_Theta;
  std::vector<std::unique_ptr<TH1D>> m_Rad_Npe;
  std::vector<std::unique_ptr<TH1D>> m_Rad_NHits;
  //Radiator-Track
  std::vector<std::unique_ptr<TH2D>> m_Theta_Mom;
  
  //OutFile;
  std::unique_ptr<TFile> m_OutRootFile;
};
