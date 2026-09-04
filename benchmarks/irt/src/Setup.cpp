#include "Setup.h"
Setup::Setup(){
  std::cout<<"IRT Setup Called\n"<<std::endl;
  ///PFRICH (mom & eta)
  m_nBinsMomPFRICH = 150;
  m_PLowerPFRICH = 0.;
  m_PUpperPFRICH = 15.0;
  m_nBinsEtaPFRICH = 50;
  m_ELowerPFRICH = -1.50;
  m_EUpperPFRICH = -3.50;
  // (Angle and Hits)
  m_nBinsThetaPFRICH = 150;
  m_ThetaLowerPFRICH = 150;
  m_ThetaUpperPFRICH = 300;
  m_nBinsNpePFRICH = 50;
  m_NpeLowerPFRICH = 0;
  m_NpeUpperPFRICH = 50;
  ///DRICH (mom & eta)
  m_nBinsMomDRICH = 500;
  m_PLowerDRICH = 0.;
  m_PUpperDRICH = 50.;
  m_nBinsEtaDRICH = 50;
  m_ELowerDRICH = 1.50;
  m_EUpperDRICH = 3.50;
  // (Angle and Hits)
  m_nBinsThetaDRICH = 100;
  m_ThetaLowerGasDRICH =0.;
  m_ThetaUpperGasDRICH =50;
  m_ThetaLowerAerogelDRICH =200.;
  m_ThetaUpperAerogelDRICH =250.;
  m_nBinsNpeDRICH = 50;
  m_NpeLowerDRICH = 0;
  m_NpeUpperDRICH = 50;
  
}
void Setup::HistoInit(){
  if(m_IrtRICH == "PFRICH"){
    std::cout<<"PFRICH Histos for: "<<m_IrtRICH<<std::endl;
    //TH1::AddDirectory(false);
    //std::unique_ptr<TH1D> m_IrtPDG;
    //std::unique_ptr<TH1D> m_IrtNpe;
    //std::unique_ptr<TH1D> m_IrtHits;
    
    //Tracks
    m_histTitle = m_IrtRICH+"_IRT_Mom; P(GeV/c)";
    m_RICH_Mom = std::make_unique<TH1D>("m_RICH_Mom",m_histTitle.c_str(),
					m_nBinsMomPFRICH, m_PLowerPFRICH,
					m_PUpperPFRICH
					);
    
    m_histTitle = m_IrtRICH+"_IRT_Eta; #eta";
    m_RICH_Eta = std::make_unique<TH1D>("m_RICH_ETA",m_histTitle.c_str(),
					m_nBinsEtaPFRICH, m_ELowerPFRICH,
					m_EUpperPFRICH
					);
    //Radiators
    m_histTitle = m_IrtRICH+"_Npe_Rad_0; #theta(mrad)";
    m_Rad_Theta.push_back(std::make_unique<TH1D>("m_Theta_Rad_0",m_histTitle.c_str(),
					       m_nBinsThetaPFRICH, m_ThetaLowerPFRICH,
					       m_ThetaUpperPFRICH
					       ));
    m_histTitle = m_IrtRICH+"_Npe_Rad_0; Npe";
    m_Rad_Npe.push_back(std::make_unique<TH1D>("m_Npe_Rad_0",m_histTitle.c_str(),
					       m_nBinsNpePFRICH, m_NpeLowerPFRICH,
					       m_NpeUpperPFRICH
					       ));
    m_histTitle = m_IrtRICH+"_NHits_Rad_0; NHits";
    m_Rad_NHits.push_back(std::make_unique<TH1D>("m_NHits_Rad_0",m_histTitle.c_str(),
						 m_nBinsNpePFRICH, m_NpeLowerPFRICH,
						 m_NpeUpperPFRICH));
    //Radiators-Track Combination
    m_histTitle = m_IrtRICH+"_ThetaMom_Rad_0; #theta(mrad)";
    m_Theta_Mom.push_back(std::make_unique<TH2D>("m_ThetaMom_Rad_0",m_histTitle.c_str(),
						 m_nBinsMomPFRICH, m_PLowerPFRICH,
						 m_PUpperPFRICH,
						 m_nBinsThetaPFRICH, m_ThetaLowerPFRICH,
						 m_ThetaUpperPFRICH
					       ));
    
    //std::unique_ptr<TFile> m_OutFile
    m_OutRootFile = std::make_unique<TFile>((m_OutFile+"_"+m_IrtRICH+".edm4hep.root").c_str(),
					    "RECREATE");
    std::cout 
    << "File: " << m_OutRootFile->GetName()
    << " writable: " << m_OutRootFile->IsWritable()
    << " zombie: " << m_OutRootFile->IsZombie()
    << std::endl;  
  }

  
  if(m_IrtRICH == "DRICH"){
    std::cout<<"DRICH Histos for: "<<m_IrtRICH<<std::endl;
    //Tracks
    m_histTitle = m_IrtRICH+"_IRT_Mom; P(GeV/c)";
    m_RICH_Mom = std::make_unique<TH1D>("m_RICH_Mom",m_histTitle.c_str(),
					m_nBinsMomDRICH, m_PLowerDRICH,
					m_PUpperDRICH
					);
    
    m_histTitle = m_IrtRICH+"_IRT_Eta; #eta";
    m_RICH_Eta = std::make_unique<TH1D>("m_RICH_ETA",m_histTitle.c_str(),
					m_nBinsEtaDRICH, m_ELowerDRICH,
					m_EUpperDRICH
					);
    //Radiators
    //Theta
    m_histTitle = m_IrtRICH+"_Theta_Rad_0; #theta(mrad)";
    m_Rad_Theta.push_back(std::make_unique<TH1D>("m_Theta_Rad_0",m_histTitle.c_str(),
						 m_nBinsThetaDRICH, m_ThetaLowerAerogelDRICH,
						 m_ThetaUpperAerogelDRICH
						 ));
    m_histTitle = m_IrtRICH+"_Theta_Rad_1; #theta(mrad)";
    m_Rad_Theta.push_back(std::make_unique<TH1D>("m_Theta_Rad_1",m_histTitle.c_str(),
						 m_nBinsThetaDRICH, m_ThetaLowerGasDRICH,
						 m_ThetaUpperGasDRICH
						 ));
    //Npe
    m_histTitle = m_IrtRICH+"_Npe_Rad_0; Npe";
    m_Rad_Npe.push_back(std::make_unique<TH1D>("m_Npe_Rad_0",m_histTitle.c_str(),
					       m_nBinsNpeDRICH, m_NpeLowerDRICH,
					       m_NpeUpperDRICH
					       ));
    m_histTitle = m_IrtRICH+"_Npe_Rad_1; Npe";
    m_Rad_Npe.push_back(std::make_unique<TH1D>("m_Npe_Rad_1",m_histTitle.c_str(),
					       m_nBinsNpeDRICH, m_NpeLowerDRICH,
					       m_NpeUpperDRICH
					       ));
    //Hits
    m_histTitle = m_IrtRICH+"_NHits_Rad_0; NHits";
    m_Rad_NHits.push_back(std::make_unique<TH1D>("m_NHits_Rad_0",m_histTitle.c_str(),
						 m_nBinsNpeDRICH, m_NpeLowerDRICH,
						 m_NpeUpperDRICH));
    m_histTitle = m_IrtRICH+"_NHits_Rad_1; NHits";
    m_Rad_NHits.push_back(std::make_unique<TH1D>("m_NHits_Rad_1",m_histTitle.c_str(),
						 m_nBinsNpeDRICH, m_NpeLowerDRICH,
						 m_NpeUpperDRICH));
    
    
   
    //Radiators-Track Combination
    m_histTitle = m_IrtRICH+"_ThetaMom_Rad_0; #theta(mrad)";
    m_Theta_Mom.push_back(std::make_unique<TH2D>("m_ThetaMom_Rad_0",m_histTitle.c_str(),
						 m_nBinsMomDRICH, m_PLowerDRICH,
						 m_PUpperDRICH,
						 m_nBinsThetaDRICH, m_ThetaLowerAerogelDRICH,
						 m_ThetaUpperAerogelDRICH
						 ));
    m_histTitle = m_IrtRICH+"_ThetaMom_Rad_1; #theta(mrad)";
    m_Theta_Mom.push_back(std::make_unique<TH2D>("m_ThetaMom_Rad_1",m_histTitle.c_str(),
						 m_nBinsMomDRICH, m_PLowerDRICH,
						 m_PUpperDRICH,
						 m_nBinsThetaDRICH, m_ThetaLowerGasDRICH,
						 m_ThetaUpperGasDRICH
						 ));
  
    
    //std::unique_ptr<TFile> m_OutFile
    m_OutRootFile = std::make_unique<TFile>((m_OutFile+"_"+m_IrtRICH+".edm4hep.root").c_str(),
					    "RECREATE");
    std::cout 
      << "File: " << m_OutRootFile->GetName()
      << " writable: " << m_OutRootFile->IsWritable()
      << " zombie: " << m_OutRootFile->IsZombie()
      << std::endl;  
  }
}
Setup::~Setup(){
  std::cout<<"IRT Setup Destroyed\n"<<std::endl;
}
