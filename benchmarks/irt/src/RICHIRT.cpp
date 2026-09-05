#include "RICHIRT.h"
#include <stdexcept>
///////////////////////////
RICHIRT::RICHIRT(){
  std::cout<<"RICH IRT Called\n"<<std::endl;
}
///////////////////////////
void RICHIRT::SetCollections(){
  m_IrtRadiatorInfo = get_irtRICH()+"IrtRadiatorInfo";
  m_IrtParticleCollection = get_irtRICH()+"IrtParticles";
  m_Tracks =  get_irtRICH()+"Tracks";
  m_MCParticles = "MCParticles";
  m_CentralCKFTracks="CentralCKFTracks";
}
/////////////////////
void RICHIRT::Process(){
  auto reader = podio::makeReader(get_inFile());
  std::cout << "Events in file: " << reader.getEvents() << std::endl;

  std::vector<int> mc_truth;
  for (size_t evt = 0; evt <reader.getEvents(); ++evt) {
    auto frame = reader.readEvent(evt);
    //std::cout << "\n=== Event " << evt << " ===" << std::endl;

    const auto& mcParts = frame.get<edm4hep::MCParticleCollection>(m_MCParticles);
    // IRT Radiators, Particles and Track Segment
    const auto& irtRads = frame.get<edm4eic::IrtRadiatorInfoCollection>(m_IrtRadiatorInfo);
    const auto& irtParts = frame.get<edm4eic::IrtParticleCollection>(m_IrtParticleCollection);
    const auto& tracks = frame.get<edm4eic::TrackSegmentCollection>(m_Tracks);
    //Central track
    const auto& cCKFTracks=frame.get<edm4eic::TrackCollection>(m_CentralCKFTracks);
    //Dump MC Part PDG
    for(auto mcpart: mcParts){
      if(mcpart.getGeneratorStatus()!=1) continue;
      mc_truth.push_back(PDGToBin(mcpart.getPDG()));
    }
    //Get IRT-Tracks
    TVector3 p3;
    for(auto trs: tracks){
      //std::cout<<"DEBUG TRACK SEGMENT:: "<<trs.getLength()<<std::endl;
      auto tr = trs.getTrack();
      if(tr.isAvailable()){
        auto p = tr.getMomentum();
        p3 = TVector3 (p.x,p.y,p.z);
	m_RICH_Mom->Fill(p3.Mag());
	m_RICH_Eta->Fill(p3.PseudoRapidity());
        //std::cout<<"DEBUG: "<<p3.Mag()<<" "<<p3.PseudoRapidity()<<std::endl;
      }else std::cerr<<"DEBUG No Track "<<std::endl;
    }//IRT-Tracks
    int idx =0;
    //IRT Particles
    for(auto part: irtParts){
      int recPDG = part.getPDG();
      int recNp = part.getNpe();
      //if(recNp <=2 )std::cout<<recPDG<<" "<<recNp<<std::endl;
      //if(recPDG!= 211 ) std::cout<<"Info: "<<recPDG<<" "<<p3.Mag()<<" "<<p3.PseudoRapidity()<<std::endl;
      int rec = PDGToBin(part.getPDG());
      int truth = mc_truth.at(idx); idx++;
      FillConfusion(rec,truth);
     
      ///// OneToMany : Radiator
      //std::cout<<m_Rad_Npe.size()<<"; "<<part.getRadiators().size()<<std::endl;
      if(m_Rad_Npe.size()!= part.getRadiators().size()){
	throw std::runtime_error("Radiator size (" + std::to_string(part.getRadiators().size())
				 + ") and histogram assignment (" + std::to_string(m_Rad_Npe.size())
				 + ") mismatch");
      }
      for (const auto& rad : part.getRadiators()) {
        auto id = rad.id();
        //if(part.getPDG() == 0) continue;
        m_Rad_Theta[id.index]->Fill(rad.getAngle());
        m_Rad_Npe[id.index]->Fill(rad.getNpe());
        m_Rad_NHits[id.index]->Fill(rad.getNhits());
	m_Theta_Mom[id.index]->Fill(p3.Mag(),rad.getAngle());
        //CKF_Angle->Fill(momCKF.Mag(),rad.getAngle());
        //MC_Angle->Fill(momMCP.Mag(),rad.getAngle());
      }//Radiators
    }//Irt-Particles
  }//events
  BookHistos();
}
void RICHIRT::BookHistos(){
  m_OutRootFile->cd();
  m_RICH_Mom->Write();
  m_RICH_Eta->Write();
  for(auto& h: m_Rad_Theta) h->Write();
  for(auto& h: m_Rad_Npe) h->Write();
  for(auto& h: m_Rad_NHits) h->Write();
  for(auto& h: m_Theta_Mom) h->Write();
  BookConfusionMatrix();
  m_OutRootFile->Write();
  m_OutRootFile->Close();

}
RICHIRT::~RICHIRT(){
   std::cout<<"RICH IRT Destroyed\n"<<std::endl;
}
