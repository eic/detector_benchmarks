#pragma once

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using RVecI       = ROOT::VecOps::RVec<int>;

// Lazy execution methods
std::pair<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>> createReconstructionPlots( RNode d1 ){

  int ePrimID = 4;
  int pBeamID = 1;
  int eBeamID = 0;
  
  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;

  d1 = d1.Define("eBeam", [eBeamID]
                            (std::vector<edm4hep::MCParticleData> h)
                            {return ROOT::Math::PxPyPzMVector (h[eBeamID].momentum.x,h[eBeamID].momentum.y,h[eBeamID].momentum.z,h[eBeamID].mass);},
                            {"MCParticles"})
        .Define("pBeam", [pBeamID]
                            (std::vector<edm4hep::MCParticleData> h)
                            {return ROOT::Math::PxPyPzMVector (h[pBeamID].momentum.x,h[pBeamID].momentum.y,h[pBeamID].momentum.z,h[pBeamID].mass);},
                            {"MCParticles"})
        .Define("primMom", [ePrimID]
                            (std::vector<edm4hep::MCParticleData> h)
                            {return ROOT::Math::PxPyPzMVector (h[ePrimID].momentum.x,h[ePrimID].momentum.y,h[ePrimID].momentum.z,h[ePrimID].mass);},
                            {"MCParticles"})
        .Define("primE","primMom.E()")
        .Define("primPx","primMom.Px()")
        .Define("primPy","primMom.Py()")
        .Define("primPz","primMom.Pz()")
        .Define("primTheta","1000*(TMath::Pi()-primMom.Theta())")
        .Define("primEta","primMom.Eta()")
        .Define("scatteredV","eBeam-primMom")
        .Define("Q2","-scatteredV.M2()")
        .Define("logQ2","log10(Q2)")
        .Define("x","Q2/(2*scatteredV.Dot(pBeam))")
        .Define("logx","log10(x)")
        .Define("primPhi","primMom.Phi()")
        .Define("W2","(pBeam+eBeam-primMom).M2()")
        .Define("reconTheta","(double)LowQ2TrackParameters[0].theta")
        .Define("reconPhi","(double)LowQ2TrackParameters[0].phi")
        .Define("reconP","(double)(LowQ2TrackParameters[0].qOverP*LowQ2TrackParameters[0].charge)")
        .Define("reconMom",[](double p, double theta, double phi){return ROOT::Math::PxPyPzMVector(p*sin(theta)*cos(phi),p*sin(theta)*sin(phi),p*cos(theta),0.000511);},{"reconP","reconTheta","reconPhi"})
        .Define("reconE","reconMom.E()")
        .Define("reconScatteredV","eBeam-reconMom")
        .Define("reconQ2","-reconScatteredV.M2()")
        .Define("reconlogQ2","log10(reconQ2)")
        .Define("reconX","reconQ2/(2*reconScatteredV.Dot(pBeam))")
        .Define("reconlogx","log10(reconX)")
        .Define("reconW2","(pBeam+eBeam-reconMom).M2()")
        .Define("reconPx","(double)reconMom.Px()")
        .Define("reconPy","(double)reconMom.Py()")
        .Define("reconPz","(double)reconMom.Pz()")
        .Define("thetaRes","reconTheta-primTheta")
        .Define("phiRes","reconPhi-primPhi")
        .Define("ERes","reconMom.E()-primE")
        .Define("Q2Res","reconQ2-Q2")
        .Define("XRes","reconX-x")
        .Define("W2Res","reconW2-W2");
        
        

  hHists2D["reconThetaVsPrimTheta"] = d1.Histo2D({"reconThetaVsPrimTheta","reconThetaVsPrimTheta;#theta_{prim} [mrad];#theta_{recon} [mrad]",100,0,100,100,0,100},"primTheta","reconTheta");
  
  hHists2D["reconPhiVsPrimPhi"] = d1.Histo2D({"reconPhiVsPrimPhi","reconPhiVsPrimPhi;#phi_{prim} [rad];#phi_{recon} [rad]",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()},"primPhi","reconPhi");

  hHists2D["reconPxVsPrimPx"] = d1.Histo2D({"reconPxVsPrimPx","reconPxVsPrimPx;p_{x,prim} [GeV];p_{x,recon} [GeV]",100,-0.5,0.5,100,-0.5,0.5},"primPx","reconPx");
  hHists2D["reconPyVsPrimPy"] = d1.Histo2D({"reconPyVsPrimPy","reconPyVsPrimPy;p_{y,prim} [GeV];p_{y,recon} [GeV]",100,-0.5,0.5,100,-0.5,0.5},"primPy","reconPy");
  hHists2D["reconPzVsPrimPz"] = d1.Histo2D({"reconPzVsPrimPz","reconPzVsPrimPz;p_{z,prim} [GeV];p_{z,recon} [GeV]",100,0,0.5,100,0,0.5},"primPz","reconPz");
  hHists2D["reconEVsPrimE"]   = d1.Histo2D({"reconEVsPrimE","reconEVsPrimE;E_{prim} [GeV];E_{recon} [GeV]",100,0,0.5,100,0,0.5},"primE","reconE");
  hHists2D["reconQ2VsPrimQ2"] = d1.Histo2D({"reconQ2VsPrimQ2","reconQ2VsPrimQ2;Q^{2}_{prim} [GeV^{2}];Q^{2}_{recon} [GeV^{2}]",100,0,0.5,100,0,0.5},"Q2","reconQ2");
  hHists2D["reconXVsPrimX"]   = d1.Histo2D({"reconXVsPrimX","reconXVsPrimX;x_{prim};x_{recon}",100,0,0.5,100,0,0.5},"x","reconX");
  hHists2D["reconW2VsPrimW2"] = d1.Histo2D({"reconW2VsPrimW2","reconW2VsPrimW2;W^{2}_{prim} [GeV^{2}];W^{2}_{recon} [GeV^{2}]",100,0,0.5,100,0,0.5},"W2","reconW2");

  hHists1D["thetaRes"] = d1.Histo1D({"thetaRes","thetaRes;#theta_{recon}-#theta_{prim} [mrad]",100,-50,50},"thetaRes");
  hHists1D["phiRes"]   = d1.Histo1D({"phiRes","phiRes;#phi_{recon}-#phi_{prim} [rad]",100,-0.1,0.1},"phiRes");
  hHists1D["ERes"]     = d1.Histo1D({"ERes","ERes;E_{recon}-E_{prim} [GeV]",100,-0.1,0.1},"ERes");
  hHists1D["Q2Res"]    = d1.Histo1D({"Q2Res","Q2Res;Q^{2}_{recon}-Q^{2}_{prim} [GeV^{2}]",100,-0.1,0.1},"Q2Res");
  hHists1D["XRes"]     = d1.Histo1D({"XRes","XRes;x_{recon}-x_{prim}",100,-0.1,0.1},"XRes");
  hHists1D["W2Res"]    = d1.Histo1D({"W2Res","W2Res;W^{2}_{recon}-W^{2}_{prim} [GeV^{2}]",100,-0.1,0.1},"W2Res");

  hHists1D["thetaResVsE"] = d1.Histo1D({"thetaResVsE","thetaResVsE;E_{prim} [GeV];#theta_{recon}-#theta_{prim} [mrad]",100,0,0.5},"primE","thetaRes");
  hHists1D["phiResVsE"]   = d1.Histo1D({"phiResVsE","phiResVsE;E_{prim} [GeV];#phi_{recon}-#phi_{prim} [rad]",100,0,0.5},"primE","phiRes");
  hHists1D["EResVsE"]     = d1.Histo1D({"EResVsE","EResVsE;E_{prim} [GeV];E_{recon}-E_{prim} [GeV]",100,0,0.5},"primE","ERes");
  hHists1D["Q2ResVsE"]    = d1.Histo1D({"Q2ResVsE","Q2ResVsE;E_{prim} [GeV];Q^{2}_{recon}-Q^{2}_{prim} [GeV^{2}]",100,0,0.5},"primE","Q2Res");

  //Plot showing where the phi resolution is less than 30 degrees in terms of E and theta
  //hHists2D["phiResVsETheta"] = d1.Histo2D({"phiResVsETheta","phiResVsETheta;E_{prim} [GeV];#theta_{prim} [mrad]",100,0,0.5,100,0,100},"primE","primTheta",[](double phiRes){return fabs(phiRes)<0.5;},{"phiRes-primPhi"});

  return {hHists1D,hHists2D};
 
}
