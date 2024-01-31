#pragma once

// Define alias
using RNode       = ROOT::RDF::RNode;
using H1ResultPtr = ROOT::RDF::RResultPtr<TH1D>;
using H2ResultPtr = ROOT::RDF::RResultPtr<TH2D>;
using H3ResultPtr = ROOT::RDF::RResultPtr<TH3D>;
using RVecI       = ROOT::VecOps::RVec<int>;

// Lazy execution methods
std::tuple<std::map<TString,H1ResultPtr>,std::map<TString,H2ResultPtr>,std::map<TString,H3ResultPtr>> createReconstructionPlots( RNode d1 ){

  int ePrimID = 2;
  int pBeamID = 1;
  int eBeamID = 0;
  
  std::map<TString,H1ResultPtr> hHists1D;
  std::map<TString,H2ResultPtr> hHists2D;
  std::map<TString,H3ResultPtr> hHists3D;

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
        .Define("primTheta","(TMath::Pi()-primMom.Theta())")
        .Define("primTheta_mrad","1000*primTheta")
        .Define("primEta","primMom.Eta()")
        .Define("scatteredV","eBeam-primMom")
        .Define("Q2","-scatteredV.M2()")
        .Define("logQ2","log10(Q2)")
        .Define("x","Q2/(2*scatteredV.Dot(pBeam))")
        .Define("logx","log10(x)")
        .Define("primPhi","primMom.Phi()")
        .Define("primPhi_deg","primPhi * TMath::RadToDeg()")
        .Define("W2","(pBeam+eBeam-primMom).M2()")
        .Define("reconTheta","(double)(TMath::Pi()-LowQ2TrackParameters[0].theta)")
        .Define("reconTheta_mrad","1000.0*reconTheta")
        .Define("reconPhi","(double)LowQ2TrackParameters[0].phi")
        .Define("reconPhi_deg","reconPhi * TMath::RadToDeg()")
        .Define("reconP","(18./10.)*(double)(1/(LowQ2TrackParameters[0].qOverP*LowQ2TrackParameters[0].charge))")
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
        .Define("thetaRes","(reconTheta-primTheta)")
        .Define("thetaRes_mrad","1000.0*thetaRes")
        .Define("phiRes","reconPhi_deg-primPhi_deg")
        .Define("ERes","(reconMom.E()-primE)/primE")
        .Define("logQ2Res","reconlogQ2-logQ2")
        .Define("XRes","reconX-x")
        .Define("W2Res","reconW2-W2");
        
      
    
  double thetaMin = 0; // mrad
  double thetaMax = 12; // mrad
  double phiMin = -180; // deg
  double phiMax = 180; // deg
  double pxMin = -0.2; // GeV
  double pxMax = 0.2; // GeV
  double pyMin = -0.2; // GeV
  double pyMax = 0.2; // GeV
  double pzMin = -18; // GeV
  double pzMax = 0; // GeV
  double eMin = 0; // GeV
  double eMax = 18; // GeV
  double logQ2Min = -8.5; // GeV^2
  double logQ2Max = 0; // GeV^2
  double logxMin = -11;
  double logxMax = 0;
  double w2Min = 0; // GeV^2
  double w2Max = 0.5; // GeV^2
  double thetaResMin = -5; // mrad
  double thetaResMax = 5; // mrad
  double phiResMin = -45; // deg
  double phiResMax = 45; // deg
  double eResMin = -0.1; // GeV
  double eResMax = 0.1; // GeV
  double q2ResMin = -0.01; // GeV^2
  double q2ResMax = 0.01; // GeV^2
  double xResMin = -0.1;
  double xResMax = 0.1;
  double w2ResMin = -0.1; // GeV^2
  double w2ResMax = 0.1; // GeV^2

  int bins1D    = 400;
  int bins2D    = 200;
  int bins3D    = 50;
  int bins3DRes = 100;

  hHists2D["reconThetaVsPrimTheta"] = d1.Histo2D({"reconThetaVsPrimTheta","reconThetaVsPrimTheta;#theta_{prim} [mrad];#theta_{recon} [mrad]",bins2D,thetaMin,thetaMax,bins2D,thetaMin,thetaMax},"primTheta_mrad","reconTheta_mrad");
  hHists2D["reconPhiVsPrimPhi"] = d1.Histo2D({"reconPhiVsPrimPhi","reconPhiVsPrimPhi;#phi_{prim} [deg];#phi_{recon} [deg]",bins2D,phiMin,phiMax,bins2D,phiMin,phiMax},"primPhi_deg","reconPhi_deg");

  hHists2D["reconPxVsPrimPx"] = d1.Histo2D({"reconPxVsPrimPx","reconPxVsPrimPx;p_{x,prim} [GeV];p_{x,recon} [GeV]",bins2D,pxMin,pxMax,bins2D,pxMin,pxMax},"primPx","reconPx");
  hHists2D["reconPyVsPrimPy"] = d1.Histo2D({"reconPyVsPrimPy","reconPyVsPrimPy;p_{y,prim} [GeV];p_{y,recon} [GeV]",bins2D,pyMin,pyMax,bins2D,pyMin,pyMax},"primPy","reconPy");
  hHists2D["reconPzVsPrimPz"] = d1.Histo2D({"reconPzVsPrimPz","reconPzVsPrimPz;p_{z,prim} [GeV];p_{z,recon} [GeV]",bins2D,pzMin,pzMax,bins2D,pzMin,pzMax},"primPz","reconPz");
  hHists2D["reconEVsPrimE"]   = d1.Histo2D({"reconEVsPrimE","reconEVsPrimE;E_{prim} [GeV];E_{recon} [GeV]",bins2D,eMin,eMax,bins2D,eMin,eMax},"primE","reconE");
  hHists2D["reconQ2VsPrimQ2"] = d1.Histo2D({"reconQ2VsPrimQ2","reconQ2VsPrimQ2;log(Q^{2}_{prim}) [GeV^{2}];log(Q^{2}_{recon}) [GeV^{2}]",bins2D,logQ2Min,logQ2Max,bins2D,logQ2Min,logQ2Max},"logQ2","reconlogQ2");
  hHists2D["reconXVsPrimX"]   = d1.Histo2D({"reconXVsPrimX","reconXVsPrimX;log(x_{prim});log(x_{recon})",bins2D,logxMin,logxMax,bins2D,logxMin,logxMax},"logx","reconlogx");
  hHists2D["reconW2VsPrimW2"] = d1.Histo2D({"reconW2VsPrimW2","reconW2VsPrimW2;W^{2}_{prim} [GeV^{2}];W^{2}_{recon} [GeV^{2}]",bins2D,w2Min,w2Max,bins2D,w2Min,w2Max},"W2","reconW2");

  hHists1D["thetaRes"] = d1.Histo1D({"thetaRes","thetaRes;#theta_{recon}-#theta_{prim} [mrad]",bins1D,thetaResMin,thetaResMax},"thetaRes_mrad");
  hHists1D["phiRes"]   = d1.Histo1D({"phiRes","phiRes;#phi_{recon}-#phi_{prim} [rad]",bins1D,phiResMin,phiResMax},"phiRes");
  hHists1D["ERes"]     = d1.Histo1D({"ERes","ERes;E_{recon}-E_{prim} [GeV]",bins1D,eResMin,eResMax},"ERes");
  hHists1D["Q2Res"]    = d1.Histo1D({"Q2Res","Q2Res;log(Q^{2}_{recon})-log(Q^{2}_{prim}) [GeV^{2}]",bins1D,q2ResMin,q2ResMax},"logQ2Res");
  hHists1D["XRes"]     = d1.Histo1D({"XRes","XRes;log(x_{recon})-log(x_{prim})",bins1D,xResMin,xResMax},"XRes");
  hHists1D["W2Res"]    = d1.Histo1D({"W2Res","W2Res;W^{2}_{recon}-W^{2}_{prim} [GeV^{2}]",bins1D,w2ResMin,w2ResMax},"W2Res");

  hHists2D["thetaResVsE"] = d1.Histo2D({"thetaResVsE", "thetaResVsE;#theta_{recon}-#theta_{prim} [mrad];E_{prim} [GeV]", bins2D, thetaResMin, thetaResMax, bins2D, eMin, eMax}, "thetaRes_mrad", "primE");
  hHists2D["phiResVsE"]   = d1.Histo2D({"phiResVsE", "phiResVsE;#phi_{recon}-#phi_{prim} [rad];E_{prim} [GeV]", bins2D, phiResMin, phiResMax, bins2D, eMin, eMax}, "phiRes", "primE");
  hHists2D["EResVsE"]     = d1.Histo2D({"EResVsE", "EResVsE;E_{recon}-E_{prim} [GeV];E_{prim} [GeV]", bins2D, eResMin, eResMax, bins2D, eMin, eMax}, "ERes", "primE");
  hHists2D["Q2ResVsE"]    = d1.Histo2D({"Q2ResVsE", "Q2ResVsE;log(Q^{2}_{recon})-log(Q^{2}_{prim}) [GeV^{2}];E_{prim} [GeV]", bins2D, q2ResMin, q2ResMax, bins2D, eMin, eMax}, "logQ2Res", "primE");

  hHists2D["phiResVsTheta"] = d1.Histo2D({"phiResVsTheta", "phiResVsTheta;#phi_{recon}-#phi_{prim} [rad];#theta_{prim} [mrad]", bins2D, phiResMin, phiResMax, bins2D, thetaMin, thetaMax}, "phiRes", "primTheta_mrad");

  //3D histograms to extract the resolution as a function of E and theta
  hHists3D["thetaResVsETheta"] = d1.Histo3D({"thetaResVsETheta","thetaResVsETheta;E_{prim} [GeV];#theta_{prim} [mrad];#theta_{recon}-#theta_{prim} [mrad]",bins3D,eMin,eMax,bins3D,thetaMin,thetaMax,bins3DRes,thetaResMin,thetaResMax},"primE","primTheta_mrad","thetaRes_mrad");
  hHists3D["phiResVsETheta"]   = d1.Histo3D({"phiResVsETheta","phiResVsETheta;E_{prim} [GeV];#theta_{prim} [mrad];#phi_{recon}-#phi_{prim} [rad]",bins3D,eMin,eMax,bins3D,thetaMin,thetaMax,bins3DRes,phiResMin,phiResMax},"primE","primTheta_mrad","phiRes");
  hHists3D["EResVsETheta"]     = d1.Histo3D({"EResVsETheta","EResVsETheta;E_{prim} [GeV];#theta_{prim} [mrad];E_{recon}-E_{prim} [GeV]",bins3D,eMin,eMax,bins3D,thetaMin,thetaMax,bins3DRes,eResMin,eResMax},"primE","primTheta_mrad","ERes");

  auto d2 = d1.Filter("reconTheta_mrad<2");

  hHists2D["thetacut/phiResVsE"]     = d2.Histo2D({"phiResVsE", "phiResVsE;#phi_{recon}-#phi_{prim} [rad];E_{prim} [GeV]", bins2D, phiResMin, phiResMax, bins2D, eMin, eMax}, "phiRes", "primE");
  hHists2D["thetacut/reconPhiVsPrimPhi"] = d2.Histo2D({"reconPhiVsPrimPhi","reconPhiVsPrimPhi;#phi_{prim} [deg];#phi_{recon} [deg]",bins2D,phiMin,phiMax,bins2D,phiMin,phiMax},"primPhi_deg","reconPhi_deg");

  //Plot showing where the phi resolution is less than 30 degrees in terms of E and theta
  //hHists2D["phiResVsETheta"] = d1.Histo2D({"phiResVsETheta","phiResVsETheta;E_{prim} [GeV];#theta_{prim} [mrad]",100,eMin,eMax,100,thetaMin,thetaMax},"primE","primTheta",[](double phiRes){return fabs(phiRes)<0.5;},{"phiRes-primPhi"});

  return {hHists1D, hHists2D, hHists3D};  
  
}

