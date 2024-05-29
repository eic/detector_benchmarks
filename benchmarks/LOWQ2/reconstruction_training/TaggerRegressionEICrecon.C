#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
 
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
 
#include "TMVA/MethodDNN.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"
  
using namespace TMVA;
 
// The training currently requires files which only contained a single DIS scattered electron to have been simulated e.g. generated using GETaLM
// The scattered electron must be the 3rd particle in the file after the two beam particles
// At least one track reconstructed by EIC algorithms in the LOWQ2 tagger is needed.

void TaggerRegressionEICrecon(
			      TString inDataNames    = "/scratch/EIC/ReconOut/qr_18x275_ab/qr_18x275_ab*_recon.edm4hep.root",
			      TString outDataName    = "/scratch/EIC/LowQ2Model/trainedData.root",
            TString dataFolderName = "LowQ2Model",
            TString mcBeamEnergy   = "18",
            TString typeName       = "LowQ2MomentumRegression",
            TString methodName     = "DNN",
			      TString inWeightName   = "dataset/weights/LowQ2Reconstruction_DNN.weights.xml"
			      )
{
 
  Bool_t loadWeights = 0;
  
  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();
  
  ROOT::EnableImplicitMT(8);
  
  // --------------------------------------------------------------------------------------------------
  // Here the preparation phase begins
  // Create a new root output file
  TFile* outputFile = TFile::Open( outDataName, "RECREATE" );
 
  // Create the factory object. Later you can choose the methods
   
  TMVA::Factory *factory = new TMVA::Factory( typeName, outputFile,
					      "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );

  ; 
  TMVA::DataLoader *dataloader=new TMVA::DataLoader(dataFolderName);
        
  // Input TrackParameters variables from EICrecon - 
  TString collectionName     = "TaggerTrackerProjectedTracks[0]";
  dataloader->AddVariable( collectionName+".loc.a", "fit_position_y", "units", 'F' );
  dataloader->AddVariable( collectionName+".loc.b", "fit_position_z", "units", 'F' );
  dataloader->AddVariable( "sin("+collectionName+".phi)*sin("+collectionName+".theta)",   "fit_vector_x",   "units", 'F' );
  dataloader->AddVariable( "cos("+collectionName+".phi)*sin("+collectionName+".theta)",   "fit_vector_y",   "units", 'F' );
  
  // Regression target particle 3-momentum, normalised to beam energy.
  // Takes second particle, in the test data this is the scattered electron
  // TODO add energy and array element information to be read directly from datafile - EMD4eic and EICrecon changes.
  TString mcParticleName = "MCParticles[MCScatteredElectrons_objIdx[0].index]";
  //TString mcParticleName = "MCParticles[0]";
  dataloader->AddTarget( mcParticleName+".momentum.x/"+mcBeamEnergy );
  dataloader->AddTarget( mcParticleName+".momentum.y/"+mcBeamEnergy );
  dataloader->AddTarget( mcParticleName+".momentum.z/"+mcBeamEnergy );
 
  std::cout << "--- TMVARegression           : Using input files: " << inDataNames << std::endl;
 
  // Register the regression tree 
  TChain* regChain = new TChain("events");
  regChain->Add(inDataNames);
  //regChain->SetEntries(8000); // Set smaller sample for tests
  
  // global event weights per tree (see below for setting event-wise weights)
  Double_t regWeight  = 1.0;
 
  // You can add an arbitrary number of regression trees
  dataloader->AddRegressionTree( regChain, regWeight );
 
  // This would set individual event weights (the variables defined in the
  // expression need to exist in the original TTree)
  // dataloader->SetWeightExpression( "1/(eE)", "Regression" ); // If MC event weights are kept use these
  // Apply additional cuts on the data
  TCut mycut = "@TaggerTrackerProjectedTracks.size()==1"; // Make sure there's one reconstructed track in event
   
  dataloader->PrepareTrainingAndTestTree(mycut,"nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:SplitSeed=1:NormMode=NumEvents:!V");

  // TODO - Optimise layout and training more
  TString layoutString("Layout=TANH|1024,TANH|128,TANH|64,TANH|32,LINEAR");
  
  TString trainingStrategyString("TrainingStrategy=");
  trainingStrategyString +="LearningRate=1e-4,Momentum=0,MaxEpochs=2000,ConvergenceSteps=200,BatchSize=64,TestRepetitions=1,Regularization=None,Optimizer=ADAM";   
  
  TString nnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:WeightInitialization=XAVIERUNIFORM:RandomSeed=1234");

  // Use GPU if possible on the machine
  TString architectureString("Architecture=GPU");

  // Transformation of data prior to training layers - decorrelate and normalise whole dataset
  TString transformString("VarTransform=D,N");
  
  nnOptions.Append(":");
  nnOptions.Append(architectureString);
  nnOptions.Append(":");
  nnOptions.Append(transformString);
  nnOptions.Append(":");
  nnOptions.Append(layoutString);
  nnOptions.Append(":");
  nnOptions.Append(trainingStrategyString);
  
  TMVA::MethodDNN* method = (MethodDNN*)factory->BookMethod(dataloader, TMVA::Types::kDL, methodName, nnOptions); // NN
  
  // If loading previous model for further training
  if(loadWeights){
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    reader->BookMVA( methodName, inWeightName );
    TMVA::MethodDNN* kl = dynamic_cast<TMVA::MethodDNN*>(reader->FindMVA(methodName));
    method = kl;
  }


  // --------------------------------------------------------------------------------------------------
  // Now you can tell the factory to train, test, and evaluate the MVAs
  
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  
  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  // --------------------------------------------------------------
  
  // Save the output
  outputFile->Close();
  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVARegression is done!" << std::endl;
  
  // delete factory;
  // delete dataloader;
  
}
