//-------------------------
//
// Hit reader to relate hits at Roman Pots to momentum vectors from MC.
//
// Input(s): output file from npsim particle gun for RP particles.
//
// Output(s): txt file with training information with px_mc, py_mc, pz_mc, x_rp, slope_xrp, y_rp, slope_yrp
//
//
// Author: Alex Jentsch
//------------------------
//Low PT preprocessing added by David Ruth

using namespace std;

void preprocess_model_training_data(TString inputFile, TString outputFile, TString outputFile_lo) {

    string fileName;
    TFile* inputRootFile;
    TTree* rootTree;

    // MC information
    TH1D* h_eta_MC = new TH1D("h_eta", ";Pseudorapidity, #eta", 100, 0.0, 15.0);
    TH1D* h_px_MC = new TH1D("px_MC", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
    TH1D* h_py_MC = new TH1D("py_MC", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
    TH1D* h_pt_MC = new TH1D("pt_MC", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
    TH1D* h_pz_MC = new TH1D("pz_MC", ";p_{z} [GeV/c]", 100, 0.0, 320.0);
    TH1D* h_theta_MC = new TH1D("theta_MC", ";#theta [mrad]", 100, 0.0, 25.0);
    TH1D* h_phi_MC = new TH1D("phi_MC", ";#phi [rad]", 100, -3.2, 3.2);

    // Roman pots
    TH1D* h_px_RomanPots = new TH1D("px_RomanPots", ";p_{x} [GeV/c]", 100, -10.0, 10.0);
    TH1D* h_py_RomanPots = new TH1D("py_RomanPots", ";p_{y} [GeV/c]", 100, -10.0, 10.0);
    TH1D* h_pt_RomanPots = new TH1D("pt_RomanPots", ";p_{t} [GeV/c]", 100, 0.0, 2.0);
    TH1D* h_pz_RomanPots = new TH1D("pz_RomanPots", ";p_{z} [GeV/c]", 100, 0.0, 320.0);

    int fileCounter = 0;
    int iEvent = 0;

    inputRootFile = new TFile(inputFile);
    if (!inputRootFile) { 
        cout << "MISSING_ROOT_FILE" << fileName << endl; 
        return;
    }

    TTree* evtTree = (TTree*)inputRootFile->Get("events");

    int numEvents = evtTree->GetEntries();

    TTreeReader tree_reader(evtTree); // !the tree reader

    ofstream outputTrainingFile;
    outputTrainingFile.open(outputFile);

    ofstream outputTrainingFile_lo;
    outputTrainingFile_lo.open(outputFile_lo);
    // MC particles
    TTreeReaderArray<float> mc_px_array = { tree_reader, "MCParticles.momentum.x" };
    TTreeReaderArray<float> mc_py_array = { tree_reader, "MCParticles.momentum.y" };
    TTreeReaderArray<float> mc_pz_array = { tree_reader, "MCParticles.momentum.z" };
    TTreeReaderArray<double> mc_mass_array = { tree_reader, "MCParticles.mass" };
    TTreeReaderArray<int> mc_pdg_array = { tree_reader, "MCParticles.PDG" };
    TTreeReaderArray<int> mc_genStatus_array = { tree_reader, "MCParticles.generatorStatus" };

    // Roman pots -- momentum vector

    // hit locations (for debugging)
    TTreeReaderArray<double> global_hit_RP_x = { tree_reader, "ForwardRomanPotHits.position.x" };
    TTreeReaderArray<double> global_hit_RP_y = { tree_reader, "ForwardRomanPotHits.position.y" };
    TTreeReaderArray<double> global_hit_RP_z = { tree_reader, "ForwardRomanPotHits.position.z" };

    cout << "file has " << evtTree->GetEntries() << " events..." << endl;

    tree_reader.SetEntriesRange(0, evtTree->GetEntries());
    while (tree_reader.Next()) {

        cout << "Reading event: " << iEvent << endl;

        if (mc_px_array.GetSize() > 1) { 
            cout << "Event has more than one particle -- skip..." << endl;
        }

        bool hasMCProton = false;
        bool hitLayerOne = false;
        bool hitLayerTwo = false;
	bool lowPT = false;

        double mcProtonMomentum[3];
        TVector3 mctrk;

        for (int imc = 0; imc < mc_px_array.GetSize(); imc++) {
            mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);

            if (mc_pdg_array[imc] == 2212 && mc_genStatus_array[imc] == 1) { //only checking for protons here -- change as desired

                mcProtonMomentum[0] = mctrk.Px();
                mcProtonMomentum[1] = mctrk.Py();
                mcProtonMomentum[2] = mctrk.Pz();
                hasMCProton = true;
            }
        }

        double hit1minZ = 25099.0;
        double hit1maxZ = 26022.0;
        double hit2minZ = 27099.0;
        double hit2maxZ = 28022.0;

        double rpHitLayerOne[3];
        double rpHitLayerTwo[3];

        double slopeXRP = 0.0;
        double slopeYRP = 0.0;

        // roman pots reco tracks
        for (int iRPPart = 0; iRPPart < global_hit_RP_x.GetSize(); iRPPart++) {

            if (global_hit_RP_z[iRPPart] > hit1minZ && global_hit_RP_z[iRPPart] < hit1maxZ) {

                rpHitLayerOne[0] = global_hit_RP_x[iRPPart];
                rpHitLayerOne[1] = global_hit_RP_y[iRPPart];
                rpHitLayerOne[2] = global_hit_RP_z[iRPPart];
                hitLayerOne = true;
            }

            if (global_hit_RP_z[iRPPart] > hit2minZ && global_hit_RP_z[iRPPart] < hit2maxZ) {

                rpHitLayerTwo[0] = global_hit_RP_x[iRPPart];
                rpHitLayerTwo[1] = global_hit_RP_y[iRPPart];
                rpHitLayerTwo[2] = global_hit_RP_z[iRPPart];
                hitLayerTwo = true;
            }
        }
        if (hasMCProton) {
		double Pt = sqrt(pow(mcProtonMomentum[0],2) + pow(mcProtonMomentum[1],2));
		if(Pt < 0.3) {
			lowPT = true;
		}
	}
        if (hasMCProton && hitLayerOne && hitLayerTwo) {

            double slope_x = (rpHitLayerTwo[0] - rpHitLayerOne[0]) / (rpHitLayerTwo[2] - rpHitLayerOne[2]);
            double slope_y = (rpHitLayerTwo[1] - rpHitLayerOne[1]) / (rpHitLayerTwo[2] - rpHitLayerOne[2]);

            outputTrainingFile << mcProtonMomentum[0] << "\t" << mcProtonMomentum[1] << "\t" << mcProtonMomentum[2] << "\t";
            outputTrainingFile << rpHitLayerTwo[0] << "\t" << slope_x << "\t" << rpHitLayerTwo[1] << "\t" << slope_y << endl;
	    if (lowPT) {
               outputTrainingFile_lo << mcProtonMomentum[0] << "\t" << mcProtonMomentum[1] << "\t" << mcProtonMomentum[2] << "\t";
               outputTrainingFile_lo << rpHitLayerTwo[0] << "\t" << slope_x << "\t" << rpHitLayerTwo[1] << "\t" << slope_y << endl;
            }
        }

        iEvent++;
    } // event loop

    inputRootFile->Close();
    outputTrainingFile.close();
    outputTrainingFile_lo.close();

    return;
}

