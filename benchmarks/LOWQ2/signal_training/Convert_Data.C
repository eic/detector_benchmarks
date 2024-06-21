#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>
#include <TChain.h>


#include "MCParticle.hpp"
#include "PixelCharge.hpp"
#include "PixelHit.hpp"


int Convert_Data(TString inName="output/Out_tpx4_genprop.root", TString outName="output/Out_Convert_genprop.root") {



    // Open the input file
    TFile* inputFile = new TFile(inName, "READ");

    TTree* MCTree          = static_cast<TTree*>(inputFile->Get("MCParticle"));
    TTree* PixelChargeTree = static_cast<TTree*>(inputFile->Get("PixelCharge"));
    TTree* PixelHitTree    = static_cast<TTree*>(inputFile->Get("PixelHit")); 

    // get branch from the chain
    TBranch* mcParticleBranch  = MCTree->GetBranch("mydetector");
    TBranch* pixelChargeBranch = PixelChargeTree->GetBranch("mydetector");
    TBranch* pixelHitBranch    = PixelHitTree->GetBranch("mydetector");

    // Bind the information to predefined vectors
    std::vector<allpix::PixelHit*> input_hits;
    pixelHitBranch->SetObject(&input_hits);

    std::vector<allpix::PixelCharge*> input_charges;
    pixelChargeBranch->SetObject(&input_charges);

    std::vector<allpix::MCParticle*> input_particles;
    mcParticleBranch->SetObject(&input_particles);

    // Create output tree and branches
    TFile* outputFile = new TFile(outName, "RECREATE");
    TTree* outputTree = new TTree("events", "events");

    // Constants
    double pitch = 0.055;

    // Limit grid size from 10x10
    int grid_size = 6;
    int grid_area = grid_size*grid_size;
    // shift grid by half the size rounded down
    int shift_grid = int((10-grid_size)/2);

    // Create branches in the output tree for the members of the allpix::MCParticle class
    double x, y, z = 0;
    double px, py, pz = 0;
    outputTree->Branch("x",    &x   );
    outputTree->Branch("y",    &y   );
    outputTree->Branch("z",    &z   );
    outputTree->Branch("px",   &px  );
    outputTree->Branch("py",   &py  );
    outputTree->Branch("pz",   &pz  );
    // outputTree->Branch("time", &time);

    // Create branches in the output tree for the members of the allpix::PixelHit class
    std::vector<int>    pixel_x, pixel_y;
    std::vector<double> charge,  time;
    outputTree->Branch("pixel_x", &pixel_x);
    outputTree->Branch("pixel_y", &pixel_y);
    outputTree->Branch("charge",  &charge );
    outputTree->Branch("time",    &time   );


    //loop over the entries in the input branches filling the output tree
    Long64_t nentries = mcParticleBranch->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        mcParticleBranch->GetEntry(i);
        pixelChargeBranch->GetEntry(i);
        pixelHitBranch->GetEntry(i);

        charge = std::vector<double>(grid_area, 0);
        time   = std::vector<double>(grid_area, 0);

        for (auto particle: input_particles) {
            auto start = particle->getGlobalStartPoint();
            auto end   = particle->getGlobalEndPoint();
            x    = start.X();
            y    = start.Y();
            z    = start.Z();

            auto x_end = end.X();
            auto y_end = end.Y();
            auto z_end = end.Z();

            px = x_end-x;
            py = y_end-y;
            pz = z_end-z;

            //Convert into unit vector components
            double norm = sqrt(px*px + py*py + pz*pz);
            px = px/norm;
            py = py/norm;
            pz = pz/norm;

            x = x/pitch;
            y = y/pitch;

            break;
        }

        for (auto hit: input_hits) {
            if(hit->getSignal()==0) continue;
            int x_index = hit->getIndex().X();
            int y_index = hit->getIndex().Y();

            //shift grid by half the size rounded down
            x_index = x_index - shift_grid;
            y_index = y_index - shift_grid;

            if(x_index < 0 || x_index >= grid_size || y_index < 0 || y_index >= grid_size) continue;
            
            //pixel_x.push_back(x_index);
            //pixel_y.push_back(y_index);
            //charge.push_back(hit->getSignal());
            //time.push_back(hit->getLocalTime());
            charge[x_index*grid_size + y_index] = hit->getSignal();
            time[x_index*grid_size + y_index] = hit->getLocalTime();
        }

        outputTree->Fill();

        //reset vectors
        pixel_x.clear();
        pixel_y.clear();
        charge.clear();
        time.clear();

    }

    // Write the output tree to the output file
    outputFile->Write();
    outputFile->Close();

    inputFile->Close();

    return 0;
}