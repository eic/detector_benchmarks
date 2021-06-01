//////////////////////////
// Common Particle Functions
// M. Scott 05/2021
//////////////////////////

#include <cctype>
#include "TH1.h"

// Returns particle pdgID and mass in [GeV]
std::tuple <int, double> extract_particle_parameters(std::string particle_name) {
    if (particle_name == "electron") return std::make_tuple(11,    0.51099895e-3);
    if (particle_name == "photon")   return std::make_tuple(22,    0.0);
    if (particle_name == "positron") return std::make_tuple(-11,   0.51099895e-3);
    if (particle_name == "proton")   return std::make_tuple(2212,  0.938272);
    if (particle_name == "muon")     return std::make_tuple(13,    0.1056583745);
    if (particle_name == "pi0")      return std::make_tuple(111,   0.1349768);
    if (particle_name == "piplus")   return std::make_tuple(211,   0.13957039);
    if (particle_name == "piminus")  return std::make_tuple(-211,  0.13957039);

    std::cout << "wrong particle name" << std::endl;
    abort();
}

// Returns Environment Variables
std::string getEnvVar( std::string const & key ){
    char * val = getenv( key.c_str() );
    return val == NULL ? std::string("") : std::string(val);
}
// Added detetcor name to title
void addDetectorName(std::string name, TH1 * inhist){
  std::string newName = inhist -> GetTitle();
  for (auto& x : name){
    x = toupper(x);
  }
  inhist -> SetTitle((name + " " + newName).c_str());
}