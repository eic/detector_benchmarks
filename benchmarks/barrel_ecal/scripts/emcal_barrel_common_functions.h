//////////////////////////
// Common Particle Functions
// M. Scott 05/2021
//////////////////////////

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