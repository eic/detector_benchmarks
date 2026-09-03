#ifndef COMMON_BENCH_PARTICLES_HH
#define COMMON_BENCH_PARTICLES_HH 1

#include <map>

namespace common_bench {

  struct ParticleData {
    int pdgCode;
    int charge;
    double mass;
    // std::string name;
  };

  using ParticleMap = std::map<int, ParticleData>;

  const ParticleMap particleMap = {
      {          11, {          11,  -1,   0.000510998928 }},  // e-
      {         -11, {         -11,   1,   0.000510998928 }},  // e+
      {          13, {          13,  -1,   0.105658357    }},  // mu-
      {         -13, {         -13,   1,   0.105658357    }},  // mu+
      {          22, {          22,   0,   0.0            }},  // gamma
      {         111, {         111,   0,   0.1349766      }},  // pi0
      {         211, {         211,   1,   0.1395701      }},  // pi+
      {        -211, {        -211,  -1,   0.1395701      }},  // pi-
      {         130, {         130,   0,   0.49767        }},  // K_L^0
      {         310, {         310,   0,   0.49767        }},  // K_S^0
      {         321, {         321,   1,   0.49360        }},  // K^+
      {        -321, {        -321,  -1,   0.49360        }},  // K^-
      {        2212, {        2212,   1,   0.93827        }},  // p+
      {       -2212, {       -2212,  -1,   0.93827        }},  // p~-
      {        2112, {        2112,   0,   0.93957        }},  // n0
      {  1000010020, {  1000010020,   1,   1.87561        }},  // 1000010020 Deuterium
      {  1000010030, {  1000010030,   1,   2.80925        }},  // 1000010030 Tritium
      {  1000020030, {  1000020030,   1,   2.80923        }},  // 1000020030 He3
      {  1000020040, {  1000020040,   1,   3.72742        }},  // 1000020040 Alpha-(He4)
  };

} // namespace common_bench

#endif
