#ifndef MT_H
#define MT_H

// Defines the number of threads to run within the ROOT analysis scripts.
// TODO: make this a file configured by the CI scripts so we can specify
//       the number of threads (and the number of processes) at a global
//       level

constexpr const int kNumThreads = 8;

#endif
