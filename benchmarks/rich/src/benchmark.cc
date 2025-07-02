// Copyright 2023, Christopher Dilks
// Subject to the terms in the LICENSE file found in the top-level directory.

#include <functional>
#include <iostream>
#include <unistd.h>
#include <spdlog/spdlog.h>

#include <TFile.h>

#include <podio/podioVersion.h>
#if podio_VERSION >= PODIO_VERSION(0, 99, 0)
#include <podio/ROOTReader.h>
#else
#include <podio/ROOTFrameReader.h>
#endif
#include <podio/Frame.h>

#include "SimHitAnalysis.h"
#include "RawHitAnalysis.h"
#include "CherenkovPIDAnalysis.h"
#include "ReconstructedParticleAnalysis.h"

using namespace benchmarks;

// -------------------------------------------------------------

// main logger
std::shared_ptr<spdlog::logger> m_log;

// get a collection, and check its validity
template<class C>
const C& GetCollection(podio::Frame& frame, std::string name) {
  const auto& coll = frame.get<C>(name);
  if(!coll.isValid())
    m_log->error("invalid collection '{}'", name);
  return coll;
}

// -------------------------------------------------------------
// MAIN
// -------------------------------------------------------------
int main(int argc, char** argv) {

  // -------------------------------------------------------------
  // setup and CLI option parsing
  // -------------------------------------------------------------

  // loggers
  m_log = spdlog::default_logger()->clone("benchmark_rich");
  auto log_level = spdlog::level::info;
  m_log->set_level(log_level);

  // default options
  int verbosity             = 0;
  long long num_events      = 0;
  std::string ana_file_name = "out_rich.root";
  std::vector<std::string> rec_files;
  std::vector<std::string> analysis_algorithms = {
    "SimHit",
    "RawHit",
    "CherenkovPID",
    "ReconstructedParticle"
  };
  auto list_algorithms = [&analysis_algorithms] () {
    for(auto analysis_algorithm : analysis_algorithms)
      std::cout << "        " << analysis_algorithm << std::endl;
    std::cout << std::endl;
  };

  // usage guide
  auto PrintUsage = [&argv, &list_algorithms, ana_file_name] () {
    std::cout << "\nUSAGE: " << argv[0] << " -i [INPUT FILES]... [OPTIONS]..." << std::endl
      << "\nREQUIRED OPTIONS:\n"
      << " -i [INPUT FILES]...\n"
      << "    input file from reconstruction output;\n"
      << "    can specify more than one (delimited by space)\n"
      << "\nOPTIONAL OPTIONS:\n"
      << " -o [OUTPUT FILE]\n"
      << "    output analysis file (default: " << ana_file_name << ")\n"
      << " -n [NUMBER OF EVENTS]\n"
      << "    number of events to process (default: process all)\n"
      << " -v: verbose output\n"
      << " -vv: more verbose output\n"
      << " -a [ALGORITHMS]...\n"
      << "    list of analysis algorithms to run (delimited by spaces)\n"
      << "    default: run all of them\n"
      << "    available algorithms:\n";
    list_algorithms();
  };
  if(argc==1) {
    PrintUsage();
    return 2;
  }

  // parse options
  int opt;
  while( (opt=getopt(argc, argv, "i:o:n:v|a:")) != -1) {
    switch(opt) {
      case 'i':
        optind--;
        for(; optind<argc && *argv[optind]!='-'; optind++)
          rec_files.push_back(std::string(argv[optind]));
        break;
      case 'o':
        ana_file_name = std::string(optarg);
        break;
      case 'n':
        num_events = std::atoll(optarg);
        break;
      case 'v':
        verbosity++;
        break;
      case 'a':
        analysis_algorithms.clear();
        optind--;
        for(; optind<argc && *argv[optind]!='-'; optind++)
          analysis_algorithms.push_back(std::string(argv[optind]));
        break;
      default:
        m_log->error("Unknown option '{}'", opt);
        PrintUsage();
        return 1;
    }
  }

  // print options
  m_log->info("{:-^50}", " User options ");
  m_log->info("  input files:");
  for(auto rec_file : rec_files)
    m_log->info("    {}", rec_file);
  m_log->info("  output file: {}", ana_file_name);
  m_log->info("  algorithms:");
  for(auto analysis_algorithm : analysis_algorithms)
    m_log->info("    {}", analysis_algorithm);
  m_log->info("  num_events: {}", num_events);
  m_log->info("  verbosity: {}", verbosity);
  m_log->info("{:-^50}", "");

  // check options
  std::vector<std::string> bad_options;
  if(rec_files.size()==0) bad_options.push_back("no [INPUT FILES] have been specified");
  if(analysis_algorithms.size()==0) bad_options.push_back("no [ALGORITHMS] have been specified");
  for(auto bad_option : bad_options) m_log->error(bad_option);
  if(bad_options.size()>0) {
    PrintUsage();
    return 1;
  }

  // set log level
  switch(verbosity) {
    case 0:  log_level = spdlog::level::info;  break;
    case 1:  log_level = spdlog::level::debug; break;
    default: log_level = spdlog::level::trace;
  }
  m_log->set_level(log_level);

  // start the output file
  auto ana_file = new TFile(TString(ana_file_name), "RECREATE");

  // struct to hold an algorithm execution processor
  struct AnalysisAlgorithm {
    std::shared_ptr<spdlog::logger>    log;     // algorithm-specific logger
    std::function<void(podio::Frame&)> process; // call AlgorithmProcess
    std::function<void()>              finish;  // call AlgorithmFinish
  };
  std::vector<AnalysisAlgorithm> algorithm_processors;

  // loop over algorithms to run, and define how to run them
  for(auto analysis_algorithm : analysis_algorithms) {
    AnalysisAlgorithm algo;

    // algorithm logger
    algo.log = spdlog::default_logger()->clone(analysis_algorithm);
    algo.log->set_level(log_level);
    algo.log->debug("{:-<50}","INITIALIZE ");

    // --------------------------------------------------------------
    // define how to run each algorithm
    // --------------------------------------------------------------

    // truth hits (photons) .........................................
    if(analysis_algorithm == "SimHit") {
      auto sim_algo = std::make_shared<SimHitAnalysis>();
      ana_file->mkdir("phot")->cd();
      sim_algo->AlgorithmInit(algo.log);
      algo.process = [sim_algo] (podio::Frame& frame) {
        sim_algo->AlgorithmProcess(
            GetCollection<edm4hep::MCParticleCollection>(frame,"MCParticles"),
            GetCollection<edm4hep::SimTrackerHitCollection>(frame,"DRICHHits")
            );
      };
      algo.finish = [sim_algo] () { sim_algo->AlgorithmFinish(); };
    }

    // digitizer ....................................................
    else if(analysis_algorithm == "RawHit") {
      auto digi_algo = std::make_shared<RawHitAnalysis>();
      ana_file->mkdir("digi")->cd();
      digi_algo->AlgorithmInit(algo.log);
      algo.process = [digi_algo] (podio::Frame& frame) {
        digi_algo->AlgorithmProcess(
            GetCollection<edm4eic::RawTrackerHitCollection>(frame,"DRICHRawHits"),
            GetCollection<edm4eic::MCRecoTrackerHitAssociationCollection>(frame,"DRICHRawHitsAssociations")
            );
      };
      algo.finish = [digi_algo] () { digi_algo->AlgorithmFinish(); };
    }

    // IRT ..........................................................
    else if(analysis_algorithm == "CherenkovPID") {
      std::vector<std::string> radiators = { "Aerogel", "Gas", "Merged" };
      std::map<std::string, std::shared_ptr<CherenkovPIDAnalysis>> irt_algos;
      for(auto radiator : radiators)
        irt_algos.insert({radiator, std::make_shared<CherenkovPIDAnalysis>()});

      for(auto& [radiator, irt_algo] : irt_algos) {
        ana_file->mkdir(("pid"+radiator).c_str())->cd();
        irt_algo->AlgorithmInit(radiator, algo.log);
      }

      algo.process = [irt_algos] (podio::Frame& frame) {
        const auto& mc_parts = GetCollection<edm4hep::MCParticleCollection>(frame,"MCParticles");
        for(auto& [radiator, irt_algo] : irt_algos) {
          const auto& cherenkov_pids = GetCollection<edm4eic::CherenkovParticleIDCollection>(
              frame,
              "DRICH" + radiator + "IrtCherenkovParticleID"
              );
          irt_algo->AlgorithmProcess(mc_parts, cherenkov_pids);
        }
      };

      algo.finish = [irt_algos] () {
        for(auto& [radiator, irt_algo] : irt_algos)
          irt_algo->AlgorithmFinish();
      };

    }

    // linking to reconstructed particles ...........................
    else if(analysis_algorithm == "ReconstructedParticle") {
      auto link_algo = std::make_shared<ReconstructedParticleAnalysis>();
      ana_file->mkdir("link")->cd();
      link_algo->AlgorithmInit(algo.log);
      algo.process = [link_algo] (podio::Frame& frame) {
        link_algo->AlgorithmProcess(
            GetCollection<edm4eic::MCRecoParticleAssociationCollection>(
              frame,
              "ReconstructedChargedParticleAssociations"
              )
            );
      };
      algo.finish = [link_algo] () { link_algo->AlgorithmFinish(); };
    }

    // unknown algorithm ............................................
    else {
      m_log->error("Unknown [ALGORITHM] '{}'; ignoring", analysis_algorithm);
      continue;
    }

    // --------------------------------------------------------------

    algorithm_processors.push_back(algo);
  }

  // -------------------------------------------------------------
  // read the input file and run all algorithms
  // -------------------------------------------------------------

  // open the input files
#if podio_VERSION >= PODIO_VERSION(0, 99, 0)
  podio::ROOTReader podioReader;
#else
  podio::ROOTFrameReader podioReader;
#endif
  m_log->warn("podio::ROOTFrameReader cannot yet support multiple files; reading only the first");
  // podioReader.openFiles(rec_files);
  podioReader.openFile(rec_files.front());

  // get the number of events to process
  const std::string tree_name = "events";
  long long num_entries = podioReader.getEntries(tree_name);
  if(num_events>0) num_entries = std::min(num_events, num_entries);

  // event loop
  m_log->info("BEGIN EVENT LOOP");
  for(long long e=0; e<num_entries; e++) {
    if(e%100==0) m_log->info("Progress: {:.3f}%", 100.0 * e/num_entries);
    auto podioFrame = podio::Frame(podioReader.readNextEntry(tree_name));
    for(auto& algo : algorithm_processors) {
      algo.log->debug("{:-<50}","PROCESS EVENT ");
      algo.process(podioFrame);
    }
  }

  // finish
  for(auto& algo : algorithm_processors) {
    algo.log->debug("{:-<50}","FINISH ");
    algo.finish();
  }

  // write output
  ana_file->Write();
  ana_file->Close();
  m_log->info("Done. Wrote analysis file:");
  m_log->info("  {}", ana_file_name);
}
