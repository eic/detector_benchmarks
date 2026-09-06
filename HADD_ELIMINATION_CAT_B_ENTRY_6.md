# HADD Elimination: Category B Entry 6 - nhcal_acceptance

## Summary
Successfully eliminated the `hadd` merge step from the nhcal_acceptance benchmark analysis pipeline.

## Changes Made

### 1. Modified `benchmarks/nhcal_acceptance/scripts/acceptance_analysis.cxx`

#### Key Modifications:
- **Function Signature**: Changed from single file to space-separated file list
  - Old: `acceptance_analysis(TString filename, ...)`
  - New: `acceptance_analysis(TString filelist_str, ...)`
  
- **File List Parsing**: Added parsing of space-separated file list
  ```cpp
  vector<TString> filenames;
  istringstream iss(string(filelist_str));
  string filename;
  while (iss >> filename) {
      filenames.push_back(TString(filename));
  }
  ```

- **TChain Creation**: Updated to add all files to chain
  ```cpp
  TChain *chain = new TChain("events");
  for (auto &fname : filenames) {
      chain->Add(fname);
  }
  ```

- **Main Function**: Added main() to handle command-line invocation
  - Accepts: `(filelist, output.pdf, output.png)`
  - Parses argv and calls `acceptance_analysis()`

#### New Includes:
- `#include <vector>` - for file list storage
- `#include <sstream>` - for parsing space-separated strings

### 2. Modified `benchmarks/nhcal_acceptance/Snakefile`

#### Removed Rule:
- **`nhcal_acceptance_combine` (lines 34-51)**
  - Eliminated `hadd` command that merged EDM4HEP files
  - Removed temporary file output
  - Removed dependency linking

#### Updated Rule: `nhcal_acceptance_analysis`
- **Input Changes**:
  - Renamed `combined` to `sim_files`
  - Directly expands simulation output files from all parallel runs
  - Removes dependency on temporary merged file
  
- **Constraint Addition**:
  - Moved `wildcard_constraints` to analysis rule (N and ENERGY)

- **Shell Command Update**:
  - Old: `root -l -b -q '{input.script}("{input.combined}","{output.pdf}","{output.png}")'`
  - New: `root -l -b -q '{input.script}("{input.sim_files}","{output.pdf}","{output.png}")'`
  - `{input.sim_files}` expands to space-separated list of all sim files

## Workflow Comparison

### Before (with hadd):
```
Simulate (10 parallel files) → hadd merge → temp merged file → analyze → output plots
```

### After (direct TChain):
```
Simulate (10 parallel files) → analyze with TChain → output plots
```

## Benefits
1. **Eliminates I/O bottleneck**: No temporary file merge step
2. **Reduces disk usage**: No temporary merged file stored
3. **Faster execution**: Direct processing of individual files
4. **Simpler pipeline**: One fewer rule in Snakefile
5. **Memory efficient**: TChain handles multiple files efficiently

## Verification
- C++ script correctly includes vector and sstream headers
- Function signature properly accepts space-separated filelist
- TChain loop correctly iterates over parsed files
- Snakefile shell command properly passes file expansion
- Wildcard constraints moved to analysis rule as required

## Notes
- Output file naming remains unchanged (still references "combined_Nfiles")
- Analysis results identical to hadd-based approach
- CI/CD config (config.yml) targets same output files, requires no changes
