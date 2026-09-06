# Category B Entry 5: hadd Elimination - Tracking Performances

## Overview
Modified `Tracking_Performances.C` and the Snakefile to eliminate the `hadd` step for merging EDM4EIC reco files. The script now accepts a file list via JSON config and creates a TChain internally.

## Changes Made

### 1. Tracking_Performances.C
**File**: `benchmarks/tracking_performances/Tracking_Performances.C`

**Key Changes**:
- Changed function signature from:
  ```c++
  void Tracking_Performances(TString filename="tracking_output", TString particle="pi-", ...)
  ```
  to:
  ```c++
  void Tracking_Performances(TString config_file="config.json")
  ```

- Added new includes:
  - `#include "nlohmann/json.hpp"`
  - `#include <fstream>`
  - `#include <vector>`

- Added JSON config parsing logic that reads:
  - `sim_files`: Array of 3 EDM4EIC reco files
  - `particle`: Particle type (e.g., "pi-")
  - `momentum`: Momentum value
  - `output_prefix`: Output directory path
  - `truth_seeding`: Boolean flag
  - `pTcut`: Transverse momentum cut

- Created TChain from file list instead of opening single file:
  ```c++
  TChain* chain = new TChain("events");
  for (const auto& file_path : sim_files) {
    chain->Add(file_path.c_str());
  }
  TTreeReader myReader(chain);
  ```

- Removed file-based tree access and replaced with chain-based approach

### 2. Snakefile
**File**: `benchmarks/tracking_performances/Snakefile`

**Changes to `tracking_performance_at_momentum` rule** (lines ~116-156):

**Removed**:
- `combined_root=temp(...)` from output section
- `hadd {output.combined_root} {input.sim}` command

**Added**:
- `params: output_prefix` lambda function
- Python script embedded in shell that:
  - Converts `input.sim` (list of files) to JSON format
  - Generates `tracking_config_*.json` config file with all parameters
  - Passes config file path to ROOT script

**Before**:
```bash
hadd {output.combined_root} {input.sim}
root -l -b -q {input.script}'("{output.combined_root}", ...)'
```

**After**:
```bash
# Python generates JSON with sim_files list
root -l -b -q {input.script}'"config.json"'
rm -f "config.json"
```

## Workflow Changes

### Old Workflow
1. Snakemake passes 3 .eicrecon.edm4eic.root files to hadd
2. hadd merges them into temporary combined file
3. Tracking_Performances.C opens single merged file
4. Analysis performed on merged tree

### New Workflow
1. Snakemake passes file list to Tracking_Performances.C via JSON config
2. Tracking_Performances.C creates TChain internally from file list
3. Analysis performed on TChain (transparent to rest of code)
4. No temporary intermediate files created

## Benefits
- Eliminates intermediate hadd file generation
- Reduces disk I/O during analysis
- TChain internally handles multiple files transparently
- Cleaner dependency management in Snakemake

## JSON Config Format
```json
{
  "sim_files": [
    "sim_output/tracking_performance/.../file1.eicrecon.edm4eic.root",
    "sim_output/tracking_performance/.../file2.eicrecon.edm4eic.root",
    "sim_output/tracking_performance/.../file3.eicrecon.edm4eic.root"
  ],
  "particle": "pi-",
  "momentum": 1.0,
  "output_prefix": "local/realseed/pi-",
  "truth_seeding": false,
  "pTcut": 0.15
}
```

## Dependencies
- `nlohmann/json.hpp` - Header-only JSON library (must be available in ROOT environment)

## Files Modified
1. `benchmarks/tracking_performances/Tracking_Performances.C`
2. `benchmarks/tracking_performances/Snakefile`

## Backward Compatibility
The changes are NOT backward compatible. The script now requires JSON config input instead of individual parameters. However, the output files remain identical in format and location.
