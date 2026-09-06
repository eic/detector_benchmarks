# HADD Elimination Category B Entry 10

**Task:** Modify basic_distribution_analysis.cxx and remove nhcal_basic_distribution_combine rule

## Summary of Changes

### 1. Modified basic_distribution_analysis.cxx
- **Changed function signature** from accepting a single filename to accepting `vector<string>` of filenames
- **Implemented file loop** to iterate over all input files and calculate total event count
- **Refactored event processing** to loop over files and process each file independently
- **Added main() function** to parse command-line arguments:
  - Usage: `./binary <out.png> <out.root> <compact.xml> <input1.root> [input2.root ...]`
  - Accepts space-separated list of input files as command-line arguments
- **Output remains unchanged**: Same analysis plots (PNG) and ROOT files

### 2. Modified Snakefile - benchmarks/nhcal_basic_distribution/Snakefile
- **Removed rule**: `nhcal_basic_distribution_combine` (lines 38-55 in old version)
  - Previously used `hadd` to merge parallel EDM4HEP sim files into temp file
  - Temp file is no longer needed

- **Updated rule**: `nhcal_basic_distribution_analysis` (now lines 38-71)
  - **Input changes**:
    - Removed: dependency on merged temp file (`combined=...`)
    - Added: direct dependency on all sim output files (`sim_files=lambda wildcards: expand(...)`)
  - **Build & Execute changes**:
    - Now compiles the script to a binary using g++ with required dependencies
    - Executes binary with file list directly: `{params.binary} "{output.png}" "{output.root}" "{params.DETECTOR_PATH}/{params.DETECTOR_CONFIG}" {input.sim_files}`
    - Cleans up binary after execution

## Workflow Change
**Before:**
```
nhcal_basic_distribution_simulate (10 parallel files)
          ↓
nhcal_basic_distribution_combine (hadd → temp file)
          ↓
nhcal_basic_distribution_analysis (process temp file)
```

**After:**
```
nhcal_basic_distribution_simulate (10 parallel files)
          ↓
nhcal_basic_distribution_analysis (process all files directly via TChain)
```

## Implementation Pattern
Follows the same pattern as `benchmarks/nhcal_sampling_fraction/scripts/sampling_fraction_analysis.cxx`:
- Accepts multiple files as command-line arguments
- Creates file vector and iterates to count total events
- Opens each file independently and processes events
- Produces merged analysis output

## Benefits
- ✅ Eliminates temporary merged files (disk space savings)
- ✅ Reduces I/O overhead (no hadd step)
- ✅ Simpler dependency graph
- ✅ Consistent with other benchmarks in the project
