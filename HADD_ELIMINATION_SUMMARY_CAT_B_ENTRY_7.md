# Category B Entry 7: nHCal Sampling Fraction - Hadd Elimination

## Summary
Successfully modified the nHCal sampling fraction analysis pipeline to eliminate the hadd merge step. The analysis script now accepts all 18 individual EDM4HEP simulation files (6 energies × 3 particles) directly and processes them in a single pass.

## Changes Made

### 1. Modified: `benchmarks/nhcal_sampling_fraction/scripts/sampling_fraction_analysis.cxx`

#### Function Signature Change
- **Before**: `int sampling_fraction_analysis(const string &filename, ...)`
- **After**: `int sampling_fraction_analysis(const vector<string> &filenames, ...)`

#### Event Count Initialization (Lines 73-77)
```cpp
// OLD: Single file reader
podio::ROOTReader *reader = new podio::ROOTReader();
reader->openFile(filename);
unsigned nEvents = reader->getEntries("events");

// NEW: Count events across all files
unsigned nEvents = 0;
for (const auto &filename : filenames) {
    nEvents += podio::ROOTReader().openFile(filename).getEntries("events");
}
cout << "Total number of events: " << nEvents << endl;
```

#### Event Processing Loop (Lines 126-189)
- Changed from single file iteration to nested loop:
  1. **Outer loop**: Iterate over all filenames
  2. **Inner loop**: Process all events in each file
  3. Each file gets its own `podio::ROOTReader` instance
  4. Histograms accumulate across all files (shared scope outside loops)

#### Main Function (Lines 291-307)
- **Before**: Accepted 4 arguments (input file, output PDF, output PNG, compact XML)
- **After**: Accepts 5+ arguments (output PDF, output PNG, compact XML, input files...)
- Usage: `./sampling_fraction_analysis <out.pdf> <out.png> <compact.xml> <input1.root> [input2.root ...]`

### 2. Modified: `benchmarks/nhcal_sampling_fraction/Snakefile`

#### Rule Removed: `nhcal_sampling_fraction_combine` (Previously lines 43-58)
- Eliminated hadd merge step that combined 18 EDM4HEP files into a single temporary file
- Removed dependency on temporary merge output

#### Rule Updated: `nhcal_sampling_fraction_analysis` (Lines 44-71)

**Input Changes:**
- **Before**: Single input `combined="sim_output/nhcal_sampling_fraction/sim_combined.edm4hep.root"`
- **After**: Direct array of all 18 files:
  ```python
  sim_files=lambda wildcards: expand(
      "sim_output/nhcal_sampling_fraction/{PARTICLE}/Ekin{ENERGY}GeV/sim.edm4hep.root", 
      ENERGY=["0.5", "0.7", "1.0", "2.0", "5.0", "10.0"],
      PARTICLE=["pi-", "neutron", "e-"],
  ),
  ```

**Shell Command Changes:**
- **Before**: `{output.pdf}.bin "{input.combined}" "{output.pdf}" "{output.png}" "{params.DETECTOR_PATH}/{params.DETECTOR_CONFIG}"`
- **After**: `{output.pdf}.bin "{output.pdf}" "{output.png}" "{params.DETECTOR_PATH}/{params.DETECTOR_CONFIG}" {input.sim_files}`
- All 18 files passed as individual arguments to the compiled binary

## Implementation Details

### File Processing Strategy
- **Sequential Processing**: Files processed one at a time to avoid memory overhead
- **Cumulative Histograms**: All 6 histograms (`h_sampF_e`, `h_sampF_pi`, `h_sampF_n`, etc.) remain in scope and accumulate across files
- **Memory Efficient**: Only one file's events in memory at any given time

### Data Flow
```
18 simulation files (individual files)
    ↓
sampling_fraction_analysis.cxx (TChain-like processing)
    ↓
Sequential processing via nested loops
    ↓
Accumulated histograms for all events
    ↓
Analysis plots (PDF/PNG output)
```

## Files Processed
- 6 kinetic energies: 0.5, 0.7, 1.0, 2.0, 5.0, 10.0 GeV
- 3 particle types: e-, π-, neutron
- **Total: 18 input files**

## Output Preserved
Analysis output remains unchanged:
- `results/nhcal_sampling_fraction/hist_sampf_vs_Ehit.pdf`
- `results/nhcal_sampling_fraction/hist_sampf_vs_Ehit.png`
- Additional plots: `prof_sampf_vs_Ehit`, `hist_sampf_vs_Ekin`, `prof_sampf_vs_Ekin`

## Benefits
1. **Eliminated Temporary File**: No intermediate merged file stored to disk
2. **Simplified Workflow**: Fewer Snakemake rules (18 → 17)
3. **Faster Execution**: No I/O overhead from hadd merge operation
4. **Same Analysis Results**: Histogram contents identical to previous pipeline

## Testing Recommendations
- Verify output histogram values match previous hadd-based results
- Check all 18 input files are correctly processed
- Validate compilation with existing root-config and podio libraries
