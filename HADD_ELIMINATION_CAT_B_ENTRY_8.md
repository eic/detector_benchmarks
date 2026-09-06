# Hadd Elimination - Category B Entry 8: LFHCAL Performance

## Task Overview
Modify LFHCAL_Performance.C and the lfhcal_at_momentum Snakefile rule to eliminate temporary hadd file merges. The script now accepts multiple input files and creates a TChain internally.

## Changes Made

### 1. Modified: `benchmarks/lfhcal/LFHCAL_Performance.C`

#### Function Signature Change
**Before:**
```cpp
void LFHCAL_Performance(TString filename="tracking_output", TString particle="pi-", double mom=0.1, ...)
```

**After:**
```cpp
void LFHCAL_Performance(const char* file_list, const char* particle="pi-", double mom=0.1, ...)
```

#### Key Changes:
1. **Added Headers**
   - `#include "TChain.h"`
   - `#include <vector>`
   - `#include <string>`
   - `#include <sstream>`

2. **File Input Handling**
   - Changed from opening a single file with `TFile::Open()` to creating a `TChain`
   - Accepts space-separated filenames as a single string parameter
   - Parses filenames using `std::stringstream`
   - Adds each file to the TChain

3. **Implementation Details**
   ```cpp
   // Parse space-separated filenames from input string
   std::string file_str(file_list);
   std::stringstream ss(file_str);
   std::string filename;
   int file_count = 0;
   
   while (ss >> filename) {
     if (debug) cout << "Adding file to chain: " << filename << endl;
     chain->Add(filename.c_str());
     file_count++;
   }
   ```

4. **Parameter Type Updates**
   - Changed `output_dir` from `TString` to `const char*`
   - Changed `particle` from `TString` to `const char*`
   - Changed `name` from `TString` to `const char*`

### 2. Modified: `benchmarks/lfhcal/Snakefile` - `lfhcal_at_momentum` rule

#### Removed
- **Temporary output**: `combined_root=temp("{CAMPAIGN}/lfhcal_sim_{MOMENTUM}_{PARTICLE}.root")`
- **hadd command**: `hadd {output.combined_root} {input.sim}`

#### Updated
- **Output section**: Now contains only the primary output (performance ROOT file)
- **Shell command**: Passes space-separated input files directly to the ROOT macro:
  ```bash
  root -l -b -q 'LFHCAL_Performance.C("{input.sim}", "{wildcards.PARTICLE}", {wildcards.MOMENTUM}, 0.15, "", "{wildcards.CAMPAIGN}")'
  ```

## Behavior

### Old Workflow
1. Snakemake expands `{input.sim}` into space-separated filenames
2. hadd command merges files into temporary file: `lfhcal_sim_{MOMENTUM}_{PARTICLE}.root`
3. LFHCAL_Performance.C opens single merged file
4. Script processes and produces output: `lfhcal_mom_{MOMENTUM}_mom_resol_{PARTICLE}.root`
5. Temporary file cleaned up

### New Workflow
1. Snakemake expands `{input.sim}` into space-separated filenames
2. Space-separated filenames passed directly to ROOT macro as single string: `"file1.root file2.root file3.root"`
3. LFHCAL_Performance.C parses the string and creates TChain from individual files
4. TChain transparently merges tree entries from all files
5. Script processes merged tree and produces output: `lfhcal_mom_{MOMENTUM}_mom_resol_{PARTICLE}.root`
6. No temporary file needed

## Benefits

1. **Eliminates hadd overhead**: No temporary file creation/merging/cleanup
2. **Simpler pipeline**: Fewer intermediate steps
3. **Performance**: Direct chain processing typically faster than hadd merge
4. **Scalability**: Can handle arbitrary number of input files without modification
5. **Backward compatible**: Output format and naming unchanged

## Testing

The modified macro compiles successfully with ROOT's C++ parser. The function signature correctly accepts space-separated filenames as the first parameter and processes them through a TChain.

## Notes

- The script maintains 100% backward compatibility with existing output files
- The TChain approach is transparent to the analysis logic
- Debug output will show each file being added to the chain
- All eta bin histograms and momentum resolution calculations remain unchanged
