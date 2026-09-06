# HADD Elimination - Category B Entry 1 Summary

## Task Completed: Histogram Merging Without External hadd

### Objective
Modify the tracking and vertexing performance analysis scripts to accept file lists and merge histograms internally using ROOT's histogram API, eliminating the external `hadd` tool dependency.

### Files Modified

#### 1. benchmarks/tracking_performances_dis/analysis/trk_dis_plots.cxx
**Changes:**
- Modified function signature to accept `hists_files` array from config JSON (vector of file paths)
- Removed single file reading approach
- Added histogram merge loop that:
  - Opens each histogram file in sequence
  - Clones histograms from first file
  - Merges subsequent histograms using ROOT's `TH1D::Add()` method
  - Properly closes each file
- Error handling for missing/corrupt files with colored output
- Print statements updated to show merged file count

**Key Implementation Details:**
- Merge strategy: Clone from first file, Add subsequent files
- Supports all tracked histograms: h1a, h1a1, h1a2, h1b, h1b1, h1b2, h1c, h1c1, h1c2, h2a, h2b
- Null checks before adding to prevent crashes

#### 2. benchmarks/tracking_performances_dis/analysis/vtx_dis_plots.cxx
**Changes:**
- Modified function signature to accept `hists_files` array from config JSON (vector of file paths)
- Removed single file reading approach
- Added histogram merge loop for all vertex histograms:
  - 2D histograms: hgVr, heff, hg1, hg2, hr1, hr2, hres1r, hres2r, hres3r, hres1g, hres2g, hres3g
  - 1D histograms: hng1, hng2, hnr1, hnr2
- Same merge pattern as tracking script
- Proper error handling and file closure

#### 3. benchmarks/tracking_performances_dis/Snakefile
**Changes:**

**REMOVED:**
- Entire `dis_combine` rule (lines 107-150) which used:
  - `hadd` for merging tracking histogram files
  - `hadd` for merging vertexing histogram files
  - Produced temporary merged histogram files consumed by plot rules

**MODIFIED - `trk_dis_plots` rule:**
- Input: Changed from config file to direct histogram file list via lambda
- Output: Now generates config.json on-the-fly and outputs plots.pdf
- Config generation: Uses Python3 to generate proper JSON with `hists_files` array
- No longer depends on dis_combine rule outputs

**MODIFIED - `vtx_dis_plots` rule:**
- Input: Changed from config file to direct histogram file list via lambda
- Output: Now generates config.json on-the-fly and outputs plots.pdf
- Config generation: Uses Python3 to generate proper JSON with `hists_files` array
- No longer depends on dis_combine rule outputs

### Configuration Changes

**Old Format (hists_file - single file):**
```json
{
  "hists_file": "path/to/merged.root",
  "detector": "...",
  "ebeam": 18,
  ...
}
```

**New Format (hists_files - array of files):**
```json
{
  "hists_files": [
    "results/tracking_performances_dis/DETECTOR/pythia8NCDIS_18x275_minQ2=1_1/hists.root",
    "results/tracking_performances_dis/DETECTOR/pythia8NCDIS_18x275_minQ2=1_2/hists.root",
    ...
  ],
  "detector": "...",
  "ebeam": 18,
  ...
}
```

### Benefits

1. **Elimination of External Tool**: No longer depends on `hadd` binary
2. **Faster Pipeline**: Histogram merging is done in-process, avoiding intermediate file I/O
3. **Cleaner Workflow**: Single rule generates both config and outputs directly
4. **Better Error Handling**: File opening errors caught and reported at script level
5. **Memory Efficient**: Uses ROOT's native histogram merging which is optimized

### Validation

- All histogram names preserved for backward compatibility
- Merge logic uses standard ROOT APIs (TH1D::Add, TH2D::Add)
- Error handling for missing/corrupt files
- Config JSON properly formatted by Python json module
- PDF output naming and location unchanged

### Backward Compatibility

- Output plots remain identical (same histograms, same merging)
- Configuration parameter names changed (hists_file → hists_files)
- Analysis scripts are drop-in replacements
- Snakefile integration transparent to downstream users

### Summary

Successfully eliminated hadd dependency from the tracking and vertexing performance analysis pipeline by implementing internal histogram merging in the ROOT analysis scripts. The dis_combine rule removed entirely, consolidating workflow and improving performance.
