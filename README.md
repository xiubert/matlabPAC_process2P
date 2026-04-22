# matlabPAC_process2P

MATLAB pipeline for processing two-photon calcium imaging data, from raw ScanImage `.tif` files through motion correction, neuropil subtraction, and cohort-level analysis.

---

## Analysis workflow

### 0. Ad hoc FRA map during experiment — `adhocFRAmap.m`

Run during or immediately after a recording session to get a quick best-frequency (BF) / frequency response area (FRA) map without running the full pipeline.

- Select the animal data directory and the tif files used for FRA mapping
- ROIs are drawn with `roiGUI` on the first map tif and reused for the rest
- Calls `FRAmap` to compute peak fluorescence responses across tone frequencies and amplitudes
- Plots FRA map and mean response curve via `plotFRAmap`
- Output saved to `adHocMap/` in the data directory

### 0b. Ad hoc trace inspection — `adhoc2P.m`

Quick look at average cell responses across a set of tifs without running the full pipeline.

- Requires ROIs drawn by `roiGUI` on the first selected tif
- Applies those ROIs to all selected tifs, computes dF/F per cell, and plots individual and population-mean traces

---

### 1. Per-animal processing — `processAnimal2P.m`

Full processing pipeline for a single animal session. Edit `dataPath` at the top and run section by section.

| Step | Section | What it does |
|------|---------|--------------|
| 1 | Tif inventory | List all tifs, assign pre/post treatment labels, flag FRA map tifs, save `_tifFileLegend.mat` |
| 2 | Condition split | Group tifs by treatment condition so each group is motion-corrected independently |
| 3 | Motion correction | Optionally split multi-channel tifs to single-channel (channel 2); concatenate per group; run NoRMCorre non-rigid correction; write corrected tifs to `NoRMCorred/` for FISSA |
| 3b | Resume block | Uncomment `%{...%}` to reload motion-corrected data without re-running NoRMCorre (e.g. when drawing ROIs for a new cell type) |
| 4–5 | ROI drawing | Interactively draw ROIs on the motion-corrected stack via `TIFcatROIgui`; repeat for each treatment condition; save per condition to `_moCorrROI_<condition>.mat` |
| 6 | ROI matching | `intersectROIfiles` keeps only ROIs present in all conditions so the same cells are compared pre and post treatment |
| 7 | Raw F extraction | Extract rawF and motion-corrected rawF per ROI per tif; save to `_moCorr_Tifs_Params.mat` |
| 8 | FISSA (Python) | Run `FISSAviaMatlab_prePostTreatment.py` in a FISSA-enabled Python environment; output written to `NoRMCorred/FISSAoutput/matlab.mat` |
| 9 | FISSA parsing | Load FISSA output; split map vs. stimulus trials; apply neuropil scaling factor (`corrected = ROI − scale × neuropil`; default 0.8); save `_tifFileList.mat` |
| 10 | Stim alignment | `stimParam2ROI` attaches stimulus parameters from `_Pulses.mat` files to the corrected traces |

**Output:** `tifFileList` struct where `tifFileList.stim(n).SCALEDfissaFroi` is an `nROI × nFrames` array of motion- and neuropil-corrected fluorescence for the nth stimulus tif. Equivalent `.map` field for FRA/BF mapping tifs.

---

### 2. Compile cohort data — `compileCohortData.m`

Aggregates per-animal outputs into a single cohort-level table. Edit the parameters block at the top.

- `compileAnmlFRA` — collects FRA maps and BF estimates across animals
- `compileAnmlROItables` — collects stimulus/response tables per animal
- Joins BF data into the response table; keeps only tone-responsive cells (`dPrime > 0`)
- Re-indexes animal and ROI IDs as sequential integers
- Saves `<cohortName>_dataTable.mat` and `<cohortName>_params.mat`

---

### 3. Plot cohort data — `plotCohortData.m`

Loads the compiled cohort table and produces analysis figures. Each analysis block is wrapped in `%{...%}` and run independently.

- Computes dF/F relative to DRC and pure-tone (PT) baseline windows
- Extracts peak dF/F responses at PT onset using `pkFcalc`
- Analyses include: PT response ratio vs. PT onset, sustained DRC response, pre/post treatment comparisons
- Statistical testing via `sigDiffCalc` (parametric or Wilcoxon) with permutation test fallback

---

## External dependencies

| Dependency | Language | Purpose |
|------------|----------|---------|
| [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) | MATLAB | Non-rigid motion correction |
| [FISSA](https://github.com/rochefort-lab/fissa) | Python | Neuropil signal separation |
| Ephus `@signalObject` library (`ephus_library`) | MATLAB | Reading `.signal` files in `inspectSignalObject` |
| ScanImage | — | Acquisition; tif files contain SI headers parsed throughout |

---

## Helper functions

```
GUIs/
  roiGUI.m            — draw and save ROIs on a single tif; outputs _roiOutput.mat
  TIFcatROIgui.m      — draw ROIs on a motion-corrected concatenated stack
  meanFluoROIvt.m     — interactive ROI selection and mean fluorescence extraction

helperFcns/
  tif/                — ScanImage tif reading (readSCIMtif, justLoadTif),
                        channel splitting (splitTifChans), writing (writeMoCorTifs)
  dFF/                — dF/F computation (dFoFcalc), peak response detection (pkFcalc)
  FRA/                — FRA map construction (FRAmap), BF extraction, d-prime calculation
  dataOrg/            — FISSA output parsing, tifFileList assembly, stimulus table builders
  ROI/                — ROI mask ↔ polygon conversion, raw F extraction from masks
  plotting/           — SEM shaded plots (fillSEMplot), regression plots (regPlot),
                        paper-style formatting, figure export
  sound/              — Ephus .signal file inspection (inspectSignalObject)
  general/            — utility functions (SEMcalc, sigDiffCalc, scaleZeroToOne, zero2nan)
```
