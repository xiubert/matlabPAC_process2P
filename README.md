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
| 2 | Condition split | Group tifs by treatment condition, then **auto-split each group by scan geometry** (`[H W zoom mFast mSlow]`, read from headers) so every motion-correction group is dimensionally *and* optically uniform; the main-zoom full frame keeps the base name, others get a `_<H>x<W>` suffix (e.g. `postZX1` + `postZX1_128x256`). **Alternate-zoom tifs** (a different FOV that can't share a group or reuse ROIs): ≤2 are dropped with a warning; >2 prompt to keep as a separate group (default **No** = drop; headless/`-batch` takes the default) |
| 3 | Motion correction | Optionally split multi-channel tifs to single-channel (channel 2); concatenate per group; run NoRMCorre non-rigid correction (params from data dims, so 256×128 corrects normally); write corrected tifs to `NoRMCorred/` for FISSA |
| 3b | Resume block | Uncomment `%{...%}` to reload motion-corrected data without re-running NoRMCorre (e.g. when drawing ROIs for a new cell type) |
| 4–5 | ROI drawing | Interactively draw ROIs on the motion-corrected stack via `TIFcatROIgui`; repeat for each **256×256** treatment condition; save per condition to `_moCorrROI_<condition>.mat` |
| 5b | ROI reuse (256×128) | **Auto, no manual lever:** detects any 256×128 condition by header and remaps the contained ROIs from the matching 256×256 condition's `_moCorrROI_` file via `remapROIfile` (centered crop, IDs preserved). No-op when no 256×128 condition exists. Errors if the geometry is not the expected zoom-matched centered crop |
| 6 | ROI matching | `intersectROIfiles` keeps only ROIs present in all conditions so the same cells are compared pre/post. **256×128 (spont) conditions are excluded** so they don't reduce the 256×256 stim ROI set |
| 7 | Raw F extraction | Extract rawF and motion-corrected rawF per ROI per tif (per condition, using that condition's ROI set); save to `_moCorr_Tifs_Params.mat` |
| 8 | FISSA (Python) | Run `FISSAviaMatlab_prePostTreatment.py [dataPath]` in a FISSA-enabled Python env. **FISSA needs a uniform ROI count per run**, so the driver auto-**groups the ROI files by count** and runs FISSA once per group. One group (no 256×128) → `FISSAoutput/matlab.mat` (legacy, unchanged). Mixed → `FISSAoutput/g<k>/matlab.mat` + a `groups.json` manifest |
| 9 | FISSA parsing | Load FISSA output; split map vs. stimulus trials; apply neuropil scaling (`corrected = ROI − scale × neuropil`; default 0.8); save `_tifFileList.mat`. If `groups.json` exists, `mergeFISSAgroups` attaches each tif's traces from its group **by name** (per-tif row count = its own ROI count); otherwise the legacy single-output path runs unchanged |
| 10 | Stim alignment | `stimParam2ROI` attaches stimulus parameters from `_Pulses.mat` files to the corrected traces; **resolves the matching ROI set per stim group by trace row-count**, so a 256×128 spont group uses its reduced (remapped) ROI set. Writes one `_raw` table per stimulus family |
| 11 | Per-stim analysis | `processAnimalStimFamilies` runs `processBPN2P` / `processCGC` for whichever families this animal has, adding dF/F and peak responses and writing the **processed** table each group step reads. Set `runPerStimProcessing = false` to do it by hand with edited parameters |

**Output:** `tifFileList` struct where `tifFileList.stim(n).SCALEDfissaFroi` is an `nROI × nFrames` array of motion- and neuropil-corrected fluorescence for the nth stimulus tif. Equivalent `.map` field for FRA/BF mapping tifs.

---

### 1b. 256×128 (10 Hz spontaneous) reuse path

Spontaneous sessions can be acquired at **256×128, 10 Hz** (double the 5 Hz frame rate) by halving `linesPerFrame` to 128 with `scanAngleMultiplierSlow = 0.5` at the **same zoom** — a centered vertical crop of the 256×256 field (drops the top/bottom 64 rows, 1:1 pixels). ROIs drawn on the 256×256 data are reused on the 256×128 spont tifs; only cells fully contained in the central crop survive. **5 Hz and 10 Hz traces are never pooled** — the 256×128 data follows the **Spont** analysis path as its own family.

The path runs through `processAnimal2P.m` with **no manual intervention**, assuming the 256×256 and 256×128 tifs share a recording folder and a treatment token (so the crop condition pairs with its 256² source):

1. **§2** auto-splits the mixed treatment group into a 256² condition and a `_128x256` condition.
2. **§3** motion-corrects each independently (NoRMCorre takes dims from the data).
3. **§4–5** draw ROIs on the 256² condition only.
4. **§5b** auto-remaps those ROIs into the 256×128 condition (`remapROItoAcq` → centered crop; errors on any zoom/shift/rotation mismatch, e.g. an accidental zoom≠1).
5. **§6** matches ROIs across the 256² conditions only; the crop condition is excluded.
6. **§7** extracts F per condition using each condition's own ROI set (256² stim tifs → full set; 128 spont tifs → reduced set).
7. **§8 FISSA — automatic per-ROI-count grouping.** FISSA's output is a fixed cells×trials grid, so a single run cannot mix an 18-ROI and an 11-ROI set. The driver groups the ROI files by count and runs FISSA once per group (explicit per-image ROI lists), writing `g<k>/matlab.mat` + a `groups.json` manifest (ordered tif names per group).
8. **§9** detects `groups.json` and calls `mergeFISSAgroups`, attaching each tif's neuropil traces from its group by name and computing `SCALEDfissaFroi` (per-tif row count = its own ROI count). No manifest → the legacy single-output path runs unchanged.
9. **§10** `stimParam2ROI` builds the per-family tables. Its spont branch resolves the **256×128 ROI set** (by trace row-count, via `resolveROIset`) and produces `<animal>_anmlROI_SpontstimTable.mat`, whose `anmlROIbyStim` rows hold **10 Hz traces from the fixed (remapped) ROIs**, ready for further spontaneous-activity analysis.

**ROI-reuse helpers** (`helperFcns/ROI/`): `remapROItoAcq.m` (geometry-validated centered crop of a `moCorROI` struct; regenerates `mask` for raw-F and `ROIcurveOrderedXY` for FISSA, preserving IDs) and `remapROIfile.m` (driver: load source ROIs + read src/tgt tif headers + remap + save the pipeline bundle). FISSA grouping: `mergeFISSAgroups.m` + the grouped `FISSAviaMatlab_prePostTreatment.py`. Tests: `tests/testRemapROItoAcq_centeredCrop.m`, `tests/testFISSAgrouping.m`, `tests/testAA0072_pipeline.m`.

---

### 2. Per-stimulus analyses — `stimulusSpecific/`

**`processAnimal2P` section 11 runs these for you** — this section is the reference for what they do, and for running one by hand with non-default parameters. Section 10 (`stimParam2ROI`) writes one stim-aligned `_raw` table per stimulus family (`<animal>_anmlROI_BPNstimTable_raw.mat`, `<animal>_anmlROI_CGCstimTable_raw.mat`). Each script below loads its family's `_raw` table, adds dF/F and peak metrics, plots, and writes the processed bundle to the same name **without** the `_raw` suffix. Each resolves `dataPath` via `uigetdir` if not already in the workspace.

**`_raw` two-stage convention (BPN and CGC).** Keeping the raw and processed artifacts distinct means re-running a `process*` script never mutates its own input, and an animal that has not been processed yet has no processed file for `aggregateStimGroup` to pick up silently. Downstream consumers read the processed file.

#### `processBPN2P.m` — band-pass noise (BPN), single animal

Two-stage: reads the `_raw` table from `stimParam2ROI`, writes the processed `<animal>_anmlROI_BPNstimTable.mat` (re-running overwrites the processed file but never the `_raw` input).

- Configurable pre-onset `baselineSec` (default 1 s); each row's dF/F is onset-normalized so stim onset lands at the same frame regardless of its recorded `BPNsOnset`
- `combineDiffOnset` merges same-stim / different-onset rows (onset-aligning the raw F cells too)
- Trial-averages dF/F per row (`dFF_avg`), then runs `pkFcalc` on the cell-average to get peak dF/F + significance
- Plots: single ROI × dB, single ROI all dB, population (between-cell SEM per dB), peak dF/F vs. dB

#### `processRLF.m` — rate-level function for a group

Interactive entry point for the BPN group plots: select a `BPN_Group<g>.mat` and it calls `plotBPNgroup`.

- Builds per-cell RLFs and dB thresholds via `tableRLF`; cells are included only with ≥ `nConsec` consecutive significant dB levels
- Plots the group-mean RLF, population dF/F re sound level, and peak dF/F re level
- Assembling animals into a group is `aggregateStimGroup`'s job (see step 2), so this script no longer carries a list of per-animal paths

#### `processCGC.m` — pure-tone-in-contrast (contrast gain control), single animal

Two-stage: reads `<animal>_anmlROI_CGCstimTable_raw.mat`, writes `<animal>_anmlROI_CGCstimTable.mat`. A pre-split animal with only the un-suffixed file is still handled — it is re-processed with its derived columns stripped first.

- dF/F referenced first to a pre-DRC baseline (`dFF_DRC`), then additively to a pre-pure-tone baseline (`dFF_PT`), matching the manuscript method
- Peak PT responses + significance from the **cell-average** trace (`dFF_PT_avg`), via `pkFcalc`. The `t >= 1` crop before `pkFcalc` is load-bearing: it makes the significance baseline equal the `F0_PT` window
- Plots: per-ROI dF/F grids (3×3), then delegates the population panels to `plotCGCgroup` so the single-animal and group cases run identical code

#### `processFRA.m` — frequency response area mapping

Operates on `tifFileList` directly (not via a per-stim table); computes per-cell FRA maps and best-frequency estimates.

---

### 3. Condition groups — `aggregateStimGroup` + the group plotters

The path for comparing **treatment/condition groups** (Group A vs B vs …).

#### End-to-end, from raw tifs to a group figure

```
PER ANIMAL, once each                     ONCE PER GROUP
──────────────────────────────────        ─────────────────────────────
processAnimal2P                           aggregateStimGroup(manifest)
  §1–9   motion correction, FISSA            → <Family>_Group<g>.mat
  §10    stimParam2ROI  → _raw tables               │
  §11    processBPN2P / processCGC                  ▼
         → processed tables (dF/F, peaks)     plotBPNgroup / plotCGCgroup
```

```matlab
% 1. per animal (repeat for every animal in the group)
processAnimal2P                     % §10 writes _raw, §11 writes processed

% 2. once per group, per family
aggregateStimGroup(manifest)        % -> BPN_GroupD.mat + BPN_GroupD_manifest.json

% 3. plot
plotBPNgroup('BPN_GroupD.mat')
plotCGCgroup('CGC_GroupD.mat')
```

Aggregation reads the **processed** table, not the `_raw` one, because dF/F and peak responses must be computed per animal — baselines and significance are per-animal quantities. `processAnimal2P` §11 does that automatically for every family the animal has. If it is skipped (`runPerStimProcessing = false`, or an older animal processed before §11 existed), `aggregateStimGroup` raises `aggregateStimGroup:notProcessed` naming the animal and the script to run, rather than building a group missing its analysis columns — the point of the [`_raw` two-stage convention](#2-per-stimulus-analyses--stimulusspecific).

**Build a group** from a manifest — a `.json` kept next to the data, so membership is reviewable and diffable:

```matlab
manifest = struct('group','D','family','BPN', ...
                  'animals',["TO0006","TO0007"], ...
                  'cohortRoot','/media/DATA/Ophys/Jinbo/TOMT', ...
                  'outDir','/media/DATA/Ophys/Jinbo/TOMT/aggregate data');
aggregateStimGroup(manifest);            % or aggregateStimGroup('BPN_GroupD_manifest.json')
```

Writes `<Family>_Group<g>.mat` containing `anmlROIbyStim` (unchanged variable name — existing scripts keep working) plus a `groupInfo` stamp recording the animals, per-animal cell counts, time axis, source files, and **the analysis convention in force**. That last field is what makes a convention mismatch between groups visible instead of silent.

**Plot a group** — both handle any group size, including a single animal, and state cells/mice on every panel:

```matlab
plotBPNgroup('BPN_GroupD.mat');    % RLF, dF/F re sound level, peak re level
plotCGCgroup('CGC_GroupD.mat');    % contrast traces, low-vs-high scatter, paired bar
```

`processRLF.m` is the interactive entry point for the BPN plots — it prompts for a group file and calls `plotBPNgroup`.

**Chase an outlier** seen in a group plot back to the cell and its individual repetitions:

```matlab
plotCGCgroup('CGC_GroupD.mat','showCells',true)   % every cell, faint, behind the mean
plotCellTrials('CGC_GroupD.mat',"TO0007","7")     % every repetition of one cell
```

`makeGroupFigures` (in `etc/`) renders the whole standard set — per-group summaries, the cross-group RLF, significant-cell counts, and the per-cell trace check — into `etc/figures/`.

**Check a group** without plotting:

```matlab
validateStimGroup('CGC_GroupD.mat')   % required columns, cell identity, time axis,
                                      % trace widths, frame rates, column order
loadStimGroup('CGC_GroupD.mat')       % load + validate + cross-check the manifest
```

Small-*n* handling is a contract enforced in `helperFcns/cohort/`, not per-plot: zero cells renders a labelled empty panel, one cell gets no SEM band, a single-animal group is stamped as carrying no across-animal inference, and a test below `minN` reports why it was skipped rather than emitting `p = NaN`.

> `combinetable.m` (in the data directory) is **retired** — it recorded no membership, validated nothing, and its "multiple regions" branch corrupted `roiID` by doing arithmetic on a string. Use `aggregateStimGroup`.

---

### 4. Compile cohort data — `compileCohortData.m`

A separate, older path that builds the flat `Tinput` table for the manuscript figures in `plotCohortData.m` / `matlabPAC_CGCplot/plotDataTable.m`. Unlike step 2 it re-derives dF/F at plot time and re-indexes IDs as sequential integers. Edit the parameters block at the top (`cohortName`, `family`).

- `compileAnmlFRA` — collects FRA maps and BF estimates across animals
- `compileAnmlROItables(family, params)` — collects stimulus/response tables per animal, validating each against its family contract before reducing columns
- Joins BF data into the response table; keeps only tone-responsive cells (`dPrime > 0`)
- Saves `<cohortName>_dataTable.mat` and `<cohortName>_params.mat`

---

### 5. Plot cohort data — `plotCohortData.m`

Loads the compiled cohort table and produces analysis figures. Each analysis block is wrapped in `%{...%}` and run independently.

- Computes dF/F relative to DRC and pure-tone (PT) baseline windows
- Extracts peak dF/F responses at PT onset using `pkFcalc`
- Analyses include: PT response ratio vs. PT onset, sustained DRC response, pre/post treatment comparisons
- Statistical testing via `sigDiffCalc` (parametric or Wilcoxon) with permutation test fallback

---


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

stimulusSpecific/
  processBPN2P.m      — per-animal band-pass noise (BPN) dF/F + peak analysis
  processCGC.m        — per-animal pure-tone-in-contrast (contrast gain control)
  processFRA.m        — frequency response area mapping from tifFileList
  plotBPNgroup.m      — BPN group plots: RLF, dF/F re level, peak re level
  plotCGCgroup.m      — CGC group plots: contrast traces, low-vs-high, paired bar
  plotCellTrials.m    — every repetition of one cell, for chasing outliers
  processRLF.m        — interactive entry point for the BPN group plots

helperFcns/
  tif/                — ScanImage tif reading (readSCIMtif, justLoadTif),
                        channel splitting (splitTifChans), writing (writeMoCorTifs)
  dFF/                — dF/F computation (dFoFcalc), peak response detection (pkFcalc)
  FRA/                — FRA map construction (FRAmap), BF extraction, d-prime calculation
  RLF/                — rate-level functions (cellRLF, tableRLF) and plotting (plotRLF)
  cohort/             — group-size-agnostic primitives shared by both family
                        plotters: counts (groupN), trace/value stacking
                        (gatherCellTraces, gatherCellValues), between-cell
                        mean/SEM (cohortMeanSEM), guarded statistics
                        (cohortStat), figure n-stamps (annotateN)
  dataOrg/            — FISSA output parsing, tifFileList assembly, stimulus table
                        builders (stimParam2ROI, combineDiffOnset); condition-group
                        aggregation (stimGroupSpec, aggregateStimGroup,
                        validateStimGroup, loadStimGroup); per-animal
                        per-family driver (processAnimalStimFamilies)
  ROI/                — ROI mask ↔ polygon conversion, raw F extraction from masks
  plotting/           — SEM shaded plots (fillSEMplot), regression plots (regPlot),
                        paper-style formatting, figure export
  sound/              — Ephus .signal file inspection (inspectSignalObject)
  general/            — utility functions (SEMcalc, sigDiffCalc, scaleZeroToOne, zero2nan)
```
