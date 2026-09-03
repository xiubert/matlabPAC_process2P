# matlabPAC_process2P

MATLAB pipeline for processing two-photon calcium imaging data, from raw ScanImage `.tif` files through motion correction, neuropil subtraction, and cohort-level analysis.

ROIs come either from **you, drawing them** (`processAnimal2P`) or from
**Cellpose, unattended** (`processAnimal2Pheadless`). Both feed the same
downstream pipeline.

📖 **[usage.md](usage.md) — what to type, in what order, for either path**, and
[where artifacts land](usage.md#where-artifacts-are-written): pass a run name and
each analysis gets its own `analysis/<run>/` folder instead of overwriting the last.
This README is the reference for what each piece *does*.

---

## Pipeline at a glance

```
 RAW ACQUISITION                      per animal, once each
 ───────────────                      ─────────────────────
 *.tif  *_Pulses.mat                  processAnimal2P.m          (interactive)
                                      processAnimal2Pheadless.m  (unattended; §1c)
        │                               §1  tif inventory        → _tifFileLegend.mat
        │                               §2  condition split      → _tifCondSplitLegend.mat
        ▼                               §3  NoRMCorre            → NoRMCorred/*.tif
                                        §4–5 ROI drawing  [GUI]  → _moCorrROI_<cond>.mat
                                        §6  ROI matching
                                        §7  raw F extraction     → _moCorr_Tifs_Params.mat
                                        §8  FISSA        [PYTHON]→ NoRMCorred/FISSAoutput/matlab.mat
                                        §9  FISSA parsing        → _tifFileList.mat
                                       ─────────────────────────────────────────────────
                                        §10 stimParam2ROI        → _anmlROI_<Fam>stimTable_raw.mat
                                        §11 processAnimalStim-   → _anmlROI_<Fam>stimTable.mat
                                            Families                (dF/F + peak responses)
                                                                 → _FRAmap.mat
                                                                 → _anmlROI_FRAtable.mat
                                                    │               (FRA has no _raw stage: it is
                                                    │                driven from _tifFileList.mat)
 ONCE PER GROUP, PER FAMILY                         ▼
 ──────────────────────────      aggregateStimGroup(manifest)
                                   validate each animal, canonicalise column order,
                                   concatenate, stamp provenance
                                        → <Fam>_Group<g>.mat   { anmlROIbyStim , groupInfo }
                                                       (FRA: anmlROIbyFRA)
                                        → <Fam>_Group<g>_manifest.json
                                                    │
                                                    ▼
                       plotBPNgroup / plotCGCgroup / plotFRAgroup / makeGroupFigures
```

### What each stage needs on disk

| To run | Needs in the animal folder | Produced by |
|---|---|---|
| §3 motion correction | `*.tif` | acquisition |
| §8 FISSA | `NoRMCorred/*.tif`, `FISSA_ROIs.npy` | §3, §4–5 |
| §9 FISSA parsing | `NoRMCorred/FISSAoutput/matlab.mat` | §8 |
| **§10 stim alignment** | **`<animal>_tifFileList.mat`**, **`<animal>_moCorrROI_*.mat`**, **`*_Pulses.mat`** | §9, §4–5, acquisition |
| §11 per-stim analysis | `<animal>_anmlROI_<Fam>stimTable_raw.mat` | §10 |
| §11 **FRA** | `<animal>_tifFileList.mat` (with a non-empty `.map`) + `*_Pulses.mat` | §9, acquisition |
| `aggregateStimGroup` | `<animal>_anmlROI_<Fam>stimTable.mat` for every animal in the manifest (`_anmlROI_FRAtable.mat` for FRA) | §11 |
| group plotters | `<Fam>_Group<g>.mat` | `aggregateStimGroup` |

**To resume at §10** an animal needs exactly three things: `_tifFileList.mat`, `_moCorrROI_*.mat`, and its `_Pulses.mat` files. `NoRMCorred/` and `FISSAoutput/` are *inputs to* §8–9 and are not read again afterwards — the fluorescence traces they produced already live inside `_tifFileList.mat`. Copying an animal without `NoRMCorred/` is therefore fine **if** `_tifFileList.mat` and `_moCorrROI_*.mat` came with it; copying only tifs and `_Pulses.mat` is not, because §4–5 (ROI drawing) is interactive and §8 (FISSA) is a Python step — unless you run the [headless path](#1c-headless-run--processanimal2pheadlessm), which automates both.

### Validate and stamp

`aggregateStimGroup` does two things beyond concatenating, and both are what make a group file trustworthy later.

**Validate** — `validateStimGroup` checks a table against its family's contract: required columns present, no duplicate `(animal, roiID, stim)` rows, one common dF/F time axis, trace widths matching that axis, scalar peak/significance columns, uniform frame rate, and (optionally) column order and time axis against a reference. It returns a report of problems graded `error` / `warning` / `info`; it never modifies anything. Each per-animal table is validated *before* concatenation and the assembled group *after*, so a bad animal is refused rather than diluted into a group. `loadStimGroup` re-runs it on every load, and it can be called standalone on any group file.

**Stamp** — a `groupInfo` struct saved alongside the table, recording:

| field | why it matters |
|---|---|
| `convention` | the analysis parameters in force (`tBasePT`, `pkPTsigSD`, …). Groups built under different conventions are not comparable; `loadStimGroup` warns when a stamp no longer matches `stimGroupSpec` |
| `animals`, `perAnimal`, `nCells` | membership and per-animal contribution, so a group dominated by one animal is visible |
| `manifest` | the membership request, cross-checked on load against the `.json` on disk |
| `sourceFiles` | which per-animal files went in, with size and timestamp |
| `timeVar`, `timeAxis` | the axis the traces share |
| `createdBy`, `created` | script and git SHA, and when |

The stamp is the direct answer to how CGC Groups C and D ended up built by two different script versions without anyone noticing: nothing recorded which recipe had produced a file. A group file without `groupInfo` is treated as legacy, not invalid.

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

*Step-by-step instructions: [usage.md § Path A](usage.md#path-a--manual-rois).*

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
| 11 | Per-stim analysis | `processAnimalStimFamilies` runs `processBPN2P` / `processCGC` / `processFRA` for whichever families this animal has, adding dF/F and peak responses and writing the **processed** table each group step reads. Each family is gated on the input it actually needs: a `_raw` table, or — for FRA, which has none — a `_tifFileList.mat` carrying map tifs. Set `runPerStimProcessing = false` to do it by hand with edited parameters |

**Output:** `tifFileList` struct where `tifFileList.stim(n).SCALEDfissaFroi` is an `nROI × nFrames` array of motion- and neuropil-corrected fluorescence for the nth stimulus tif. Equivalent `.map` field for FRA/BF mapping tifs.

---

### 1a. Tif compression — why motion-corrected tifs are written uncompressed

ScanImage writes some acquisitions **LZW-compressed**. `writeMoCorTifs` copies
the source acquisition's tags onto its output, and it used to copy
`Compression` along with them — so the motion-corrected tifs inherited LZW.
MATLAB reads those back correctly, which is why it went unnoticed; no other
reader does, and each failure surfaces far from the cause:

| Consumer | Symptom |
|---|---|
| `ScanImageTiffReader` (used by `flattenTif`) | **garbled pixels, no error** — a mean of int16 data in the range 20–1700 came back as ±8000 |
| FISSA / `tifffile` | `COMPRESSION.LZW requires the 'imagecodecs' package` |

`writeMoCorTifs` now sets `Compression.None` explicitly and no longer inherits
the tag, matching `writeTifWithHeader`. `flattenTif` additionally falls back to
MATLAB's `Tiff` class when `ScanImageTiffReader`'s frame 1 matches neither
orientation, rather than assuming it must be transposed.

**`NoRMCorred/` folders written before this fix still contain LZW tifs.**
`imagecodecs` is installed in the `env_fissa` conda environment, so FISSA reads
them; `ScanImageTiffReader` still cannot, so `flattenTif` pays a slower
streaming re-read. Regenerate those folders if you want the fast path.

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

### 1c. Headless run — `processAnimal2Pheadless.m`

*Step-by-step instructions: [usage.md § Path B](usage.md#path-b--auto-rois-headless).*

A **second entry point** that runs the same stages unattended: every dialog
becomes a config field, and `TIFcatROIgui` (§4–5) is replaced by Cellpose
segmentation. `processAnimal2P.m` remains the path for hand-drawn ROIs and for
anything you want to watch — this does not replace it, and both write the same
artefacts, so you can start an animal in one and finish it in the other.

```matlab
% one animal, everything
out = processAnimal2Pheadless('/data/TO0003');

% a cohort, resuming past motion correction
for a = ["TO0006","TO0007"]
    processAnimal2Pheadless(fullfile(root,a), 'stages',[4 11], ...
        'treatmentName','ZX1', 'preTifs','_0003[4-9]_');
end
```

Stage numbers match `processAnimal2P`'s sections (5 is folded into 4). **Each
stage is skipped when its artefact already exists** unless `'overwrite',true`,
so re-running to re-tune segmentation never repeats motion correction;
`'stages'` takes a subset or a `[from to]` range.

#### Every dialog, as a config field

`headlessConfig` documents all of them; the ones you are most likely to set:

| Interactive prompt | Config field | Default |
|---|---|---|
| `uigetdir` / animal ID | `dataPath`, `animal` | animal from the folder name |
| treatment name + pre-tif picker | `treatmentName`, `preTifs` | `''` → every tif `'none'` |
| FRA map tif picker | `mapTifs` | `'auto'` — see below |
| alternate-zoom prompt | `altZoomPolicy` | `'drop'` (the interactive default) |
| treatment filters | `condFilters` | derived from the treatments present |
| functional channel | `funcChan` | `2` |
| **ROI drawing** | `roi` (→ `cellposeROIset`) | consensus, `minVotes` 2, `dilatePx` 2 |
| FISSA scale factor | `fissaScaleFactor` | `0.8` |
| — (manual Python step) | `fissaCmd` | runs the repo script in conda env `env_fissa` |

`preTifs` and `mapTifs` accept **indices, a logical mask, file names, or a
regular expression**. A selector that matches *nothing* is an error, not an
empty group — a typo'd regexp would otherwise relabel a whole session as
post-treatment silently.

**`mapTifs` defaults to reading the pulse files**, not file size: a tif is an
FRA map when its `pulseSet` contains `'map'` (case-insensitive), the same
`contains(pulseSet,…)` idiom `stimParam2ROI` uses for BPN, PTinContrast and
spont. The `'bytes'` rule is available as a fallback for animals with no pulse
files, but it is only the *hint* the interactive dialog prints before a human
chooses: on TO0003 it flags all 8 real map tifs **and 5 stim tifs**, because
the 85.6 MB BPN tifs sit just under the 94.2 MB maps — which silently drops
the whole BPN family from the run.

#### ROI detection — `cellposeROIset`

Two things have to be true of the result: the masks must be cells, and an ID
must mean the *same* cell in every treatment condition. A human keeps the
second true by eye; a segmenter cannot, so one shared ROI set is segmented in a
common reference frame and placed into each condition
(`registerConditionMeans`, translation only, refusing a shift above
NoRMCorre's own `max_shift`). §6 `intersectROIfiles` then becomes an invariant
check rather than the thing doing the matching.

Detection is a **vote across per-tif segmentations** (`consensusROIsets`), not
one pass over the session mean. Cellpose-SAM's detection here depends on the
granular noise texture of a single tif's projection; averaging tifs removes it
and detection collapses non-monotonically (TO0003: 1 tif → 33 cells, 2 → 19,
4 → 0, all 29 → 0), while each individual tif agrees well with hand-drawn ROIs.
Voting keeps the regime where the segmenter works. `minVotes` spans the obvious
rules — `1` is the union, `numel(tifs)` the strict inner join (which keeps
nothing here, being governed by the worst tif) — and every cell carries a
**detection count**, which a single pass cannot provide.

Calibrated against TO0003's 18 hand-drawn ROIs: at `minVotes` 2 and
`dilatePx` 2 it recovers **18/18** at median IoU 0.71 with mask areas 0.82× the
hand-drawn convention, finding 102 ROIs in total. `diameter` must be set
explicitly (default 15) — cpsam's automatic sizing fails outright on a smoothed
mean. `adhoc/calibrateCellposeROI.m` re-runs that calibration on any animal
that has a hand-drawn set.

Cellpose runs **locally by default** in conda env `suite2p`
(`~/.cellpose/models/cpsam`); `cellposeSegment`'s `'backend','podman-exec'`
targets a container instead, for machines where that is the only cellpose.

#### Known limits

- The multi-condition registration path has unit coverage with synthetic
  shifts but has not yet run on a real pre/post animal — run the first one with
  `'verbose',true` and sanity-check the reported shift.
- `stimParam2ROI` and the `process*` scripts re-derive the animal ID from the
  **folder name**, so a folder that does not contain it fails at stage 10;
  `headlessConfig` warns about this up front.
- `type` on a Cellpose ROI is `'cellpose'`, which `TIFcatROIgui`'s load path
  does not yet understand — you cannot currently open a headless ROI set in the
  GUI for manual touch-up.

---

### 2. Per-stimulus analyses — `stimulusSpecific/`

**`processAnimal2P` section 11 runs these for you** — this section is the reference for what they do, and for running one by hand with non-default parameters. Section 10 (`stimParam2ROI`) writes one stim-aligned `_raw` table per stimulus family (`<animal>_anmlROI_BPNstimTable_raw.mat`, `<animal>_anmlROI_CGCstimTable_raw.mat`). Each script below loads its family's `_raw` table, adds dF/F and peak metrics, plots, and writes the processed bundle to the same name **without** the `_raw` suffix. Each resolves `dataPath` via `uigetdir` if not already in the workspace.

**`_raw` two-stage convention (BPN and CGC).** Keeping the raw and processed artifacts distinct means re-running a `process*` script never mutates its own input, and an animal that has not been processed yet has no processed file for `aggregateStimGroup` to pick up silently. Downstream consumers read the processed file.

**FRA is the exception:** it has no `_raw` stage, because pure-tone mapping trials are never stim-aligned by `stimParam2ROI` — `processFRA` reads `_tifFileList.mat` directly. Its `stimGroupSpec` entry therefore has an empty `suffixRaw`, and anything that gates on a `_raw` file must branch on that (`processAnimalStimFamilies` does).

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

#### `processFRA.m` — frequency response area (FRA) mapping, single animal

**The one family with no `_raw` stage.** Pure-tone mapping trials are not stim-aligned by `stimParam2ROI`; `processFRA` reads `<animal>_tifFileList.mat` directly (its `.map` field) and the `*_Pulses.mat` files beside it, and delegates to `FRAmap`. It writes two artifacts:

| File | Contains | Read by |
|---|---|---|
| `<animal>_FRAmap.mat` | `FRAmap` struct: `dBFreqMap` (nDB × nFreq cell), `CellPkRespLinDBfreq`, `CellSigPkLinDBfreq`, `BFuDB`, `dPrime`, `params` | `plotFRAmap`, `compileAnmlFRA` |
| `<animal>_anmlROI_FRAtable.mat` | `anmlROIbyFRA`: long-form, one row per (ROI, frequency, level) | the group path |

The table exists so FRA can use the **same generic group machinery** as BPN and CGC — `aggregateStimGroup`, `validateStimGroup`, `groupN` all work on tables — rather than needing a parallel aggregator. `FRAmap2table` does the conversion and refuses a struct with no `.params` field (see below).

**Significance convention.** One test per (ROI, frequency, level), on the **trial-averaged** onset-aligned dF/F trace — the same rule `processBPN2P` and `processCGC` follow, and what `pkFcalc` requires. `sigPkDFF` is a logical mask; `CellSigPkLinDBfreq` holds the peak where significant and `NaN` elsewhere.

**Baseline.** Pure-tone pulses are contiguous segments of one recording, so the baseline reaches **backwards across the pulse boundary** into the preceding pulse's silence rather than using only the 2–3 pre-onset frames inside the window. Its length is computed automatically as the longest that still clears the previous tone's response (7 frames on the TOMT data, up from 3), and never shorter than the within-pulse baseline. `pkFcalc` gained an optional `baseIDX` argument for this, defaulting to its previous window so **BPN and CGC are bit-identical**.

**Onset** resolves via `find(t >= onset, 1)`, not `onset*fs` — frame *k* spans `(k-1)/fs`, so the latter lands one frame early.

> **`params` is the version marker.** A saved `_FRAmap.mat` **without** a `.params` field predates this convention: it holds per-trial significance and a peak-**squared** response map. Regenerate it rather than comparing against it. `params` records `sigMethod`, `baseline`, `onsetCol`, `nBaselineFrames` and the analysis parameters.

---

### 3. Condition groups — `aggregateStimGroup` + the group plotters

*Step-by-step instructions: [usage.md § After either path](usage.md#after-either-path--condition-groups). Note `requireRun` — see [usage.md § Provenance](usage.md#provenance).*

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
  §11    processFRA (no _raw stage)            plotFRAgroup
         → _FRAmap.mat + _anmlROI_FRAtable.mat
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
plotFRAgroup('FRA_GroupD.mat');    % BF, threshold, bandwidth, FRA colormaps
```

`processRLF.m` and `processFRAgroup.m` are the interactive entry points for the BPN and FRA plots — each prompts for a group file and calls its plotter.

#### FRA group plots — `plotFRAgroup`

Four population panels, built from per-cell tuning descriptors that `tableFRAmetrics` derives from the group table (the FRA analogue of `tableRLF`):

| Panel | Shows |
|---|---|
| `bf` | best-frequency distribution, occurrence (%) |
| `threshold` | threshold distribution, occurrence (%) |
| `bw` | bandwidth at threshold + one level, plus bandwidth-vs-level |
| `fra` | group-mean FRA colormap + a grid of single-cell receptive fields |

Two measurement limits are surfaced in the output rather than left implicit, because both are properties of how the data was acquired, not of the code:

- **Bandwidth is BW20, not BW10.** Levels are 30/50/70 dB SPL, 20 dB apart, so threshold + 10 dB is never a level that was presented. Bandwidth is measured at the next sampled level and labelled `BW20` everywhere. `bwByLevel` additionally gives bandwidth at every level.
- **Every run prints a noise floor.** `tableFRAmetrics` re-runs `FRAmap`'s significance test on a *silent* late window of the same width, where no tone-evoked response can exist, and reports `realRate / shamRate`. A ratio near 1 means the significance mask carries little stimulus information, and `plotFRAgroup` warns below 1.2. On the TOMT groups it is **A 1.73, B 0.84, C 1.63, D 1.08** — so the threshold and bandwidth panels for groups B and D should not be over-interpreted.

`threshold` also requires a contiguous band of `minBand` frequencies (default 2) at a level before that level counts. With 28 frequencies tested per level, isolated significant frequencies are common by chance; without the band requirement every cell's threshold collapsed to the lowest level.

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

For **additional CGC figures** this path does not produce — the FRA/best-frequency join, pre→post treatment comparisons, PV/SOM cell-type splits, bandwidth tuning and the manuscript figures — refer to [`stimulusSpecific/extraCGC/`](#appendix--additional-cgc-figures-stimulusspecificextracgc). It is a separate route that this one does not replace.

---

## External dependencies

| Dependency | Language | Purpose |
|------------|----------|---------|
| [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) | MATLAB | Non-rigid motion correction |
| [FISSA](https://github.com/rochefort-lab/fissa) | Python | Neuropil signal separation (conda env `env_fissa`; needs `imagecodecs` for pre-fix LZW tifs, see §1a) |
| [Cellpose](https://github.com/MouseLand/cellpose) (`cpsam`) | Python | ROI segmentation for the headless path only (conda env `suite2p`, model at `~/.cellpose/models/cpsam`) |
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
  processFRA.m        — per-animal FRA mapping, straight from tifFileList
                        (no _raw stage); writes _FRAmap.mat + _anmlROI_FRAtable.mat
  plotBPNgroup.m      — BPN group plots: RLF, dF/F re level, peak re level
  plotCGCgroup.m      — CGC group plots: contrast traces, low-vs-high, paired bar
  plotFRAgroup.m      — FRA group plots: BF, threshold, BW20, FRA colormaps
  plotCellTrials.m    — every repetition of one cell, for chasing outliers
  processRLF.m        — interactive entry point for the BPN group plots
  processFRAgroup.m   — interactive entry point for the FRA group plots
  runFRA.m            — non-interactive FRA for one animal (function form of
                        processFRA; used to recompute FRA on its own)
  extraCGC/           — the cohort-table path (see the appendix):
                        compileCohortData.m, plotCohortData.m

helperFcns/
  tif/                — ScanImage tif reading (readSCIMtif, justLoadTif),
                        channel splitting (splitTifChans), writing (writeMoCorTifs)
  dFF/                — dF/F computation (dFoFcalc), peak response detection
                        (pkFcalc; optional baseIDX sets the significance
                        baseline — FRA passes a strictly pre-onset window)
  FRA/                — FRA map construction (FRAmap), BF extraction (anmlFRA2BF),
                        d-prime (anmlFRA2dPrime), struct→long-form table
                        (FRAmap2table), per-cell BF/threshold/bandwidth and the
                        sham noise-floor control (tableFRAmetrics)
  RLF/                — rate-level functions (cellRLF, tableRLF) and plotting (plotRLF)
  cohort/             — group-size-agnostic primitives shared by both family
                        plotters: counts (groupN), trace/value stacking
                        (gatherCellTraces, gatherCellValues), between-cell
                        mean/SEM (cohortMeanSEM), guarded statistics
                        (cohortStat), figure n-stamps (annotateN)
  dataOrg/            — artifact layout resolution (animalPaths) and optional
                        migration (migrateAnimalArtifacts);
                        FISSA output parsing, tifFileList assembly, stimulus table
                        builders (stimParam2ROI, combineDiffOnset); condition-group
                        aggregation (stimGroupSpec, aggregateStimGroup,
                        validateStimGroup, loadStimGroup); per-animal
                        per-family driver (processAnimalStimFamilies);
                        headless run configuration (headlessConfig)
  ROI/                — ROI mask ↔ polygon conversion, raw F extraction from masks;
                        256×128 ROI reuse (remapROItoAcq, remapROIfile);
                        and the headless detection path — Cellpose CLI wrapper
                        (cellposeSegment), label image → moCorROI
                        (labelImg2moCorROI), per-tif vote (consensusROIsets),
                        cross-condition alignment (registerConditionMeans),
                        the orchestrator (cellposeROIset), plus IoU comparison
                        of two ROI sets (compareROIsets) and streamed mean
                        projection (condMeanImg)
  plotting/           — SEM shaded plots (fillSEMplot), regression plots (regPlot),
                        paper-style formatting, figure export
  sound/              — Ephus .signal file inspection (inspectSignalObject)
  general/            — utility functions (SEMcalc, sigDiffCalc, scaleZeroToOne, zero2nan)
```

---

## Appendix — additional CGC figures (`stimulusSpecific/extraCGC/`)

**Where to look for CGC analyses the condition-group path does not cover.** A
second, older route: `compileCohortData` builds one flat `Tinput` table across
animals, and `plotCohortData` produces figures from it. Every section is
CGC-specific. It is **not superseded** — it is the only place several analyses
exist — but it is not the default either.

**Use the condition-group path ([section 3](#3-condition-groups--aggregatestimgroup--the-group-plotters)) for**
comparing treatment groups: RLF, dF/F re level, contrast traces, low-vs-high
peaks, significant-cell counts. It consumes the per-animal *processed* tables,
so every group inherits one analysis convention and carries a provenance stamp.

**Refer to `extraCGC/` for** figures that path has no equivalent of:

| Analysis | `plotCohortData.m` section |
|---|---|
| PT response ratio vs. PT onset | `PN \| PRE: PT onset` |
| sustained DRC response, pre and post | `PN \| PRE: dFF_DRC`, `PN \| POST: dFF_DRC` |
| peak scatter + traces re contrast | `PN \| PRE: dFF_PT \| scatterplot and traces re contrast` |
| bandwidth / octave tuning | `PN \| PRE: dFF_PT \| BW, OCT` |
| **pre → post treatment** | `PN \| PRE --> POST: dFF_PT` |
| DRC offset *(supplemental)* | `PN \| PRE: DRC OFFSET` |
| contrast-change dF/F *(supplemental)* | `PN \| PRE: DRC contrast change dFF` |
| ZX1 re tuning *(supplemental)* | `PN \| PRE --> POST: ZX1 re tuning` |
| ZnT3 pre → post response ratio | `ZnT3?? \| PRE --> POST: dFF_PT RATIO` |
| **PV & SOM cell types** | `PV & SOM \| PRE OR POST: ...` (two sections) |

Each block is wrapped in `%{...%}` and run independently. It also uniquely
provides the **FRA / best-frequency join into the CGC table** (`compileAnmlFRA` →
`BFuDB`, keeping only tone-responsive cells with `dPrime > 0`): the group path
now has its own FRA route ([`plotFRAgroup`](#fra-group-plots--plotfragroup)), but
that produces FRA figures in their own right and does not join BF onto CGC
responses. `matlabPAC_CGCplot/plotDataTable.m` is a 5-section curated extract of
these 14, for the manuscript figures.

> `compileAnmlFRA` calls `FRAmap` itself rather than loading `_FRAmap.mat`, so it
> picks up the current convention automatically — but any cohort table built
> before 2026-09-01 carries per-trial significance and a peak-squared response
> map. Rebuild rather than compare.

### Things to know before using it

- **The dF/F method matches `processCGC` as of 2026-09-01.** `dFF_PT` is the
  additive `dFF_DRC − mean(dFF_DRC over tBasePT)`, the peak is taken on the
  `t >= tBasePT(1)` crop so `pkFcalc`'s baseline is the pre-tone window, and one
  `pkPTframeBin` is used throughout. It previously used a *divisive*
  `(F − F0_PT)/F0_PT` off raw fluorescence, with two different frame bins in
  different sections. **Figures made before that change are not comparable with
  figures made after** — on the CaMKII sample cohort, peak dF/F median moved
  0.0542 → 0.0611 (r = 0.95 old vs new) and the per-trial significance rate
  48.8% → 69.0%.
- The constants (`tBaseDRC`, `tBasePT`, `pkPTsigSD`, `pkPTframeBin`) are still
  declared locally rather than read from `stimGroupSpec`, so the *method* agrees
  with the group path but the values could still drift apart.
- dF/F is recomputed at plot time from the canonical `F` column, not reused from
  the per-animal processed tables.
- `compileCohortData.m` re-indexes `animal` and `roiID` to sequential integers,
  so cell identity is **not** comparable with group files.
- Both scripts prompt rather than carrying hardcoded paths. `plotCohortData.m`
  asks for the `*_dataTable.mat` and derives the matching `*_params.mat` from
  the same prefix; preset `cohortDataFile` in the workspace to skip the dialog.
  Missing parameters are filled per field, so an older params file still works.
- `compileCohortData.m`'s call into `compileAnmlROItables` was **broken until
  2026-08-31** (wrong argument order), so anything produced before then came
  from an older revision.

### A1. Compile the cohort table — `stimulusSpecific/extraCGC/compileCohortData.m`

Builds the flat `Tinput` table that A2 and `matlabPAC_CGCplot/plotDataTable.m`
consume. Edit the parameters block at the top (`cohortName`, `family`).

- `compileAnmlFRA` — collects FRA maps and BF estimates across animals
- `compileAnmlROItables(family, params)` — collects stimulus/response tables per animal, validating each against its family contract before reducing columns
- Joins BF data into the response table; keeps only tone-responsive cells (`dPrime > 0`)
- Saves `<cohortName>_dataTable.mat` and `<cohortName>_params.mat`

### A2. Plot the cohort table — `stimulusSpecific/extraCGC/plotCohortData.m`

Produces the figures listed above. Prompts for the data table, then run whichever
`%{...%}` section you need.

- `dFF` section computes `dFF_DRC` then the additive `dFF_PT`, and peak PT responses via `pkFcalc`
- Statistical testing via `sigDiffCalc` (parametric or Wilcoxon) with a permutation-test fallback
- Figures save to `params.figSaveDir` (defaults to the data table's folder) when a section sets `figSave = true`
