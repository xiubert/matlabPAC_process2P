# Usage — running an animal through the pipeline

Two ways to get ROIs, one pipeline behind them. This is the task-oriented
guide: what to type, in what order, and what to check. For what each function
*is*, see [README.md](README.md).

| | **Manual ROIs** | **Auto ROIs (Cellpose)** |
|---|---|---|
| entry point | `processAnimal2P.m` (script, run section by section) | `processAnimal2Pheadless.m` (function, one call) |
| ROIs come from | you, drawing in `TIFcatROIgui` | Cellpose, voted across per-tif segmentations |
| needs a display | yes | no — runs under `matlab -batch` |
| cells found | what you draw | ~3× more on the TOMT cohort, at the same response quality |
| ROI review | you saw them as you drew | `<animal>_ROIoverlay_<cond>.png`, written automatically |
| where artifacts go | `analysis/<run>/` — set `runName` at the top of the script | `analysis/<run>/` — set `'run'`, or let it be derived |
| when to use it | few animals; a field you want to judge yourself; anything you'll compare to previously published counts | cohorts; re-runs; anything where consistency between animals matters more than your judgement on any one |

Both write the **same artefacts under the same names**, so you can start an
animal in one and finish it in the other. To keep two analyses side by side,
give each a run name — see
[Where artifacts are written](#where-artifacts-are-written).

---

## Before either path

```matlab
addpath(genpath('/path/to/matlabPAC_process2P'))
```

| need | where |
|---|---|
| NoRMCorre | on the MATLAB path |
| FISSA | conda env `env_fissa` (also needs `imagecodecs` for pre-2026-09 `NoRMCorred/` folders, see [README §1a](README.md)) |
| Cellpose + `cpsam` | conda env `suite2p`, model at `~/.cellpose/models/cpsam` — **auto path only** |

The animal folder must contain the raw `*.tif` and their `*_Pulses.mat`, and
**the folder name must contain the animal ID** (`TO0003`, `AA0067`) — several
downstream functions re-derive it from the path.

---

## Path A — manual ROIs

Open `processAnimal2P.m` and run it section by section. It prompts for
everything.

**Set `runName` near the top** (just after the animal ID). It names the folder
this analysis writes to — `analysis/<runName>/` — exactly as the headless
path's `run` does, so both paths produce the same structure. It defaults to
`manualROI_<date>`; set it to `''` for the old flat layout. Section 8 prints
the FISSA command with this run's folders already filled in.

| § | what you do |
|---|---|
| 1 | pick the treatment name, the pre-treatment tifs, and the FRA map tifs |
| 2 | give the treatment filters that split motion-correction groups |
| 3 | NoRMCorre runs (slow; the one step you never want to repeat) |
| **4–5** | **draw ROIs in `TIFcatROIgui`, once per 256×256 condition, then run §5 to save** |
| 5b | automatic — remaps ROIs onto any 256×128 condition |
| 6 | `intersectROIfiles` keeps only cells present in every condition |
| 7 | raw F per ROI per tif |
| **8** | **stop. Run FISSA outside MATLAB (below), then continue** |
| 9 | parse FISSA output → `_tifFileList.mat` |
| 10–11 | stimulus alignment and per-family dF/F + peak responses |

The Python step at §8 — **run the section and it prints the command**, with
this run's ROI and output folders filled in:

```bash
conda run -n env_fissa python FISSAviaMatlab_prePostTreatment.py /data/TO0003 \
    --roi-dir   /data/TO0003/analysis/manualROI_20260903 \
    --tiff-folder /data/TO0003/NoRMCorred \
    --out-dir   /data/TO0003/analysis/manualROI_20260903/FISSAoutput
```

The tifs are shared by every analysis of the animal; the ROIs and the output
belong to this run.

**Resuming.** §3b reloads the motion-corrected data without re-running
NoRMCorre — use it when you come back to draw ROIs for another cell type.

---

## Path B — auto ROIs (headless)

```matlab
out = processAnimal2Pheadless('/data/TO0003');
```

That runs every stage, including FISSA. Nothing is prompted. Common variants:

```matlab
% pre/post treatment, pre-tifs chosen by a regexp on the file names
processAnimal2Pheadless('/data/TO0006', ...
    'treatmentName','ZX1', 'preTifs','_0003[4-9]_');

% a cohort, resuming past motion correction
for a = ["TO0006","TO0007"]
    processAnimal2Pheadless(fullfile(root,a), 'stages',[4 11]);
end

% re-tune segmentation only; stage 3 is skipped because its output exists
processAnimal2Pheadless('/data/TO0003', 'stages',[4 11], 'overwrite',true, ...
    'roi',struct('dilatePx',1,'minVotes',3));
```

Stage numbers match `processAnimal2P`'s sections. **A stage is skipped when its
artefact already exists** unless `'overwrite',true`, so re-running never
repeats motion correction.

Every option is documented in `help headlessConfig`. The ones you will actually
set:

| option | default | note |
|---|---|---|
| `treatmentName`, `preTifs` | `''` → all `'none'` | which tifs are pre-treatment — see [Choosing tifs](#choosing-tifs-pretifs-and-maptifs) |
| `mapTifs` | `'auto'` | reads each tif's `pulseSet` for `map` (case-insensitive). **Do not use `'bytes'`** unless there are no pulse files — it mis-splits map from stim |
| `roi` | see below | forwarded to `cellposeROIset` |
| `excludeNeg` | `true` | movement screening; BPN only (CGC is single-pulse and never screened) |
| `run` | derived | folder under `analysis/`, e.g. `cellpose_20260903`; also becomes the provenance stamp — see [Where artifacts are written](#where-artifacts-are-written) |
| `isolateRun` | `true` | `false` writes into the animal folder instead (the pre-2026-09 flat layout) |
| `runLabel` | follows `run` | the provenance stamp, when you want it to differ from the folder name |
| `stages`, `overwrite` | `1:11`, `false` | `stages` takes a subset or a `[from to]` range |

A selector that matches **nothing** is an error, not an empty group — a typo'd
regexp would otherwise silently relabel a whole session as post-treatment.

### ROI settings worth knowing

Defaults, calibrated against hand-drawn ROIs on this data:

```matlab
'roi', struct('mode','consensus', 'minVotes',2, 'dilatePx',1, ...
              'autoDiameter',true, 'edgeMarginPx',3)
```

- **`autoDiameter`** picks Cellpose's `--diameter` per condition. Leave it on.
  The wrong diameter finds *nothing at all* rather than something worse — on
  one animal, 10 found 47 cells while 12, 15, 18, 20 and 25 each found 0 or 1.
- **`dilatePx = 1`** makes masks match the hand-drawn convention (area ratio
  1.03 across 13 animals). `2` gives masks 23 % too large, which inflates
  neuropil contamination.
- **`minVotes = 2`** means a cell must appear in ≥2 tifs. `1` is the union,
  `numel(tifs)` the strict inner join — the join keeps nothing on real data,
  because the worst tif governs.

### Choosing tifs (`preTifs` and `mapTifs`)

`treatmentName` names the treatment; `preTifs` says which tifs came **before**
it. Everything not selected becomes post-treatment. So with
`'treatmentName','ZX1'`, selected tifs are labelled `preZX1` and the rest
`postZX1`. Leave `treatmentName` empty and every tif is `'none'` — `preTifs` is
then ignored.

`mapTifs` takes the **same four forms**, but you rarely need them: its default
reads the pulse files, which is authoritative.

All examples below are real, run against `TO0007` (33 tifs, named
`TO0007AAAA_00031_00001.tif` … `TO0007AAAA_00070_00001.tif`).

**1. Indices** — positions in the tif list, in `dir()` order (alphabetical,
which for these names is acquisition order):

```matlab
'preTifs', 1:5          % the first five tifs -> _00031 … _00035
'preTifs', [1 3 7]      % specific positions
```

Positions, *not* acquisition numbers — `1:5` means the first five files, not
tifs 1–5. An out-of-range index is an error.

**2. A regular expression** — matched against the **file name** (not the path).
Usually the clearest option, because it keys on the acquisition number:

```matlab
'preTifs', '_0003[1-4]_'      % _00031 _00032 _00033 _00034   (4 tifs)
'preTifs', '_000(31|45)_'     % _00031 and _00045 only        (2 tifs)
'preTifs', '_0003'            % every tif whose number starts 0003
```

The `_`…`_` wrapping matters: `'_0003[1-4]_'` anchors to the acquisition-number
field, whereas a bare `'3'` would match almost everything.

**3. File names** — explicit and unambiguous, good when the split is irregular:

```matlab
'preTifs', {'TO0007AAAA_00031_00001.tif','TO0007AAAA_00032_00001.tif'}
```

A name that isn't in the folder is an error, so a typo can't silently shrink
the group.

**4. A logical mask** — one element per tif, for when you compute the split:

```matlab
m = false(33,1);  m([2 4 6]) = true;
'preTifs', m
```

#### Check before you commit to a run

`headlessConfig` resolves the selectors without writing anything, so you can
see exactly what you'd get:

```matlab
cfg = headlessConfig('/data/TO0007','treatmentName','ZX1','preTifs','_0003[1-4]_');
cfg.tifNames(cfg.preIDX)     % the pre-treatment tifs
cfg.tifNames(cfg.mapIDX)     % the FRA map tifs
nnz(cfg.preIDX)              % how many
```

**A selector that matches nothing is an error**
(`headlessConfig:badSelector`), not an empty group — a typo'd regexp would
otherwise relabel the whole session as post-treatment and nothing would say so.
The one case that only *warns* is `treatmentName` set with `preTifs` left
empty, since "everything is post-treatment" is occasionally what you mean.

### Always look at the overlay

Each condition gets `<animal>_ROIoverlay_<cond>.png`: every ROI outlined on the
image it was segmented from, **coloured by how many tifs detected it**. This is
the headless path's substitute for having watched someone draw. Regenerate for
a finished animal with:

```matlab
plotROIoverlay('/data/TO0003','TO0003')
```

**What to look for.** Outlines sitting on somata, and a plausible vote spread.
If a whole region is bare, or nearly every cell has the minimum votes, the
segmentation is weak — check the diameter it chose (printed during the run, and
in `cellposeParams` inside the ROI file).

### When it goes wrong

| symptom | cause | fix |
|---|---|---|
| `Every ROI set is empty` | wrong diameter, or a genuinely low-contrast field | confirm `autoDiameter` is on; widen `diameterLadder`; if the field is blurry, use manual ROIs for that animal |
| very few ROIs vs. neighbours | few tifs, so `minVotes=2` is a strict join | lower `minVotes`, or accept it |
| FISSA: `COMPRESSION.LZW requires 'imagecodecs'` | `NoRMCorred/` predates the compression fix | install `imagecodecs` in `env_fissa`, or regenerate the folder |
| `lacks moCorTifNames` | ROI file predates the grouped FISSA driver | `ensureROIfileMeta` backfills it; stage 6 calls it automatically |

---

## After either path — condition groups

Identical for both. One manifest per group per family:

```matlab
manifest = struct('group','A', 'family','BPN', ...
                  'animals',{{'TO0001','TO0002','TO0003','TO0014'}}, ...
                  'cohortRoot','/data/cohort', ...
                  'outDir','/data/cohort/aggregate');
info = aggregateStimGroup(manifest);

plotBPNgroup(info.outFile)      % or plotCGCgroup / plotFRAgroup
```

Animals whose tables live in a subfolder (a `Region1/`, a separately
motion-corrected `BPN/`) need explicit paths instead of `animals` +
`cohortRoot`:

```matlab
manifest.files = {'/data/TO0004/Region1/TO0004_anmlROI_BPNstimTable.mat', ...};
```

---

## Provenance

Every run writes the same filenames into the same animal folder. If a stim
family *fails*, the previous run's table stays on disk, and aggregation has no
way to tell it apart — which silently produced a group file mixing hand-drawn
and segmented ROIs.

So: **pass a `runLabel`, and require it when aggregating.**

```matlab
processAnimal2Pheadless(dp, 'runLabel','cellpose_2026-09');
% ... or, running the families by hand:
processAnimalStimFamilies(dp, 'runLabel','cellpose_2026-09');

aggregateStimGroup(manifest, 'requireRun','cellpose_2026-09');
```

`aggregateStimGroup` then refuses any animal whose stamp names a different run,
or carries none at all, and records what it found in `groupInfo.sourceRuns`
either way.

⚠️ **Running a `process*` script by hand leaves an unstamped table.** That is
fine until you aggregate with `requireRun`, which will reject it. Go through
`processAnimalStimFamilies` with a `runLabel` if the result is destined for a
group.

---

## Where artifacts are written

Everything the pipeline touches falls into one of three groups, and they have
different lifetimes:

| | what | lifetime |
|---|---|---|
| **raw** | `*.tif`, `*_Pulses.mat`, `*_PulseParams.mat` | never written by the pipeline |
| **shared** | both legends, `NoRMCorred/*_NoRMCorre.tif`, `_NoRMCorreParams.mat` | one per animal — expensive, identical for every analysis |
| **per-run** | `_moCorrROI_*`, `FISSAoutput/`, `_moCorr_Tifs_Params`, `_tifFileList`, `_pulseLegend2P`, `_stimGroupIDX`, `_anmlROI_*`, `_FRAmap`, QC images | one per analysis |

Pass a **run name** and the per-run group goes into its own folder:

```matlab
processAnimal2Pheadless('/data/TO0007','run','cellpose_20260903');
```

```
TO0007/
  TO0007AAAA_00031_00001.tif           raw
  TO0007AAAA_00031_00001_Pulses.mat    raw
  TO0007_tifFileLegend.mat             shared
  TO0007_tifCondSplitLegend.mat        shared
  NoRMCorred/                          shared — never duplicated
    *_NoRMCorre.tif
    TO0007_NoRMCorreParams.mat
  analysis/
    handdrawn_20260903/
      TO0007_moCorrROI_all.mat
      TO0007_tifFileList.mat
      TO0007_anmlROI_BPNstimTable.mat
      TO0007_ROIoverlay_all.png
      FISSAoutput/
    cellpose_20260903/
      …the same set, independently…
```

Two analyses of the same animal now cannot touch each other. `NoRMCorred/`
stays shared because motion correction depends only on the acquisition — it
would be wasteful and slow to redo per run.

`animalPaths` is the single source of truth for this; no function builds an
artifact path by hand.

```matlab
P = animalPaths('/data/TO0007','run','cellpose_20260903');
P.artifacts   % .../TO0007/analysis/cellpose_20260903
P.fissaDir    % .../analysis/cellpose_20260903/FISSAoutput
P.moCorrDir   % .../TO0007/NoRMCorred                    (shared)
P.legend      % .../TO0007/TO0007_tifFileLegend.mat      (shared)
```

### Runs are isolated by default

Omit `run` and a name is **derived** from what the run is — the ROI source and
the date, e.g. `cellpose_20260903` or `savedROI_20260903`. You never get two
analyses on top of each other just because nobody named them.

To reproduce an old flat run in place, opt out explicitly:

```matlab
processAnimal2Pheadless(dp,'isolateRun',false);   % straight into the animal folder
```

That is the pre-2026-09 behaviour, and the one thing to know about it is what
it cannot do: a second run overwrites the first, silently. A headless run
overwrites hand-drawn `_moCorrROI_*.mat`; a re-run with different parameters
overwrites the previous tables; and if a stim family *fails*, its table is
simply left behind from the previous run — the dangerous case, because the
folder then holds a mixture and looks perfectly normal.

### If a run cannot find its inputs

A run only sees its own folder, so stage 10 will not pick up a `_tifFileList`
that another run wrote. The error says which runs do have it:

```
Stage 10 needs TO0003_tifFileList.mat, which run "cellpose_20260903" has not written
Other run(s) do have it: legacy
```

Either work in that run (`'run','legacy'`) or run the earlier stages so this
one produces its own. If the artifact is loose in the animal folder instead,
the folder predates the layout and the error points you at
`migrateAnimalArtifacts`.

### Run name, provenance and group files

The run name does double duty — it names the folder **and** stamps the tables,
so the two cannot drift apart:

```matlab
processAnimal2Pheadless(dp,'run','cellpose_20260903');       % runLabel follows
aggregateStimGroup(manifest,'requireRun','cellpose_20260903');
```

Point aggregation at a run with `tableDir`, and keep each run's group files
apart with `outDir`:

```matlab
manifest.tableDir = fullfile('analysis','cellpose_20260903');
manifest.outDir   = '/data/cohort/aggregate_cellpose_20260903';
```

### Migrating an existing animal

A folder written before this layout has its artifacts loose at the top level,
where a named run will not look. Migrate it once:

```matlab
migrateAnimalArtifacts('/data/TO0007')                % dry run, lists what would move
migrateAnimalArtifacts('/data/TO0007','apply',true)   % -> analysis/legacy/
```

It moves only the per-run group; raw, the legends and `NoRMCorred/*.tif` stay
put, and `NoRMCorred/FISSAoutput/` moves into the run. **Back up first** — it
moves rather than copies, so anything still expecting the flat layout will stop
finding its inputs.

Afterwards that run is addressable like any other: `'run','legacy'`.

### Comparing two ROI sources

With run names this is just two runs; nothing needs cleaning between them:

```matlab
processAnimal2Pheadless(dp,'run','handdrawn_20260903','stages',[6 11]);
processAnimal2Pheadless(dp,'run','cellpose_20260903');
```

Run the hand-drawn path first only if its ROI files live in the flat folder —
a headless run without a run name would overwrite them.
