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
| when to use it | few animals; a field you want to judge yourself; anything you'll compare to previously published counts | cohorts; re-runs; anything where consistency between animals matters more than your judgement on any one |

Both write the **same artefacts under the same names**, so you can start an
animal in one and finish it in the other. What they cannot do is *coexist* —
see [Both paths write the same filenames](#both-paths-write-the-same-filenames).

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

The Python step at §8:

```bash
conda run -n env_fissa python FISSAviaMatlab_prePostTreatment.py /path/to/TO0003
```

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
| `treatmentName`, `preTifs` | `''` → all `'none'` | `preTifs` takes indices, a logical mask, file names, or a regexp |
| `mapTifs` | `'auto'` | reads each tif's `pulseSet` for `map` (case-insensitive). **Do not use `'bytes'`** unless there are no pulse files — it mis-splits map from stim |
| `roi` | see below | forwarded to `cellposeROIset` |
| `excludeNeg` | `true` | movement screening; BPN only (CGC is single-pulse and never screened) |
| `runLabel` | `''` | stamps each processed table; see [provenance](#provenance) |
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

## Both paths write the same filenames

`<animal>_moCorrROI_*.mat`, `_tifFileList.mat`, `_anmlROI_*.mat` and the rest
are **not** namespaced by ROI source. A headless run **overwrites hand-drawn
ROI files**.

If you want to keep both, or compare them:

1. Back up the animal folder's `.mat` artefacts first.
2. Run one path to completion, then copy its artefacts somewhere per-run.
3. Clean the animal folder and run the other.

Treat the animal folder as scratch space, not as the record of any particular
run. `etc/runTOMTpipeline.m` is a worked example of exactly this discipline
(clean → run → harvest to `runArtifacts/<runName>/`), if you need to do it over
a cohort.
