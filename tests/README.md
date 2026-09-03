# Tests

```
matlab -batch "addpath('tests'); runTests"                 % everything
matlab -batch "addpath('tests'); runTests('unit')"         % no data drive needed
matlab -batch "addpath('tests'); runTests('Filter','FRA')" % one area
```

`runTests` prints one line per test and exits non-zero if anything fails, so it
works as-is in CI. Each test's own output is captured and replayed only when
that test fails; pass `'Verbose',true` to stream it instead.

## Layout

| Folder | Runs in the suite | Needs recorded data |
|---|---|---|
| `unit/` | yes | no |
| `integration/` | yes | yes — skips cleanly without it |
| `investigations/` | **no** | yes |

**`unit/`** — synthetic fixtures only. Runs anywhere the repo does, in about a
second. If one of these breaks, the cause is in the repo, not in the data.

**`integration/`** — reads recorded datasets from the lab drive. A test whose
fixture is missing reports `SKIP`, not `FAIL`, so the suite still means
something on a machine that has the repo but not the data.

**`investigations/`** — diagnostic scripts with **no assertions**. They print
what a function does; they do not check it. Kept because they are the record of
how a question was settled, not because they guard anything. `runTests` does not
run them — run them by hand when revisiting that behaviour.

## Fixtures

The recorded data lives outside the repo. `testConfig` resolves it, and two
environment variables override the defaults for a machine that mounts it
elsewhere:

| Variable | Default | Used by |
|---|---|---|
| `MATLABPAC_TESTDATA` | `/media/DATA/Ophys/Jinbo` | every integration test |
| `MATLABPAC_NORMCORRE` | `~/Documents/MATLAB/NoRMCorre` | `testAA0072_pipeline` only |

Under `MATLABPAC_TESTDATA`:

| `testConfig` field | Path | What it is |
|---|---|---|
| `.aggregateDir` | `TOMT/aggregate data` | cohort group files (`<Family>_Group<g>.mat`) plus the `bkmat/` pre-fix backups |
| `.animalsRoot` | `TOMT/animals` | per-animal folders `TO0001`…`TO0014` |
| `.exampleAnimalDir` | `TO0003` | the TO0003 working copy the BPN and remap tests read |
| `.testdata10HzDir` | `Testdata_10Hz/AA0072` | mixed 5 Hz / 10 Hz acquisition |

`.exampleAnimalDir` and `animalsRoot/TO0003` are two near-identical copies of
the same animal on disk. The BPN and remap tests were written against the first,
the FRA tests against the second. Both are kept so neither set silently starts
reading the other's artifacts.

## Writing a test

Start with the bootstrap — two levels up from a test file is `tests/`, which is
where `testConfig` lives:

```matlab
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
```

`testConfig` also puts the repo, `helperFcns` and `stimulusSpecific` on the
path, so no test needs its own `addpath` lines.

Then:

- **`requireFixture(path, label)`** for anything read off the data drive. It
  throws `matlabPAC:test:fixtureMissing`, which the runner reports as `SKIP`.
  Any other error is a real failure.
- **`testSandbox(name)`** for scratch space:
  ```matlab
  [tmp,tmpCleanup] = testSandbox('myTest');   %#ok<ASGLU>
  ```
  It returns an empty folder and an `onCleanup` that removes it — including
  when the test dies partway, which a trailing `rmdir` does not. Hold onto the
  second output or the folder is deleted immediately.

**Never write into a recorded animal folder.** Copy what you need into a
sandbox first. `testProcessBPN2P_saveAndPlots` shows the shape of this, and why
it matters: `processBPN2P` does a bare `load` of the `_raw` bundle, and that
bundle carries a `dataPath` variable of its own, so `dataPath` is silently
reset to the recording folder partway through the script. Assertions there use
a separately pinned `sandboxDir`, and the test fingerprints the real folder's
processed file to prove it was not touched.

Tests may be functions or scripts; the runner detects which. Scripts run inside
a disposable function workspace, so they cannot clobber the runner's own
variables.

## What is covered

**unit/**

| Test | Covers |
|---|---|
| `testResolveROIset` | membership-first ROI-set resolution when two conditions share an ROI count |
| `testFISSAgrouping` | per-ROI-count FISSA merge; a session mixing 256×256 and 256×128 keeps each tif's own group |
| `testLabelImg2moCorROI` | segmentation label image → `moCorROI`: the field set `TIFcatROIgui` writes, closed **traced** outlines (an angular sort folds on a concave mask, and that polygon is what FISSA turns into a neuropil annulus), row-major ID order, and the area/edge/dilate QC filters |
| `testConsensusROIsets` | voting across per-tif ROI sets: `minVotes` = 1 is the union and = nSets the strict inner join, per-cell detection counts, no merging of adjacent cells, and errors for an impossible threshold or mismatched frame sizes |
| `testHeadlessConfig` | every `processAnimal2Pheadless` selector form (indices, names, regexp, logical mask, `'auto'` bytes heuristic), stage ranges, the calibrated ROI defaults, and — the point of the test — that a selector matching **nothing** errors instead of silently relabelling the whole session |

**integration/**

| Test | Covers |
|---|---|
| `testCohortPhase01` | `stimGroupSpec` / `validateStimGroup` against current group files and the pre-fix `bkmat` backups; degenerate-*n* contract |
| `testAggregateStimGroup` | aggregator round-trip — regenerating a group reproduces the current table exactly, adding only `groupInfo` |
| `testPlotBPNgroup` | BPN group plotter; RLF small-*n* unification; reproduces `RLF_GroupD_nConsec2.mat` |
| `testPlotCGCgroup` | CGC group plotter, including the single-animal zero-significant-cell group that crashed the old path |
| `testCGCtwoStage` | CGC `_raw` split: `processCGC` reads `_raw`, writes processed, never mutates its input |
| `testProcessAnimalStimFamilies` | `processAnimal2P` §11 auto-run of the per-stim scripts |
| `testPhase6Cleanup` | group manifests, retired `combinetable`, `compileCohortData` fix |
| `testProcessBPN2P_e2e` | BPN compute path through `pkFcalc`; per-row structural alignment |
| `testProcessBPN2P_robust` | `baselineSec` sweep, guard error, single-onset table, sub-frame pulse |
| `testProcessBPN2P_saveAndPlots` | `processBPN2P` under `-batch`; saved bundle contents; `_raw` left untouched |
| `testFRAsignificance` | FRA significance rework — `pkFcalc` back-compat, logical mask (no peak²), trial-count invariance, onset column, `ceil` window, linear-index ordering, real TO0003 values |
| `testFRAgroup` | FRA-by-group path — `stimGroupSpec` FRA entry, `FRAmap2table` round-trip, `tableFRAmetrics` BF/bandwidth/threshold, sham noise floor |
| `testRemapROItoAcq_centeredCrop` | remapping 256×256 ROIs onto a centered-crop 256×128 acquisition, against a bit-exact fabricated ground truth |
| `testAA0072_pipeline` | the 256×128 ROI-reuse path end to end, as far as the interactive GUIs and the Python FISSA step allow |

`testFRAgroup` and `testFRAsignificance` guard their real-data sections
internally rather than skipping the whole file, so without the data drive they
pass on their synthetic checks alone. Read their output, not just the status,
when the drive is absent.

## Known gap pinned by a test

`testProcessAnimalStimFamilies` asserts that FRA is reported as **skipped**.
FRA is registered in `stimGroupSpec` but has no `_raw` stage (`suffixRaw` is
`''`), so `processAnimalStimFamilies`'s `isfile([animal suffixRaw])` check can
never hit and it can never drive `processFRA` — the per-animal FRA table has to
come from running `processFRA` directly. The assertion pins current behaviour,
it does not endorse it. If `processAnimalStimFamilies` grows a no-`_raw` branch,
that check is the one to revisit.
