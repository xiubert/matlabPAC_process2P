# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:54:08 2019
@author: PAC | pac94@pitt.edu

Neuropil correction via FISSA for the 2P pipeline.

FISSA's output is a rectangular cells x trials structure, so a single FISSA
run requires the SAME ROI count in every trial. Sessions that mix 256x256
(full ROI set) and 256x128 (reduced, remapped set) therefore cannot be one
run. This driver groups the ROI files by ROI count and runs FISSA once per
group, using EXPLICIT per-image ROI lists (no reliance on folder sort order).

  - Single count-group (the usual case, incl. no 256x128): behaves like the
    original driver -> writes NoRMCorred/FISSAoutput/matlab.mat, no manifest.
    Fully backward compatible.
  - Multiple count-groups: each group runs separately ->
    NoRMCorred/FISSAoutput/g<k>/matlab.mat, plus a manifest groups.json that
    maps each group output to its ordered tif basenames (trial order) and ROI
    count. processAnimal2P section 9 reads the manifest and merges per tif.

Usage:
    python FISSAviaMatlab_prePostTreatment.py [animalDataPath]
(animalDataPath defaults to the hardcoded path below if no arg is given.)
"""

import os
import sys
import glob
import json
from scipy.io import loadmat
import numpy as np
import fissa

# ---- animal data path: CLI arg overrides the default ------------------------
DEFAULT_PATH = '/media/user/sutter2P_data/Data/AA0343'
animalDataPath = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATH

tiff_folder = os.path.join(animalDataPath, 'NoRMCorred')
if not os.path.isdir(tiff_folder):
    raise NotADirectoryError('NoRMCorred folder not found: %s' % tiff_folder)
output_folder = os.path.join(tiff_folder, 'FISSAoutput')
os.makedirs(output_folder, exist_ok=True)


def parse_roi_file(path):
    """Return (rois, nTif, names, count) for one *_moCorrROI_*.mat file.

    rois  : list of ROIs, each ROI = [x_array, y_array] (FISSA polygon format,
            matching the legacy driver's extraction of ROIcurveOrderedXY).
    names : list of moCorr tif basenames this ROI set applies to (or None).
    """
    m = loadmat(path)
    field = np.asarray(m['moCorROI']['ROIcurveOrderedXY']).ravel()  # one entry/ROI
    rois = []
    for i in range(field.size):
        c = np.asarray(field[i])          # 2 x nPts : row0 = x, row1 = y
        rois.append([c[0], c[1]])
    nTif = int(np.asarray(m['nTifs']).ravel()[0])
    names = None
    if 'moCorTifNames' in m:
        names = [str(np.atleast_1d(x).ravel()[0])
                 for x in np.asarray(m['moCorTifNames']).ravel()]
    return rois, nTif, names, len(rois)


# ---- build per-image (path, roi-set) pairs, keyed by ROI count -------------
roi_files = sorted(glob.glob(os.path.join(animalDataPath, '*_moCorrROI_*.mat')))
if not roi_files:
    raise FileNotFoundError('No *_moCorrROI_*.mat in %s' % animalDataPath)

groups = {}   # count -> list of (image_path, roi_set, basename)
for rf in roi_files:
    rois, nTif, names, count = parse_roi_file(rf)
    if names is None:
        raise KeyError(
            '%s lacks moCorTifNames; re-save ROI files with the updated '
            'pipeline (processAnimal2P section 5/5b) so per-group image lists '
            'can be built.' % os.path.basename(rf))
    if len(names) != nTif:
        print('WARNING: %s nTifs=%d but %d names' % (os.path.basename(rf), nTif, len(names)))
    for nm in names:
        groups.setdefault(count, []).append((os.path.join(tiff_folder, nm), rois, nm))

# sort each group's images by basename so trial order is deterministic and
# matches the sorted tif order section 9 assumes for the legacy single-group path
for count in groups:
    groups[count].sort(key=lambda t: t[2])

print('ROI-count groups: %s' % {c: len(v) for c, v in groups.items()})


def run_group(image_paths, roi_list, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    exp = fissa.Experiment(image_paths, roi_list, out_dir)
    exp.separate(redo_prep=True)
    cwd = os.getcwd()
    os.chdir(out_dir)
    try:
        exp.save_to_matlab()           # writes matlab.mat in out_dir
    finally:
        os.chdir(cwd)


if len(groups) == 1:
    # ---- backward-compatible single-group path: legacy matlab.mat, no manifest
    count = next(iter(groups))
    pairs = groups[count]
    images = [p for p, _, _ in pairs]
    rois_list = [r for _, r, _ in pairs]
    run_group(images, rois_list, output_folder)
    print('DONE (single group, %d ROIs) -> FISSAoutput/matlab.mat' % count)
else:
    # ---- multi-group: one FISSA run per ROI count + manifest ---------------
    manifest = []
    for k, count in enumerate(sorted(groups, reverse=True)):  # largest count first
        pairs = groups[count]
        images = [p for p, _, _ in pairs]
        rois_list = [r for _, r, _ in pairs]
        gdir = os.path.join(output_folder, 'g%d' % k)
        run_group(images, rois_list, gdir)
        manifest.append({
            'matfile': os.path.join('g%d' % k, 'matlab.mat'),
            'nROI': int(count),
            'tifNames': [os.path.basename(p) for p, _, _ in pairs],  # trial order
        })
        print('  group %d: %d ROIs, %d tifs -> %s' % (k, count, len(images), manifest[-1]['matfile']))
    with open(os.path.join(output_folder, 'groups.json'), 'w') as fh:
        json.dump(manifest, fh, indent=2)
    print('DONE (%d groups) -> FISSAoutput/groups.json + g*/matlab.mat' % len(manifest))

# in 'result', for a given cell and trial there is a n x numTraceFrames double;
# row 1 is ROI trace, rows 2->n are neuropil-region traces around the ROI.
