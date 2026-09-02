#!/usr/bin/env python
"""
Suite2p ROI DETECTION only, for the matlabPAC_process2P headless path.

The MATLAB pipeline already does its own motion correction (NoRMCorre) and its
own neuropil correction (FISSA), so this runs suite2p with registration and
deconvolution OFF and uses it purely as a cell detector: registered tifs in,
a uint16 label image out, which labelImg2moCorROI turns into the usual
moCorROI struct.

Why suite2p alongside Cellpose: 'sparsery' (the default) and 'sourcery' decide
what a cell is from the movie's TEMPORAL structure -- sparse, localised,
co-modulating pixel groups -- rather than from morphology in a static
projection. That is the criterion the downstream dF/F analysis actually cares
about, and it does not depend on how noisy a mean image happens to be.

Usage
    python suite2pDetect.py --tif-list files.txt --out-dir OUT \
        --fs 5.0 --tau 1.0 --diameter 14 --algorithm sparsery

Writes into OUT:
    <prefix>_labels.tif   uint16 label image (0 = background)
    <prefix>_s2p.json     per-ROI summary + the settings actually used
    suite2p/              suite2p's own output tree (stat.npy, ops.npy, ...)

Overlapping suite2p masks are flattened by assigning each contested pixel to
whichever ROI weights it highest (lam), because the rest of the pipeline
extracts F with disjoint logical masks.
"""

import argparse
import json
import os
import sys

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--tif-list', help='text file, one tif path per line, in order')
    src.add_argument('--tif-dir', help='folder of tifs to use (sorted by name)')
    p.add_argument('--out-dir', required=True)
    p.add_argument('--prefix', default='s2p')

    p.add_argument('--fs', type=float, required=True, help='frame rate (Hz)')
    p.add_argument('--tau', type=float, default=1.0,
                   help='indicator decay (s): GCaMP6f ~0.7, 6m ~1.0, 6s ~1.25')
    p.add_argument('--diameter', type=float, default=0,
                   help='expected cell diameter in px; 0 = suite2p default')
    p.add_argument('--algorithm', default='sparsery',
                   choices=['sparsery', 'sourcery', 'cellpose'])
    p.add_argument('--threshold-scaling', type=float, default=1.0,
                   help='lower finds more ROIs, higher fewer')
    p.add_argument('--highpass-neuropil', type=float, default=None,
                   help='sparsery: spatial high-pass px, ~3x cell diameter')
    p.add_argument('--max-overlap', type=float, default=0.75)
    p.add_argument('--max-rois', type=int, default=5000)
    p.add_argument('--spatial-scale', type=int, default=0)
    p.add_argument('--denoise', action='store_true')
    p.add_argument('--use-classifier', action='store_true',
                   help="apply suite2p's built-in cell classifier and keep only iscell")
    p.add_argument('--keep-bin', action='store_true',
                   help='keep suite2p data.bin (same size as the input tifs)')
    p.add_argument('--device', default=None, help="'cuda' or 'cpu'; default auto")
    return p.parse_args()


def resolve_tifs(args):
    if args.tif_list:
        with open(args.tif_list) as fh:
            tifs = [ln.strip() for ln in fh if ln.strip()]
    else:
        tifs = sorted(os.path.join(args.tif_dir, f)
                      for f in os.listdir(args.tif_dir) if f.endswith('.tif'))
    missing = [t for t in tifs if not os.path.isfile(t)]
    if missing:
        raise FileNotFoundError('tif(s) not found: %s' % ', '.join(missing[:5]))
    if not tifs:
        raise ValueError('no tifs to run on')
    # suite2p takes one data_path plus a list of basenames
    dirs = {os.path.dirname(os.path.abspath(t)) for t in tifs}
    if len(dirs) != 1:
        raise ValueError('all tifs must live in one folder; got %d folders' % len(dirs))
    return dirs.pop(), [os.path.basename(t) for t in tifs]


def build_settings(args):
    import suite2p
    s = suite2p.default_settings()

    s['fs'] = args.fs
    s['tau'] = args.tau
    if args.diameter and args.diameter > 0:
        s['diameter'] = [float(args.diameter), float(args.diameter)]
    if args.device:
        s['torch_device'] = args.device

    # detection only: NoRMCorre already registered these, and spikes are not
    # part of this pipeline
    s['run']['do_registration'] = 0
    s['run']['do_regmetrics'] = False
    s['run']['do_detection'] = True
    s['run']['do_deconvolution'] = False

    s['io']['delete_bin'] = not args.keep_bin
    s['io']['save_mat'] = False
    s['io']['combined'] = False

    d = s['detection']
    d['algorithm'] = args.algorithm
    d['threshold_scaling'] = args.threshold_scaling
    d['max_overlap'] = args.max_overlap
    d['denoise'] = bool(args.denoise)
    d['sparsery_settings']['max_ROIs'] = args.max_rois
    d['sparsery_settings']['spatial_scale'] = args.spatial_scale
    if args.highpass_neuropil is not None:
        d['sparsery_settings']['highpass_neuropil'] = args.highpass_neuropil
    elif args.diameter and args.diameter > 0:
        # the docs' rule of thumb: ~3x the cell diameter in pixels
        d['sparsery_settings']['highpass_neuropil'] = float(round(3 * args.diameter))

    s['classification']['use_builtin_classifier'] = bool(args.use_classifier)
    return s


def stat_to_label(stat, iscell, shape):
    """Flatten suite2p's weighted, possibly overlapping masks to a label image.

    Each pixel goes to the ROI that weights it highest; ties keep the first.
    """
    Ly, Lx = shape
    label = np.zeros((Ly, Lx), dtype=np.uint16)
    best = np.zeros((Ly, Lx), dtype=np.float64)

    keep = [i for i in range(len(stat)) if iscell[i]]
    for n, i in enumerate(keep, start=1):
        ypix = np.asarray(stat[i]['ypix']).ravel()
        xpix = np.asarray(stat[i]['xpix']).ravel()
        lam = np.asarray(stat[i]['lam']).ravel().astype(np.float64)
        if ypix.size == 0:
            continue
        # normalising per ROI keeps a big dim cell from losing every contested
        # pixel to a small bright one
        lam = lam / max(lam.max(), 1e-12)
        take = lam > best[ypix, xpix]
        label[ypix[take], xpix[take]] = n
        best[ypix[take], xpix[take]] = lam[take]
    return label, keep


def main():
    args = parse_args()
    data_path, tif_names = resolve_tifs(args)
    os.makedirs(args.out_dir, exist_ok=True)

    import suite2p
    import tifffile

    settings = build_settings(args)
    db = {
        'data_path': [data_path],
        'tiff_list': tif_names,
        'save_path0': args.out_dir,
        'look_one_level_down': False,
    }

    print('suite2p detection: %d tifs from %s' % (len(tif_names), data_path), flush=True)
    print('  algorithm=%s fs=%g tau=%g diameter=%s threshold_scaling=%g'
          % (args.algorithm, args.fs, args.tau, settings['diameter'],
             args.threshold_scaling), flush=True)

    suite2p.run_s2p(db=db, settings=settings)

    plane = os.path.join(args.out_dir, 'suite2p', 'plane0')
    stat = np.load(os.path.join(plane, 'stat.npy'), allow_pickle=True)
    ops = np.load(os.path.join(plane, 'ops.npy'), allow_pickle=True).item()

    iscell_path = os.path.join(plane, 'iscell.npy')
    if args.use_classifier and os.path.isfile(iscell_path):
        iscell = np.load(iscell_path)[:, 0].astype(bool)
    else:
        iscell = np.ones(len(stat), dtype=bool)

    shape = (int(ops['Ly']), int(ops['Lx']))
    label, kept = stat_to_label(stat, iscell, shape)

    label_path = os.path.join(args.out_dir, '%s_labels.tif' % args.prefix)
    tifffile.imwrite(label_path, label)

    rois = []
    for n, i in enumerate(kept, start=1):
        st = stat[i]
        rois.append({
            'label': n,
            'npix': int(st.get('npix', len(st['ypix']))),
            'med': [float(v) for v in np.asarray(st.get('med', [np.nan, np.nan])).ravel()[:2]],
            'compact': float(st.get('compact', np.nan)),
            'aspect_ratio': float(st.get('aspect_ratio', np.nan)),
            'radius': float(st.get('radius', np.nan)),
            'npix_in_label': int((label == n).sum()),
        })

    summary = {
        'labelTif': label_path,
        'nDetected': int(len(stat)),
        'nKept': int(len(kept)),
        'nInLabelImage': int(label.max()),
        'shape': list(shape),
        'algorithm': args.algorithm,
        'suite2pVersion': settings.get('version', ''),
        'settings': {
            'fs': args.fs, 'tau': args.tau, 'diameter': settings['diameter'],
            'threshold_scaling': args.threshold_scaling,
            'highpass_neuropil': settings['detection']['sparsery_settings']['highpass_neuropil'],
            'max_overlap': args.max_overlap,
            'spatial_scale': args.spatial_scale,
            'denoise': bool(args.denoise),
            'use_classifier': bool(args.use_classifier),
        },
        'tifs': tif_names,
        'rois': rois,
    }
    json_path = os.path.join(args.out_dir, '%s_s2p.json' % args.prefix)
    with open(json_path, 'w') as fh:
        json.dump(summary, fh, indent=2)

    print('DONE: %d ROIs detected, %d written to %s'
          % (len(stat), label.max(), label_path), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
