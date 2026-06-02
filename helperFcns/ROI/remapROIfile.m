function moCorROI = remapROIfile(srcROIpath, srcTif, tgtTif, opts)
% remapROIfile  Pipeline driver: remap a saved moCorROI file onto a 256x128 acq.
%
%   moCorROI = remapROIfile(srcROIpath, srcTif, tgtTif)
%   moCorROI = remapROIfile(..., 'outPath',P, 'nTifs',N, ...
%                           'tifIDXinAllTifList',IDX, 'moCorSeqN',S, ...
%                           'neuropilMarginPx',M)
%
%   Loads the source ROI set (drawn on a 256x256 acquisition), reads the scan
%   geometry from a representative source tif and target tif, and remaps the
%   ROIs into the target (256x128 / 10 Hz) frame via remapROItoAcq -- keeping
%   only ROIs fully contained in the centered crop. Intended to replace the
%   TIFcatROIgui "draw ROIs" step (processAnimal2P section 4) for a condition
%   whose cells were already drawn on the full 256x256 field.
%
%   Inputs
%     srcROIpath  path to <animal>_moCorrROI_<srcCond>.mat (contains moCorROI)
%     srcTif      a representative tif from the SOURCE (256x256) condition
%     tgtTif      a representative tif from the TARGET (256x128) condition
%
%   Name-value (required together only when saving)
%     outPath              if given, write a pipeline-format bundle here
%     nTifs                number of tifs in the target condition (FISSA reads
%                          this); REQUIRED when outPath is set
%     tifIDXinAllTifList   logical index of target tifs within the full tif
%                          list; REQUIRED when outPath is set
%     moCorSeqN            optional condition index, stored for provenance
%     neuropilMarginPx     forwarded to remapROItoAcq (default 0)
%
%   Output
%     moCorROI    remapped struct array (target coordinates), also saved to
%                 outPath (as separate variables moCorROI / nTifs /
%                 tifIDXinAllTifList / moCorSeqN) when requested.
%
%   Errors (via remapROItoAcq) if the source/target geometry is not the
%   expected centered, zoom-matched, unrotated 1:1 crop.
%
%   Example (inside processAnimal2P, as "section 4b" for a reused condition):
%     srcCond = moCorN{srcIdx};  tgtCond = moCorN{tgtIdx};
%     nT  = numel(tifList.(tgtCond));
%     idx = ismember({tifFiles.name}', {tifList.(tgtCond).name}');
%     remapROIfile( ...
%        fullfile(dataPath,[animal '_moCorrROI_' srcCond '.mat']), ...
%        fullfile(tifList.(srcCond)(1).folder, tifList.(srcCond)(1).name), ...
%        fullfile(tifList.(tgtCond)(1).folder, tifList.(tgtCond)(1).name), ...
%        'outPath', fullfile(dataPath,[animal '_moCorrROI_' tgtCond '.mat']), ...
%        'nTifs', nT, 'tifIDXinAllTifList', idx, 'moCorSeqN', tgtIdx);
%
%   See also remapROItoAcq, intersectROIfiles, readSCIMtif, TIFcatROIgui

arguments
    srcROIpath (1,:) char
    srcTif     (1,:) char
    tgtTif     (1,:) char
    opts.outPath            (1,:) char = ''
    opts.nTifs              double     = []
    opts.tifIDXinAllTifList            = []
    opts.moCorSeqN          double     = []
    opts.neuropilMarginPx   (1,1) double {mustBeNonnegative} = 0
end

S = load(srcROIpath, 'moCorROI');
if ~isfield(S, 'moCorROI')
    error('remapROIfile:noROI', '%s does not contain moCorROI.', srcROIpath);
end

[~, srcHeader] = readSCIMtif(srcTif, 'metaOnly');
[~, tgtHeader] = readSCIMtif(tgtTif, 'metaOnly');

moCorROI = remapROItoAcq(S.moCorROI, srcHeader, tgtHeader, opts.neuropilMarginPx);

fprintf('remapROIfile: kept %d of %d ROIs (target %dx%d)\n', ...
    numel(moCorROI), numel(S.moCorROI), tgtHeader.imHeight, tgtHeader.imWidth);

% --- optional save in pipeline format -------------------------------------
if ~isempty(opts.outPath)
    if isempty(opts.nTifs) || isempty(opts.tifIDXinAllTifList)
        error('remapROIfile:missingMeta', ...
            ['Saving requires nTifs and tifIDXinAllTifList (FISSA and ' ...
             'intersectROIfiles read them). Pass both name-value pairs.']);
    end
    out.moCorROI           = moCorROI;
    out.nTifs              = opts.nTifs;
    out.tifIDXinAllTifList = opts.tifIDXinAllTifList;
    if ~isempty(opts.moCorSeqN)
        out.moCorSeqN = opts.moCorSeqN;
    end
    save(opts.outPath, '-struct', 'out');
    fprintf('remapROIfile: saved %s\n', opts.outPath);
end
end
