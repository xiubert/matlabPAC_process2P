function P = animalPaths(dataPath,opts)
% animalPaths  Resolve where an animal's inputs and artifacts live.
%
%   P = animalPaths(dataPath)                     % legacy flat layout
%   P = animalPaths(dataPath,'run','cellpose_20260903')
%   P = animalPaths(dataPath,'artifactDir',D)     % explicit override
%
%   Single source of truth for the layout, so no function has to build a path
%   by hand. Everything the pipeline reads or writes falls into one of three
%   groups, and they have different lifetimes:
%
%     RAW        the acquisition: *.tif, *_Pulses.mat, *_PulseParams.mat.
%                Never written by the pipeline.
%     SHARED     depends only on the acquisition and the motion correction:
%                the tif inventory and condition-split legends, NoRMCorred/
%                and its NoRMCorreParams. Expensive to produce and identical
%                for every analysis of this animal, so it is NEVER duplicated
%                per run.
%     PER-RUN    depends on the ROI set and the analysis parameters:
%                _moCorrROI_*, FISSAoutput/, _moCorr_Tifs_Params,
%                _tifFileList, _pulseLegend2P, _stimGroupIDX, _anmlROI_*,
%                _FRAmap, and the QC images.
%
%   With a run name, per-run artifacts go to <dataPath>/analysis/<run>/ and
%   two analyses of the same animal cannot overwrite each other. Without one,
%   they go to <dataPath> exactly as they always have, so existing folders and
%   scripts keep working untouched.
%
%   Inputs
%     dataPath  animal data folder (holds the raw tifs).
%
%   Name-value
%     'run'          run name; becomes the folder under analysis/. Must be a
%                    valid folder name. Default '' = legacy flat layout.
%     'artifactDir'  put per-run artifacts here instead of deriving the path.
%                    Wins over 'run'. For callers that already know the folder.
%     'create'       create the directories. Default false, so a read-only
%                    caller never leaves empty folders behind.
%
%   Output (struct)
%     .root        the animal folder (raw inputs, and where legends live)
%     .artifacts   per-run artifacts
%     .moCorrDir   NoRMCorred/ (shared)
%     .fissaDir    FISSA output for THIS run
%     .qcDir       QC images for this run
%     .run         the run name ('' when flat)
%     .isFlat      true when artifacts == root (the legacy layout)
%     .legend, .condSplit, .ncParams, .moCorrTifs, .tifFileList
%                  fully resolved paths to the named artifacts
%
%   Legacy folders stay readable: an animal processed before this layout has
%   no run, resolves flat, and finds its FISSA output under NoRMCorred/ as
%   before. A NAMED run always owns its own FISSAoutput -- it never falls back
%   to the shared one, because that is the collision this layout prevents.
%
%   See also processAnimal2Pheadless, migrateAnimalArtifacts, headlessConfig

arguments
    dataPath        (1,:) char
    opts.run        (1,:) char = ''
    opts.artifactDir(1,:) char = ''
    opts.create     (1,1) logical = false
end

if ~isfolder(dataPath)
    error('animalPaths:noDataPath','Not a folder: %s',dataPath);
end

animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');

P.root      = dataPath;
P.run       = opts.run;
P.moCorrDir = fullfile(dataPath,'NoRMCorred');

if ~isempty(opts.artifactDir)
    P.artifacts = opts.artifactDir;
    if isempty(P.run)
        [~,P.run] = fileparts(P.artifacts);
    end
elseif ~isempty(opts.run)
    validateRunName(opts.run);
    P.artifacts = fullfile(dataPath,'analysis',opts.run);
else
    P.artifacts = dataPath;
end
P.isFlat = strcmp(P.artifacts,P.root);

P.qcDir = fullfile(P.artifacts,'QC');
if P.isFlat
    P.qcDir = P.artifacts;      % nothing to separate from in a flat folder
end

% FISSA output belongs to the ROI set that produced it, so it is per-run.
% A NAMED RUN always owns its own FISSA output. Falling back to the shared
% NoRMCorred/FISSAoutput just because the per-run one does not exist yet would
% send every run's output to the same place -- which is the collision this
% layout exists to prevent. The legacy location is used only in flat mode,
% which is what a pre-layout animal resolves to anyway.
if P.isFlat
    P.fissaDir = fullfile(P.moCorrDir,'FISSAoutput');
else
    P.fissaDir = fullfile(P.artifacts,'FISSAoutput');
end

% shared: these describe the acquisition and the motion correction, and are
% the same for every analysis of this animal
P.legend    = fullfile(P.root,[animal '_tifFileLegend.mat']);
P.condSplit = fullfile(P.root,[animal '_tifCondSplitLegend.mat']);
P.ncParams  = fullfile(P.moCorrDir,[animal '_NoRMCorreParams.mat']);

% per-run
P.moCorrTifs  = fullfile(P.artifacts,[animal '_moCorr_Tifs_Params.mat']);
P.tifFileList = fullfile(P.artifacts,[animal '_tifFileList.mat']);

if opts.create
    for f = {P.artifacts,P.qcDir,P.fissaDir}
        if ~isfolder(f{1}); mkdir(f{1}); end
    end
end
end

% ------------------------------------------------------------------------
function validateRunName(r)
if ~isequal(r,matlab.lang.makeValidName(strrep(r,'-','_')))
    % allow hyphens and dots, which makeValidName would mangle, but nothing
    % that would climb out of the folder or break a glob
    if ~isempty(regexp(r,'[\\/:*?"<>|]|\.\.','once'))
        error('animalPaths:badRunName', ...
            ['Run name "%s" is not usable as a folder name. Use letters, ' ...
             'digits, hyphens, dots and underscores.'],r);
    end
end
end
