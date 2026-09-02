function cfg = headlessConfig(dataPath,varargin)
% headlessConfig  Build and validate the config for processAnimal2Pheadless.
%
%   cfg = headlessConfig(dataPath)
%   cfg = headlessConfig(dataPath,'treatmentName','ZX1','preTifs',1:8,...)
%   cfg = headlessConfig(dataPath,existingCfgStruct)
%
%   Every dialog processAnimal2P puts on screen becomes a field here, so the
%   same pipeline can run unattended. Defaults reproduce what a user would
%   normally click, and each is stated below rather than buried in the script.
%
%   Inputs
%     dataPath - animal data folder.
%
%   Name-value -- identity and tif inventory (§1)
%     'animal'        animal ID. Default: [A-Z]{2}\d{4} from dataPath.
%     'tifPattern'    which tifs are the session. Default '<animal>*.tif'.
%     'treatmentName' e.g. 'ZX1'. Default '' = every tif is 'none'.
%     'preTifs'       which tifs are PRE-treatment; the rest are post. Only
%                     used when treatmentName is set. Accepts numeric indices
%                     into the sorted tif list, a cellstr/string of file
%                     names, or a char regular expression matched on names.
%                     Default [] = every tif is post, which is almost never
%                     what you want, so it warns.
%     'mapTifs'       which tifs are FRA/BF maps. Same selector forms, plus
%                     'auto' (default): bytes > mapBytesThreshold, the same
%                     heuristic the interactive dialog prints as a hint.
%                     Pass [] for none.
%     'mapBytesThreshold'  Default 11e6.
%
%   Name-value -- condition split (§2)
%     'condFilters'   treatment substrings that define motion-correction
%                     groups, e.g. {'preZX1','postZX1'}. Default: derived
%                     from the treatments actually present. Geometry
%                     sub-splitting happens automatically on top of this.
%     'altZoomPolicy' what to do with tifs at a non-standard zoom:
%                     'drop' (default, matching the interactive default),
%                     'keep' (own geometry group) or 'error'.
%
%   Name-value -- motion correction (§3)
%     'funcChan'      functional channel to keep from multi-channel tifs
%                     (1 = red, 2 = green). Default 2.
%     'normcorre'     struct of NoRMCorreSetParms overrides. Default matches
%                     the interactive script: grid_size [32 32], mot_uf 4,
%                     bin_width 200, max_shift 15, max_dev 3, us_fac 50,
%                     init_batch 200.
%
%   Name-value -- ROI detection (§4-5), forwarded to cellposeROIset
%     'roi'           struct of cellposeROIset options. Defaults:
%                     mode 'consensus', minVotes 2, dilatePx 2,
%                     edgeMarginPx 3, cellposeArgs {diameter 15,
%                     cellprobThreshold -1}. These are the values calibrated
%                     against hand-drawn ROIs on TO0003.
%
%   Name-value -- FISSA (§8-9)
%     'fissaCmd'      shell command template for the Python step; '%s' is
%                     replaced by dataPath. Default runs the repo's
%                     FISSAviaMatlab_prePostTreatment.py in conda env
%                     'env_fissa'. Set '' to skip the call and assume the
%                     output already exists.
%     'fissaScaleFactor'  neuropil subtraction scale. Default 0.8.
%
%   Name-value -- per-stim analysis (§10-11)
%     'runStimFamilies'  Default true.
%     'runFRAmap'        run processFRA when the animal has map tifs.
%                        Default true.
%     'stimScriptVars'   struct forwarded to processAnimalStimFamilies.
%
%   Name-value -- control
%     'stages'      which stages to run, as a numeric vector or a [from to]
%                   range. Default 1:11. Stage numbers match
%                   processAnimal2P's section numbers.
%     'overwrite'   redo a stage whose artefact already exists. Default
%                   false, so a re-run resumes rather than repeating motion
%                   correction.
%     'verbose'     Default true.
%
%   Output
%     cfg  validated struct, with .selectors resolved against the tifs on
%          disk (.preIDX, .mapIDX) so the caller can see what was chosen
%          before anything is written.
%
%   See also processAnimal2Pheadless, cellposeROIset, processAnimal2P

if ~isfolder(dataPath)
    error('headlessConfig:noDataPath','Not a folder: %s',dataPath);
end

%allow headlessConfig(dataPath, cfgStruct) as well as name-value pairs
if isscalar(varargin) && isstruct(varargin{1})
    nv = struct2nv(varargin{1});
else
    nv = varargin;
end

animalGuess = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'animal',animalGuess,@(x)ischar(x)||isstring(x));
addParameter(p,'tifPattern','',@(x)ischar(x)||isstring(x));
addParameter(p,'treatmentName','',@(x)ischar(x)||isstring(x));
addParameter(p,'preTifs',[]);
addParameter(p,'mapTifs','auto');
addParameter(p,'mapBytesThreshold',11e6,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'condFilters',{},@(x)isempty(x)||iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(p,'altZoomPolicy','drop',@(s)any(strcmpi(s,{'drop','keep','error'})));
addParameter(p,'funcChan',2,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'normcorre',struct(),@isstruct);
addParameter(p,'roi',struct(),@isstruct);
addParameter(p,'fissaCmd','default',@(x)ischar(x)||isstring(x));
addParameter(p,'fissaScaleFactor',0.8,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'runStimFamilies',true,@islogical);
addParameter(p,'runFRAmap',true,@islogical);
addParameter(p,'stimScriptVars',struct(),@isstruct);
addParameter(p,'stages',1:11,@isnumeric);
addParameter(p,'overwrite',false,@islogical);
addParameter(p,'verbose',true,@islogical);
parse(p,nv{:});
cfg = p.Results;

cfg.dataPath = char(dataPath);
cfg.animal   = char(cfg.animal);
if isempty(cfg.animal)
    error('headlessConfig:noAnimal',...
        'Could not derive an animal ID from %s. Pass ''animal'' explicitly.',dataPath);
end
if isempty(char(cfg.tifPattern))
    cfg.tifPattern = [cfg.animal '*.tif'];
else
    cfg.tifPattern = char(cfg.tifPattern);
end
cfg.treatmentName = char(cfg.treatmentName);

%--- stages ---------------------------------------------------------------
st = cfg.stages(:)';
if numel(st) == 2 && st(2) > st(1) + 1
    st = st(1):st(2);      % [from to] range
end
cfg.stages = unique(st);
if any(cfg.stages < 1 | cfg.stages > 11)
    error('headlessConfig:badStages',...
        'stages must lie in 1..11 (got %s).',mat2str(cfg.stages));
end

%--- NoRMCorre parameters -------------------------------------------------
ncDefault = struct('grid_size',[32 32],'mot_uf',4,'bin_width',200,...
    'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
cfg.normcorre = mergeStruct(ncDefault,cfg.normcorre);

%--- ROI detection defaults ----------------------------------------------
roiDefault = struct('mode','consensus','minVotes',2,'dilatePx',2,...
    'minAreaPx',20,'edgeMarginPx',3,'maxShiftPx',15,'saveMeanTif',true,...
    'cellposeArgs',{{'diameter',15,'cellprobThreshold',-1}});
cfg.roi = mergeStruct(roiDefault,cfg.roi);

%--- FISSA command --------------------------------------------------------
if strcmp(char(cfg.fissaCmd),'default')
    repoRoot = fileparts(fileparts(mfilename('fullpath')));   % helperFcns/..
    repoRoot = fileparts(repoRoot);
    pyScript = fullfile(repoRoot,'FISSAviaMatlab_prePostTreatment.py');
    cfg.fissaCmd = sprintf(...
        'bash -lc ''source ~/miniconda3/bin/activate env_fissa && python "%s" "%%s"''',...
        pyScript);
else
    cfg.fissaCmd = char(cfg.fissaCmd);
end

%--- resolve tif selectors against what is on disk -----------------------
tifFiles = dir(fullfile(cfg.dataPath,cfg.tifPattern));
tifFiles = tifFiles(~[tifFiles.isdir]);
cfg.nTifsFound = numel(tifFiles);
if isempty(tifFiles)
    warning('headlessConfig:noTifs',...
        'No tifs matching %s in %s (fine only if resuming past stage 1).',...
        cfg.tifPattern,cfg.dataPath);
    cfg.preIDX = false(0,1);
    cfg.mapIDX = false(0,1);
    cfg.tifNames = {};
    return
end
cfg.tifNames = {tifFiles.name}';

%pre/post
if isempty(cfg.treatmentName)
    cfg.preIDX = false(numel(tifFiles),1);
else
    cfg.preIDX = resolveSelector(cfg.preTifs,tifFiles,'preTifs');
    if ~any(cfg.preIDX)
        warning('headlessConfig:noPreTifs',...
            ['treatmentName ''%s'' is set but no tif matched preTifs, so every '...
             'tif would be labelled post%s. Check the selector.'],...
            cfg.treatmentName,cfg.treatmentName);
    end
end

%FRA map tifs
if ischar(cfg.mapTifs) && strcmpi(cfg.mapTifs,'auto')
    cfg.mapIDX = [tifFiles.bytes]' > cfg.mapBytesThreshold;
    cfg.mapSelector = sprintf('auto (bytes > %g)',cfg.mapBytesThreshold);
else
    cfg.mapIDX = resolveSelector(cfg.mapTifs,tifFiles,'mapTifs');
    cfg.mapSelector = 'explicit';
end

if cfg.verbose
    fprintf('headlessConfig: %s, %d tifs\n',cfg.animal,numel(tifFiles));
    if ~isempty(cfg.treatmentName)
        fprintf('  treatment %s: %d pre / %d post\n',cfg.treatmentName,...
            nnz(cfg.preIDX),nnz(~cfg.preIDX));
    end
    fprintf('  FRA map tifs: %d [%s]\n',nnz(cfg.mapIDX),cfg.mapSelector);
    fprintf('  stages %s | overwrite %d\n',mat2str(cfg.stages),cfg.overwrite);
end
end

% ------------------------------------------------------------------------
function idx = resolveSelector(sel,tifFiles,name)
%accept indices, names, or a regular expression; always return a logical
%index over tifFiles, and refuse anything that silently selects nothing it
%should have
n = numel(tifFiles);
idx = false(n,1);
if isempty(sel)
    return
end
names = {tifFiles.name}';

if isnumeric(sel) || islogical(sel)
    if islogical(sel)
        if numel(sel) ~= n
            error('headlessConfig:badSelector',...
                '%s logical index has %d elements but there are %d tifs.',...
                name,numel(sel),n);
        end
        idx = sel(:);
    else
        bad = sel(sel < 1 | sel > n | sel ~= round(sel));
        if ~isempty(bad)
            error('headlessConfig:badSelector',...
                '%s contains out-of-range indices: %s (1..%d valid).',...
                name,mat2str(bad),n);
        end
        idx(sel) = true;
    end
    return
end

if iscellstr(sel) || isstring(sel) %#ok<ISCLSTR>
    want = cellstr(sel);
    hit = ismember(names,want);
    missing = want(~ismember(want,names));
    if ~isempty(missing)
        error('headlessConfig:badSelector',...
            '%s names not found among the tifs: %s',name,strjoin(missing',', '));
    end
    idx = hit;
    return
end

if ischar(sel)
    hit = ~cellfun(@isempty,regexp(names,sel,'once'));
    if ~any(hit)
        error('headlessConfig:badSelector',...
            '%s regular expression ''%s'' matched none of the %d tifs.',...
            name,sel,n);
    end
    idx = hit;
    return
end

error('headlessConfig:badSelector',...
    '%s must be indices, a logical mask, file names, or a regexp (got %s).',...
    name,class(sel));
end
% ------------------------------------------------------------------------
function s = mergeStruct(s,over)
f = fieldnames(over);
for k = 1:numel(f)
    s.(f{k}) = over.(f{k});
end
end
% ------------------------------------------------------------------------
function nv = struct2nv(s)
f = fieldnames(s);
nv = cell(1,2*numel(f));
for k = 1:numel(f)
    nv{2*k-1} = f{k};
    nv{2*k}   = s.(f{k});
end
end
