function info = aggregateStimGroup(manifest,varargin)
% aggregateStimGroup  Build one condition-group table from per-animal tables.
%
%   info = aggregateStimGroup(manifest)
%   info = aggregateStimGroup(manifestJsonPath)
%   info = aggregateStimGroup(...,'dryRun',true,'verbose',true)
%
%   Family-agnostic replacement for combinetable.m. Resolves each animal's
%   PROCESSED per-stim table, validates it, canonicalizes column order,
%   concatenates, stamps provenance, and saves the group file.
%
%   Aggregation happens AFTER per-animal dF/F and peak detection, because
%   baselines and significance are per-animal quantities. This function never
%   recomputes them.
%
%   Inputs:
%     manifest - struct, or path to a .json holding the same fields:
%                  .group      group suffix, e.g. 'D'
%                  .family     'BPN' | 'CGC'
%                  .animals    string array of animal IDs, IN THE ORDER they
%                              should be concatenated (row order follows it)
%                  .cohortRoot parent folder containing <animal>/ subfolders
%                  .outDir     where the group file is written
%                  .tableDir   (optional) subfolder under each animal folder;
%                              default '.'
%                  .files      (optional) explicit per-animal file paths,
%                              used instead of cohortRoot/animals resolution
%
%   Name/Value:
%     'dryRun'  - resolve, load and validate but write nothing. Default false.
%     'verbose' - print progress. Default true.
%     'refVars' - canonical column order. Default: the first animal's order.
%     'outFile' - override the output path.
%     'trimTraces' - trim every animal's dF/F columns to the shortest common
%                 frame count before concatenating. Default true. Recordings
%                 differ in length across animals, so without this a group can
%                 be ragged and unusable: any cohort plot that stacks per-cell
%                 traces fails, and validateStimGroup refuses the group with
%                 raggedTimeAxis. Only the dF/F axis columns (timeVar +
%                 traceVars) are trimmed; raw F columns sit on t_total and are
%                 left alone. What was trimmed is recorded in groupInfo.
%
%   Output (struct) - also saved into the group file as `groupInfo`:
%     .schemaVersion .group .family .created .createdBy
%     .manifest      the manifest as given
%     .animals .nAnimals .nCells .nRows .perAnimal
%     .timeVar .timeAxis .traceVars .cellAvgVar
%     .convention    analysis parameters in force (from stimGroupSpec)
%     .sourceFiles   table: animal, file, bytes, modified
%     .validation    validateStimGroup report for the concatenated table
%     .trim          whether a common-axis trim was applied, and to how many frames
%     .outFile
%
%   The group file contains the table under its family's usual variable name
%   (`anmlROIbyStim`) plus `groupInfo`, so existing scripts that `load` a
%   group file and reference `anmlROIbyStim` keep working unchanged. A group
%   file without `groupInfo` is legacy, not invalid.
%
%   `convention` is the field that would have made the CGC Group C/D
%   divergence visible immediately: two script versions produced incompatible
%   dFF_PT columns under identical filenames, and nothing recorded which
%   recipe had been used.
%
%   A manifest .json is written next to the group file so membership is
%   reviewable and diffable; the same content is embedded in groupInfo so the
%   group file is self-describing if separated from it.
%
%   See also loadStimGroup, validateStimGroup, stimGroupSpec.

p = inputParser;
addRequired(p,'manifest',@(x) isstruct(x) || ischar(x) || isstring(x));
addParameter(p,'dryRun',false,@islogical);
addParameter(p,'verbose',true,@islogical);
addParameter(p,'refVars',{},@(x) isempty(x) || iscellstr(x) || isstring(x)); %#ok<ISCLSTR>
addParameter(p,'outFile','',@(x) ischar(x)||isstring(x));
addParameter(p,'trimTraces',true,@islogical);
parse(p,manifest,varargin{:});
dryRun  = p.Results.dryRun;
verbose = p.Results.verbose;
refVars = p.Results.refVars;
outFile = char(p.Results.outFile);
trimTraces = p.Results.trimTraces;
if ~isempty(refVars); refVars = cellstr(refVars); end

%% ---- manifest ----
if ~isstruct(manifest)
    mPath = char(manifest);
    if ~isfile(mPath)
        error('aggregateStimGroup:manifestNotFound','Manifest not found: %s', mPath);
    end
    manifest = jsondecode(fileread(mPath));
end
manifest = normalizeManifest(manifest);
spec = stimGroupSpec(manifest.family);

%% ---- resolve per-animal files ----
[files, animals] = resolveFiles(manifest, spec);
if verbose
    fprintf('aggregateStimGroup: %s Group%s, %d animal(s)\n', ...
        manifest.family, manifest.group, numel(animals));
end

%% ---- load + validate each animal ----
per = cell(numel(animals),1);
srcAnimal = strings(0,1); srcFile = strings(0,1);
srcBytes = zeros(0,1); srcMod = strings(0,1);

for a = 1:numel(animals)
    f = files{a};
    if ~isfile(f)
        % Distinguish "never inventoried" from "inventoried but not yet
        % processed". The second is the common mistake: processAnimal2P ends
        % at stimParam2ROI, which writes only the _raw table, and aggregation
        % reads the PROCESSED one because that is what carries dF/F and peak
        % responses. The _raw convention exists precisely so this is an error
        % rather than a group silently missing its analysis columns.
        rawFile = fullfile(fileparts(f), [char(animals(a)) spec.suffixRaw]);
        if isfile(rawFile)
            error('aggregateStimGroup:notProcessed', ...
                ['%s has been through stimParam2ROI but not through %s.\n' ...
                 '    found:   %s\n' ...
                 '    missing: %s\n' ...
                 'Run %s on %s before aggregating -- the processed table is ' ...
                 'what carries dF/F and peak responses.'], ...
                animals(a), spec.processScript, rawFile, f, ...
                spec.processScript, animals(a));
        end
        error('aggregateStimGroup:fileNotFound', ...
            ['Processed table for %s not found: %s\n' ...
             'Expected <cohortRoot>/<animal>/%s. Has processAnimal2P run for ' ...
             'this animal?'], animals(a), f, ['<animal>' spec.suffixProcessed]);
    end
    S = load(f, spec.varname);
    if ~isfield(S, spec.varname)
        error('aggregateStimGroup:missingVar', ...
            '%s has no variable "%s" (has stimParam2ROI/process%s been run?)', ...
            f, spec.varname, manifest.family);
    end
    Ta = S.(spec.varname);

    rep = validateStimGroup(Ta, manifest.family, 'verbose',false);
    if ~rep.ok
        fprintf('  %s: VALIDATION FAILED\n', animals(a));
        disp(rep.problems(rep.problems.severity=="error",:));
        error('aggregateStimGroup:animalInvalid', ...
            '%s did not validate; fix the per-animal table before aggregating.', animals(a));
    end

    if isempty(refVars); refVars = Ta.Properties.VariableNames; end
    extra = setxor(Ta.Properties.VariableNames, refVars);
    if ~isempty(extra)
        error('aggregateStimGroup:schemaMismatch', ...
            ['%s has a different variable set than the reference ' ...
             '(differing: %s). Animals must share a schema.'], ...
            animals(a), strjoin(extra,', '));
    end
    Ta = Ta(:,refVars);                       % canonical column order

    d = dir(f);
    srcAnimal(end+1,1) = animals(a); %#ok<AGROW>
    srcFile(end+1,1)   = string(f);  %#ok<AGROW>
    srcBytes(end+1,1)  = d.bytes;    %#ok<AGROW>
    srcMod(end+1,1)    = string(datetime(d.datenum,'ConvertFrom','datenum', ...
        'Format','uuuu-MM-dd HH:mm:ss')); %#ok<AGROW>

    per{a} = Ta;
    if verbose
        Na = groupN(Ta);
        fprintf('  %-8s %4d rows, %3d cells\n', animals(a), height(Ta), Na.nCells);
    end
end

%% ---- trim to a common dF/F axis ----
% Recording length differs across animals, so the dF/F columns can be ragged
% even when every animal is individually valid. Trim to the shortest before
% concatenating -- the same repair that was applied by hand to the historical
% group files (etc/fixBPNgroupTraceLength.m). Raw F columns are on the
% t_total axis and are not touched.
trimInfo = struct('applied',false,'nFrames',NaN,'perAnimalBefore',[]);
lens = cellfun(@(T) unique(cellfun(@numel,T.(spec.timeVar))), per, 'uni',0);
allLens = unique([lens{:}]);
if trimTraces && numel(allLens) > 1
    nKeep = min(allLens);
    % every animal must share the first nKeep samples, or trimming would
    % silently align traces that are not on the same clock
    ref = per{1}.(spec.timeVar){1}(1:nKeep);
    for a = 1:numel(per)
        tt = cell2mat(cellfun(@(c) reshape(c(1:nKeep),1,[]), per{a}.(spec.timeVar), 'uni',0));
        if max(abs(tt - ref),[],'all') > 1e-9
            error('aggregateStimGroup:incompatibleTimeAxis', ...
                ['%s does not share the first %d samples of the common dF/F ' ...
                 'axis, so trimming would misalign it.'], animals(a), nKeep);
        end
    end
    for a = 1:numel(per)
        per{a}.(spec.timeVar) = cellfun(@(c) c(1:nKeep), per{a}.(spec.timeVar), 'uni',0);
        for v = 1:numel(spec.traceVars)
            vn = spec.traceVars{v};
            if ismember(vn, per{a}.Properties.VariableNames)
                per{a}.(vn) = cellfun(@(c) c(:,1:nKeep), per{a}.(vn), 'uni',0);
            end
        end
    end
    trimInfo = struct('applied',true,'nFrames',nKeep,'perAnimalBefore',allLens);
    if verbose
        fprintf('  trimmed dF/F columns %s -> %d frames (common axis)\n', ...
            mat2str(allLens), nKeep);
    end
end

%% ---- concatenate + validate the group ----
T = vertcat(per{:});
rep = validateStimGroup(T, manifest.family, 'refVars',refVars, 'verbose',false);
if ~rep.ok
    disp(rep.problems(rep.problems.severity=="error",:));
    error('aggregateStimGroup:groupInvalid', ...
        'Concatenated group did not validate (see problems above).');
end

N = groupN(T);

%% ---- provenance stamp ----
info.schemaVersion = 1;
info.group   = manifest.group;
info.family  = manifest.family;
info.created = string(datetime('now','Format','uuuu-MM-dd HH:mm:ss'));
info.createdBy = struct('script','aggregateStimGroup', ...
    'matlab',string(version('-release')), 'gitSHA',repoSHA());
info.manifest   = manifest;
info.animals    = N.animals;
info.nAnimals   = N.nAnimals;
info.nCells     = N.nCells;
info.nRows      = N.nRows;
info.perAnimal  = N.perAnimal;
info.timeVar    = spec.timeVar;
info.timeAxis   = T.(spec.timeVar){1};
info.traceVars  = spec.traceVars;
info.cellAvgVar = spec.cellAvgVar;
info.convention = spec.convention;
info.sourceFiles = table(srcAnimal,srcFile,srcBytes,srcMod, ...
    'VariableNames',{'animal','file','bytes','modified'});
info.validation = rep;
info.trim = trimInfo;

if isempty(outFile)
    outFile = fullfile(manifest.outDir, sprintf(spec.groupPattern, manifest.group));
end
info.outFile = string(outFile);

if verbose
    fprintf('  -> %d rows, %d cells, %d animals | %s\n', ...
        info.nRows, info.nCells, info.nAnimals, spec.timeVar);
end

if dryRun
    if verbose; fprintf('  dryRun: nothing written.\n'); end
    return
end

%% ---- write ----
if ~isfolder(manifest.outDir); mkdir(manifest.outDir); end
S = struct(spec.varname, T, 'groupInfo', info);
save(outFile, '-struct','S', '-v7.3');

% manifest .json alongside, for review/diff
[outDir,outName] = fileparts(outFile);
jsonPath = fullfile(outDir,[outName '_manifest.json']);
fid = fopen(jsonPath,'w');
if fid > 0
    fprintf(fid,'%s', jsonencode(manifest,'PrettyPrint',true));
    fclose(fid);
end

if verbose
    fprintf('  wrote %s\n  wrote %s\n', outFile, jsonPath);
end
end

%% ================= helpers =================
function m = normalizeManifest(m)
req = {'group','family'};
missing = req(~isfield(m,req));
if ~isempty(missing)
    error('aggregateStimGroup:badManifest','Manifest is missing: %s', strjoin(missing,', '));
end
m.group  = char(m.group);
m.family = char(m.family);
if ~isfield(m,'tableDir') || isempty(m.tableDir); m.tableDir = '.'; end
if isfield(m,'animals'); m.animals = string(m.animals(:))'; end
if isfield(m,'files') && ~isempty(m.files)
    if ischar(m.files); m.files = {m.files}; end
    m.files = cellstr(string(m.files));
end
if ~isfield(m,'outDir'); m.outDir = pwd; end
m.outDir = char(m.outDir);
end

function [files, animals] = resolveFiles(m, spec)
if isfield(m,'files') && ~isempty(m.files)
    files = m.files(:)';
    if isfield(m,'animals') && numel(m.animals) == numel(files)
        animals = m.animals;
    else
        animals = strings(1,numel(files));
        for k = 1:numel(files)
            [~,nm] = fileparts(files{k});
            tok = regexp(nm,'[A-Z]{2}\d{4}','match','once');
            if isempty(tok); tok = nm; end
            animals(k) = string(tok);
        end
    end
    return
end
if ~isfield(m,'animals') || isempty(m.animals)
    error('aggregateStimGroup:badManifest','Manifest needs either animals or files.');
end
if ~isfield(m,'cohortRoot')
    error('aggregateStimGroup:badManifest','Manifest needs cohortRoot when animals are given.');
end
animals = m.animals;
files = cell(1,numel(animals));
for k = 1:numel(animals)
    a = char(animals(k));
    files{k} = fullfile(char(m.cohortRoot), a, m.tableDir, [a spec.suffixProcessed]);
end
end

function sha = repoSHA()
sha = "";
try
    here = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    [st,out] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', here));
    if st == 0; sha = string(strtrim(out)); end
catch
end
end
