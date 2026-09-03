function out = processAnimalStimFamilies(dataPath,varargin)
% processAnimalStimFamilies  Run the per-stim process* script for every family
%                            this animal has, turning _raw tables into
%                            processed ones.
%
%   out = processAnimalStimFamilies(dataPath)
%   out = processAnimalStimFamilies(dataPath,'showPlots',true)
%   out = processAnimalStimFamilies(dataPath,'families',{'BPN'})
%
%   stimParam2ROI writes one <animal>_anmlROI_<Family>stimTable_raw.mat per
%   stimulus family. Those carry stimulus-aligned traces but no dF/F and no
%   peak responses; the process* script for each family adds those and writes
%   the processed table that aggregateStimGroup reads. This runs whichever of
%   those scripts the animal actually needs, so the per-animal path completes
%   in one pass instead of relying on the user to remember the step.
%
%   dF/F and peak responses are computed HERE, per animal, because baselines
%   and significance are per-animal quantities -- they cannot be deferred to
%   the grouping step.
%
%   Not every family has a _raw stage. FRA is driven straight from the tif
%   inventory (processFRA -> FRAmap), so its stimGroupSpec entry has an empty
%   suffixRaw and its precondition is <animal>_tifFileList.mat containing map
%   tifs, not a _raw table. The gate below checks whichever applies.
%
%   Inputs:
%     dataPath - animal data folder (the same one processAnimal2P used).
%
%   Name/Value:
%     'families'  - families to consider. Default: every family in
%                   stimGroupSpec whose input for this animal is present --
%                   a _raw table, or map tifs for a family without one.
%     'showPlots' - leave each script's QC figures open. Default false;
%                   figures the scripts create are closed again afterwards so
%                   a batch run does not accumulate dozens of windows.
%     'overwrite' - re-run a family whose processed table already exists.
%                   Default true (the _raw input is never modified, so
%                   re-running is safe and repeatable).
%     'scriptVars'- struct of variables to define in the script's workspace
%                   before it runs, for per-animal parameters the scripts read
%                   from the workspace. e.g.
%                     struct('PTfreqSelect',6484)
%                   for an animal recorded with several pure-tone frequencies.
%     'runLabel'  - name of the analysis run. When given, a runInfo struct is
%                   appended to every processed table this call writes,
%                   recording the run, the ROI set it came from (count and
%                   whether those ROIs were hand-drawn or segmented), the time
%                   and the repo SHA. aggregateStimGroup reads it, so a table
%                   left behind by an earlier run can be caught instead of
%                   being aggregated as if it were current. Default '' = no
%                   stamp.
%     'verbose'   - print progress. Default true.
%
%   Output (struct array), one element per family considered:
%     .family .ok .skipped .reason .rawFile .procFile .script
%
%   A family whose script errors is reported and the others still run, so one
%   bad family does not cost the rest. Nothing downstream is weakened by
%   that: aggregateStimGroup still refuses an animal whose processed table is
%   missing (aggregateStimGroup:notProcessed).
%
%   The scripts run inside this function's workspace, so their many
%   intermediate variables do not leak into the caller's.
%
%   NOTE: each process* script has an editable parameter block at its top
%   (baselineSec, tBasePT, pkPTsigSD, ...). Running them from here uses those
%   defaults, which are the conventions stamped by stimGroupSpec. To analyse
%   an animal with different parameters, run that script by hand instead.
%
%   See also processAnimal2P, stimParam2ROI, processBPN2P, processCGC,
%   aggregateStimGroup, stimGroupSpec.

p = inputParser;
addRequired(p,'dataPath',@(x) ischar(x)||isstring(x));
addParameter(p,'families',{},@(x) isempty(x)||iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
addParameter(p,'showPlots',false,@islogical);
addParameter(p,'overwrite',true,@islogical);
addParameter(p,'scriptVars',struct(),@isstruct);
addParameter(p,'runLabel','',@(x) ischar(x)||isstring(x));
addParameter(p,'verbose',true,@islogical);
parse(p,dataPath,varargin{:});
dataPath  = char(p.Results.dataPath);
families  = p.Results.families;
showPlots = p.Results.showPlots;
overwrite = p.Results.overwrite;
scriptVars= p.Results.scriptVars;
runLabel  = char(p.Results.runLabel);
verbose   = p.Results.verbose;

if ~isfolder(dataPath)
    error('processAnimalStimFamilies:noDataPath','Not a folder: %s', dataPath);
end
if isempty(families); families = stimGroupSpec(); else; families = cellstr(families); end

animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
if isempty(animal)
    error('processAnimalStimFamilies:noAnimalID', ...
        'Could not read an animal ID (e.g. AA0001) from: %s', dataPath);
end

out = struct('family',{},'ok',{},'skipped',{},'reason',{}, ...
    'rawFile',{},'procFile',{},'script',{});

for i = 1:numel(families)
    spec = stimGroupSpec(families{i});
    proc = fullfile(dataPath,[animal spec.suffixProcessed]);

    % A family with an empty suffixRaw has no _raw stage: its process* script
    % is driven straight from the tif inventory (FRA -> FRAmap). Gating those
    % on a _raw file builds the path <dataPath>/<animal>, which never exists,
    % so the family was silently skipped for every animal AND told the caller
    % it had no tifs of that family -- which for an animal with map tifs is
    % simply false. Check the precondition each family actually has.
    if isempty(spec.suffixRaw)
        srcFile = fullfile(dataPath,[animal '_tifFileList.mat']);
        [haveSrc,whyNot] = localHasMapTifs(srcFile);
    else
        srcFile = fullfile(dataPath,[animal spec.suffixRaw]);
        haveSrc = isfile(srcFile);
        whyNot  = 'no _raw table (this animal has no tifs of this family)';
    end

    rec = struct('family',spec.family,'ok',false,'skipped',true,'reason','', ...
        'rawFile',srcFile,'procFile',proc,'script',spec.processScript);

    if ~haveSrc
        rec.reason = whyNot;
        out(end+1) = rec; %#ok<AGROW>
        if verbose; fprintf('  %-4s skipped: %s\n', spec.family, rec.reason); end
        continue
    end
    if isfile(proc) && ~overwrite
        rec.reason = 'processed table already exists (overwrite=false)';
        out(end+1) = rec; %#ok<AGROW>
        if verbose; fprintf('  %-4s skipped: %s\n', spec.family, rec.reason); end
        continue
    end

    rec.skipped = false;
    if verbose; fprintf('  %-4s running %s ...\n', spec.family, spec.processScript); end

    figsBefore = findall(0,'Type','figure');
    if ~showPlots
        vis0 = get(0,'DefaultFigureVisible');
        set(0,'DefaultFigureVisible','off');
    end
    try
        runFamilyScript(spec.processScript, dataPath, animal, scriptVars);
        rec.ok = isfile(proc);
        if ~rec.ok
            rec.reason = sprintf('%s ran but wrote no %s', spec.processScript, proc);
        elseif ~isempty(runLabel)
            stampRunInfo(proc, runLabel, dataPath, animal, spec.family);
        end
    catch ME
        rec.reason = sprintf('%s failed: %s', spec.processScript, ME.message);
        warning('processAnimalStimFamilies:familyFailed', ...
            '%s: %s\nOther families still run; this one has no processed table.', ...
            spec.family, rec.reason);
    end
    if ~showPlots
        close(setdiff(findall(0,'Type','figure'), figsBefore));
        set(0,'DefaultFigureVisible',vis0);
    end

    out(end+1) = rec; %#ok<AGROW>
    if verbose && rec.ok
        fprintf('  %-4s wrote %s\n', spec.family, [animal spec.suffixProcessed]);
    end
end

if verbose
    nOK = sum([out.ok]);
    fprintf('%s: %d family/families processed', animal, nOK);
    bad = out(~[out.ok] & ~[out.skipped]);
    if ~isempty(bad)
        fprintf(', %d FAILED (%s)', numel(bad), strjoin({bad.family},', '));
    end
    fprintf('\n');
end
end

%% ---- helper ----
function [ok,why] = localHasMapTifs(listFile)
% Precondition for a family with no _raw stage: the tif inventory exists and
% actually contains map tifs. Reported separately from "file missing" so the
% skip reason says which of the two is true.
ok = false;
if ~isfile(listFile)
    [~,nm,ext] = fileparts(listFile);
    why = sprintf('no tif inventory (%s not found)',[nm ext]);
    return
end
S = load(listFile,'tifFileList');
if ~isfield(S,'tifFileList') || ~isfield(S.tifFileList,'map') || isempty(S.tifFileList.map)
    why = 'tif inventory has no map tifs (this animal has none of this family)';
    return
end
ok = true; why = '';
end

function runFamilyScript(scriptName, dataPath, animal, scriptVars) %#ok<INUSD>
% Run the script in this isolated workspace. dataPath and animal are the two
% variables the process* scripts look for; everything else they create stays
% local to this function and never reaches the caller. scriptVars fields are
% defined here too, so a script that reads a per-animal parameter from the
% workspace (processCGC's PTfreqSelect) picks it up.
fn = fieldnames(scriptVars);
for k = 1:numel(fn)
    eval([fn{k} ' = scriptVars.(fn{k});']); %#ok<EVLEQ>
end
scriptPath = which(scriptName);
if isempty(scriptPath)
    error('processAnimalStimFamilies:scriptNotFound', ...
        '%s is not on the MATLAB path.', scriptName);
end
run(scriptPath);
end


function stampRunInfo(procFile, runLabel, dataPath, animal, family)
% Record which analysis run wrote this table, and off which ROI set.
%
% Every run writes the same filenames into the same animal folder, so a family
% that fails leaves the PREVIOUS run's table in place -- and aggregation,
% having no way to tell, treats it as current. On the TOMT cohort that quietly
% produced a group file mixing hand-drawn and segmented ROIs. This stamp is
% what lets aggregateStimGroup('requireRun',...) refuse that.
info = struct();
info.runLabel  = char(runLabel);
info.family    = char(family);
info.animal    = char(animal);
info.dataPath  = char(dataPath);
info.when      = string(datetime('now','Format','uuuu-MM-dd HH:mm:ss'));
info.roiSource = 'unknown';
info.nROI      = NaN;
info.roiFile   = '';
info.repoSHA   = "";

d = dir(fullfile(dataPath,[animal '_moCorrROI_*.mat']));
if ~isempty(d)
    info.roiFile = d(1).name;
    try
        S = load(fullfile(d(1).folder,d(1).name));
        if isfield(S,'moCorROI'); info.nROI = numel(S.moCorROI); end
        if isfield(S,'roiParams') && isfield(S.roiParams,'source')
            info.roiSource = char(S.roiParams.source);
        elseif isfield(S,'cellposeParams')
            info.roiSource = 'cellpose';
        else
            info.roiSource = 'handDrawn';
        end
    catch
    end
end

try
    here = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    [st,out] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null',here));
    if st == 0; info.repoSHA = string(strtrim(out)); end
catch
end

runInfo = info; %#ok<NASGU>
save(procFile,'runInfo','-append');
end
