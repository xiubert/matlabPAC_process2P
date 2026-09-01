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
%   Inputs:
%     dataPath - animal data folder (the same one processAnimal2P used).
%
%   Name/Value:
%     'families'  - families to consider. Default: every family in
%                   stimGroupSpec that has a _raw table in dataPath.
%     'showPlots' - leave each script's QC figures open. Default false;
%                   figures the scripts create are closed again afterwards so
%                   a batch run does not accumulate dozens of windows.
%     'overwrite' - re-run a family whose processed table already exists.
%                   Default true (the _raw input is never modified, so
%                   re-running is safe and repeatable).
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
addParameter(p,'verbose',true,@islogical);
parse(p,dataPath,varargin{:});
dataPath  = char(p.Results.dataPath);
families  = p.Results.families;
showPlots = p.Results.showPlots;
overwrite = p.Results.overwrite;
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
    spec    = stimGroupSpec(families{i});
    rawFile = fullfile(dataPath,[animal spec.suffixRaw]);
    proc    = fullfile(dataPath,[animal spec.suffixProcessed]);
    rec = struct('family',spec.family,'ok',false,'skipped',true,'reason','', ...
        'rawFile',rawFile,'procFile',proc,'script',spec.processScript);

    if ~isfile(rawFile)
        rec.reason = 'no _raw table (this animal has no tifs of this family)';
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
        runFamilyScript(spec.processScript, dataPath, animal);
        rec.ok = isfile(proc);
        if ~rec.ok
            rec.reason = sprintf('%s ran but wrote no %s', spec.processScript, proc);
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
function runFamilyScript(scriptName, dataPath, animal) %#ok<INUSD>
% Run the script in this isolated workspace. dataPath and animal are the two
% variables the process* scripts look for; everything else they create stays
% local to this function and never reaches the caller.
scriptPath = which(scriptName);
if isempty(scriptPath)
    error('processAnimalStimFamilies:scriptNotFound', ...
        '%s is not on the MATLAB path.', scriptName);
end
run(scriptPath);
end
