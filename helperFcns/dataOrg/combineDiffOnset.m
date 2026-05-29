function mergedTable = combineDiffOnset(anmlROIbyStim, baselineSec)
% combineDiffOnset  Merge anmlROIbyStim rows that share all stim params
%                   except onset, vertcat'ing per-trial cell columns so
%                   the merged row holds all repetitions across onsets.
%
% Usage:
%   mergedTable = combineDiffOnset(anmlROIbyStim)               % baselineSec = 1
%   mergedTable = combineDiffOnset(anmlROIbyStim, baselineSec)
%
% Inputs:
%   anmlROIbyStim - long-form table from anmlROIbyStimTable. Must include:
%                     * an *sOnset column (e.g. BPNsOnset). Exact case.
%                     * frameRate (scalar per row)
%                     * raw F cell columns from this set as available:
%                       rawFroi, moCorRawFroi, fissaFroi, SCALEDfissaFroi
%                   Optional columns produced by an upstream per-row
%                   dFoFcalc: t_dFF, dFF. If present, they are trimmed
%                   to common minLength before the merge.
%   baselineSec   - pre-onset baseline length in seconds. Used to
%                   onset-align the raw F cells so frame 1 of each row
%                   sits at "baselineSec before that row's own onset",
%                   matching the dFF coordinate system upstream
%                   dFoFcalc produces. Default 1.
%
% Output:
%   mergedTable - one row per unique (non-onset) stim-param combination.
%                 *Froi cells hold (nReps * nMergedOnsets) x nFrames
%                 matrices, onset-aligned. If dFF / t_dFF were present
%                 on input, they are length-equalized and vertcat'd
%                 too. The retained scalar *sOnset is the smallest of
%                 the merged onsets; treat it as metadata only — both
%                 dFF and raw F are in onset-normalized coords where
%                 the stim sits at frame frameRate*baselineSec for
%                 every trial.
%
% Notes:
%   - Group key = every non-cell column except the *sOnset column,
%     'Pulse', and 'stimID' (those split otherwise-identical rows).
%   - frameRate is taken per-row; the function does not assume a
%     common frameRate across the table.
%   - Errors on missing or ambiguous *sOnset column. Legacy tables
%     with lowercase 'sonset' should be migrated via
%     migrateBPNStimTableFields first.

arguments
    anmlROIbyStim table
    baselineSec   (1,1) double {mustBePositive} = 1
end

allVars = anmlROIbyStim.Properties.VariableNames;
isCell  = varfun(@iscell, anmlROIbyStim, 'OutputFormat', 'uniform');

% Locate the *sOnset column (BPNsOnset, FRAsOnset, etc.) — exact case.
onsetVars = allVars(endsWith(allVars,'sOnset'));
if isempty(onsetVars)
    error('combineDiffOnset:NoOnsetColumn', ...
        'No *sOnset column found in input table.');
elseif numel(onsetVars) > 1
    error('combineDiffOnset:AmbiguousOnsetColumn', ...
        'Multiple *sOnset columns found (%s); ambiguous.', strjoin(onsetVars,', '));
end
onsetVar = onsetVars{1};

%% (a) trim t_dFF / dFF to common minLength — only if present
if ismember('t_dFF', allVars) && ismember('dFF', allVars)
    t_lengths = cellfun(@length, anmlROIbyStim.t_dFF);
    minLen_dFF = min(t_lengths);
    anmlROIbyStim.t_dFF = cellfun(@(x) x(1:minLen_dFF), ...
        anmlROIbyStim.t_dFF, 'UniformOutput', false);
    anmlROIbyStim.dFF = cellfun(@(x) x(:, 1:minLen_dFF), ...
        anmlROIbyStim.dFF, 'UniformOutput', false);
end

%% (b) onset-align raw F cells, then trim to common minLength
rawFvars = intersect({'rawFroi','moCorRawFroi','fissaFroi','SCALEDfissaFroi'}, ...
                    allVars, 'stable');
if ~isempty(rawFvars)
    nRows = height(anmlROIbyStim);
    startFrames = zeros(nRows,1);
    for r = 1:nRows
        nFrames = size(anmlROIbyStim.(rawFvars{1}){r}, 2);
        t = (1:nFrames) / anmlROIbyStim.frameRate(r);
        sf = find(t > anmlROIbyStim.(onsetVar)(r) - baselineSec, 1, 'first');
        if isempty(sf), sf = 1; end
        startFrames(r) = sf;
    end
    for v = 1:numel(rawFvars)
        varName = rawFvars{v};
        anmlROIbyStim.(varName) = arrayfun(@(r) ...
            anmlROIbyStim.(varName){r}(:, startFrames(r):end), ...
            (1:nRows)', 'UniformOutput', false);
    end
    lens = cellfun(@(x) size(x,2), anmlROIbyStim.(rawFvars{1}));
    minLen_raw = min(lens);
    for v = 1:numel(rawFvars)
        varName = rawFvars{v};
        anmlROIbyStim.(varName) = cellfun(@(x) x(:, 1:minLen_raw), ...
            anmlROIbyStim.(varName), 'UniformOutput', false);
    end
end

%% (c) group by non-cell, non-onset, non-Pulse, non-stimID; vertcat per group
nonCellVars = allVars(~isCell);
groupVars = setdiff(nonCellVars, {onsetVar,'Pulse','stimID'}, 'stable');

[uniqueGroups, ~, groupIdx] = unique(anmlROIbyStim(:, groupVars), 'rows');
numGroups = size(uniqueGroups, 1);

mergedTable = table();
for i = 1:numGroups
    subTable = anmlROIbyStim(groupIdx == i, :);

    % keep the row with the smallest onset as base
    [~, minIdx] = min(subTable.(onsetVar));
    newRow = subTable(minIdx, :);

    % vertcat multi-row cell contents (raw F + dFF) across the group
    cellVars = allVars(isCell);
    for j = 1:numel(cellVars)
        varName = cellVars{j};
        sampleData = subTable.(varName){1};
        if size(sampleData, 1) > 1
            combinedData = vertcat(subTable.(varName){:});
            newRow.(varName){1} = combinedData;
        end
    end

    mergedTable = [mergedTable; newRow]; %#ok<AGROW>
end
end
