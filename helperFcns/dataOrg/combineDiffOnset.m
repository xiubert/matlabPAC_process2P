function mergedTable = combineDiffOnset(anmlROIbyStim)
% Combine rows in anmlROIbyStim that have the same sound but different onsets
% stim params should all be scalar variables
% fluorescence traces-related data should all be cell variables
% Run this before calculating peak response

t_lengths = cellfun(@(x) length(x), anmlROIbyStim.t_dFF);
minLength = min(t_lengths);
anmlROIbyStim.t_dFF = cellfun(@(x) x(1:minLength), ...
    anmlROIbyStim.t_dFF, 'UniformOutput', false);
anmlROIbyStim.dFF = cellfun(@(x) x(:, 1:minLength), ...
    anmlROIbyStim.dFF, 'UniformOutput', false);

% Define grouping variables (all non-cell columns except the *sOnset column, 'Pulse' and 'stimID')
allVars = anmlROIbyStim.Properties.VariableNames;
isCell = varfun(@iscell, anmlROIbyStim, 'OutputFormat', 'uniform');
nonCellVars = allVars(~isCell);

% Locate the *sOnset column (BPNsOnset, FRAsOnset, etc.). Polymorphic so
% this function works across stim families that follow the
% <StimFamily><unit><Quantity> naming convention.
onsetVars = allVars(endsWith(allVars,'sOnset'));
if isempty(onsetVars)
    error('combineDiffOnset:NoOnsetColumn', ...
        'No *sOnset column found in input table.');
elseif numel(onsetVars) > 1
    error('combineDiffOnset:AmbiguousOnsetColumn', ...
        'Multiple *sOnset columns found (%s); ambiguous.', strjoin(onsetVars,', '));
end
onsetVar = onsetVars{1};

groupVars = setdiff(nonCellVars, {onsetVar,'Pulse','stimID'}, 'stable');

% Find unique groups and their row indices
[uniqueGroups, ~, groupIdx] = unique(anmlROIbyStim(:, groupVars), 'rows');
numGroups = size(uniqueGroups, 1);

% Initialize a new table to store merged results
mergedTable = table();

for i = 1:numGroups
    % Find all rows belonging to this specific group
    rows = find(groupIdx == i);
    subTable = anmlROIbyStim(rows, :);
    
    % Find the row with the smallest onset
    [~, minIdx] = min(subTable.(onsetVar));

    % Create the base row for this group (using the smallest-onset row)
    newRow = subTable(minIdx, :);
    
    % Combine the cell columns for all rows in this group
    cellVars = allVars(isCell);
    for j = 1:numel(cellVars)
        varName = cellVars{j};
        % Inspect the first available entry in this group
        sampleData = subTable.(varName){1};
        
        % Combine ONLY if the matrix has multiple rows (rawFroi,moCorRawFroi,fissaFroi,SCALEDfissaFroi,dFF)
        if size(sampleData, 1) > 1
            combinedData = vertcat(subTable.(varName){:});
            newRow.(varName){1} = combinedData;
        end
    end
    
    mergedTable = [mergedTable; newRow]; %#ok<AGROW>
end