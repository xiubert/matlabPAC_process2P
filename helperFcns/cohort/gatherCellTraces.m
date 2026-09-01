function [M,cells,info] = gatherCellTraces(T,rowMask,traceVar,timeVar,varargin)
% gatherCellTraces  Stack per-cell traces for a row selection into a matrix.
%
%   [M,cells,info] = gatherCellTraces(T,rowMask,traceVar,timeVar)
%   [...] = gatherCellTraces(...,'idVars',{'animal','roiID'})
%
%   Turns "the rows of this group table matching some condition" into an
%   nCells x nFrames matrix ready for cohortMeanSEM, handling the three ways
%   that step goes wrong: empty selections, ragged trace lengths, and cells
%   contributing more than one row.
%
%   Inputs:
%     T        - group table.
%     rowMask  - logical mask or row indices selecting the condition
%                (e.g. T.BPNdBAmpl == 50). Empty selects nothing.
%     traceVar - trace column to gather. Entries may be nReps x nFrames
%                (per-trial) or 1 x nFrames (already cell-averaged); per-trial
%                entries are averaged across reps so the result is always one
%                row per cell.
%     timeVar  - time-axis column, used to return info.t and to verify the
%                axis is common across the selection.
%
%   Name/Value:
%     'idVars' - columns identifying a unique cell. Default {'animal','roiID'}.
%
%   Outputs:
%     M     - nCells x nFrames. Empty (0 x 0) when nothing was selected.
%     cells - table of cell identities, one row per row of M, in the same order.
%     info  - struct:
%               .t              1 x nFrames time vector (empty if no rows)
%               .ok             false when M is empty
%               .reason         why M is empty, '' otherwise
%               .nRowsSelected  rows matched by rowMask
%               .nDuplicateCells cells that contributed >1 row and were averaged
%
%   Empty selections are returned as ok=false with a reason rather than
%   raised, because a legitimate group can have zero cells passing a
%   significance filter -- 4 of 13 TOMT animals have zero both-contrast
%   significant CGC cells. Callers render a labeled empty panel. Passing such
%   a selection to splitapply instead errors with "Group numbers must be a
%   vector of positive integers", which is the failure this replaces.
%
%   Ragged trace widths DO raise: they mean the group was never validated,
%   and silently trimming would misalign cells in time.
%
%   See also cohortMeanSEM, groupN, validateStimGroup.

p = inputParser;
addParameter(p,'idVars',{'animal','roiID'},@(x) iscellstr(x) || isstring(x)); %#ok<ISCLSTR>
parse(p,varargin{:});
idVars = cellstr(p.Results.idVars);

need = [idVars, {traceVar, timeVar}];
missing = setdiff(need, T.Properties.VariableNames);
if ~isempty(missing)
    error('gatherCellTraces:missingVars','Table is missing column(s): %s', strjoin(missing,', '));
end

if islogical(rowMask); rows = find(rowMask); else; rows = rowMask(:); end
info.nRowsSelected   = numel(rows);
info.nDuplicateCells = 0;
info.t      = [];
info.ok     = false;
info.reason = '';

if isempty(rows)
    M = []; cells = T([],idVars);
    info.reason = 'no rows matched the selection';
    return
end

S = T(rows,:);

% ---- time axis must be common across the selection ----
tw = cellfun(@numel, S.(timeVar));
if numel(unique(tw)) > 1
    error('gatherCellTraces:raggedTime', ...
        ['%s widths differ across the selection (%s). The group is not on a ' ...
         'common time base -- run validateStimGroup.'], ...
        timeVar, mat2str(unique(tw)'));
end
tAll = cell2mat(cellfun(@(c) reshape(c,1,[]), S.(timeVar), 'uni',0));
if max(abs(tAll - tAll(1,:)),[],'all') > 1e-9
    error('gatherCellTraces:timeAxisMismatch', ...
        '%s values differ across the selection -- run validateStimGroup.', timeVar);
end
info.t = tAll(1,:);

% ---- trace widths must match the time axis ----
trw = cellfun(@(c) size(c,2), S.(traceVar));
if numel(unique(trw)) > 1
    error('gatherCellTraces:raggedTrace', ...
        '%s widths differ across the selection (%s).', traceVar, mat2str(unique(trw)'));
end
if trw(1) ~= tw(1)
    error('gatherCellTraces:traceTimeMismatch', ...
        '%s is %d frames but %s is %d frames.', traceVar, trw(1), timeVar, tw(1));
end

% ---- one row per cell (average reps, then average duplicate rows) ----
key = string(S.(idVars{1}));
for k = 2:numel(idVars)
    key = strcat(key,'_',string(S.(idVars{k})));
end
[uKey,ia] = unique(key,'stable');
nCell = numel(uKey);

M = nan(nCell, trw(1));
for c = 1:nCell
    idx = find(key == uKey(c));
    if numel(idx) > 1; info.nDuplicateCells = info.nDuplicateCells + 1; end
    % average reps within each row, then across duplicate rows for this cell
    perRow = cell2mat(cellfun(@(x) mean(x,1,'omitnan'), S.(traceVar)(idx), 'uni',0));
    M(c,:) = mean(perRow,1,'omitnan');
end

cells = S(ia, idVars);
info.ok = true;
end
