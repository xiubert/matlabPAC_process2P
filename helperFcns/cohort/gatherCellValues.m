function [v,cells,info] = gatherCellValues(T,rowMask,valueVar,varargin)
% gatherCellValues  One scalar per cell for a row selection.
%
%   [v,cells,info] = gatherCellValues(T,rowMask,valueVar)
%   [...] = gatherCellValues(...,'reduce','any','idVars',{'animal','roiID'})
%
%   Scalar counterpart to gatherCellTraces, for per-cell quantities such as
%   peak response and significance flags. Same contract: empty selections
%   report a reason rather than raising, and the result is always one value
%   per cell.
%
%   Inputs:
%     T        - group table.
%     rowMask  - logical mask or row indices selecting the condition.
%     valueVar - column holding one value per row. May be numeric, logical,
%                or a cell column of scalars (pkResp / sig / pkPT / sigPk are
%                all stored as cells of scalars).
%
%   Name/Value:
%     'reduce' - how to combine multiple rows for the same cell:
%                  'mean'  (default) numeric mean, NaN-omitting
%                  'any'   true if any row is true      (significance: OR)
%                  'all'   true only if every row is true (sig in ALL conditions)
%                  'first' first row in table order
%     'idVars' - columns identifying a unique cell. Default {'animal','roiID'}.
%
%   Outputs:
%     v     - nCells x 1. Empty when nothing was selected.
%     cells - table of cell identities, one row per element of v.
%     info  - struct: .ok .reason .nRowsSelected .nDuplicateCells .reduce
%
%   'any' vs 'all' matters: a cell "responds to this level" if any row at that
%   level is significant, but "responds in both contrasts" requires all. Making
%   the caller name it prevents the two being silently interchanged.
%
%   See also gatherCellTraces, cohortMeanSEM, groupN.

p = inputParser;
addParameter(p,'reduce','mean',@(x) any(strcmpi(x,{'mean','any','all','first'})));
addParameter(p,'idVars',{'animal','roiID'},@(x) iscellstr(x) || isstring(x)); %#ok<ISCLSTR>
parse(p,varargin{:});
red    = lower(p.Results.reduce);
idVars = cellstr(p.Results.idVars);

need = [idVars, {valueVar}];
missing = setdiff(need, T.Properties.VariableNames);
if ~isempty(missing)
    error('gatherCellValues:missingVars','Table is missing column(s): %s', strjoin(missing,', '));
end

if islogical(rowMask); rows = find(rowMask); else; rows = rowMask(:); end
info = struct('ok',false,'reason','','nRowsSelected',numel(rows), ...
    'nDuplicateCells',0,'reduce',red);

if isempty(rows)
    v = []; cells = T([],idVars);
    info.reason = 'no rows matched the selection';
    return
end

S = T(rows,:);

% ---- unwrap to a numeric column ----
col = S.(valueVar);
if iscell(col)
    n = cellfun(@numel,col);
    if any(n > 1)
        error('gatherCellValues:notScalar', ...
            '%s holds %d non-scalar entr(y/ies) in the selection (sizes %s).', ...
            valueVar, sum(n>1), mat2str(unique(n)'));
    end
    raw = nan(numel(col),1);
    hasVal = n == 1;                      % empty cells stay NaN
    raw(hasVal) = cellfun(@(x) double(x(1)), col(hasVal));
else
    raw = double(col(:));
end

% ---- one value per cell ----
key = string(S.(idVars{1}));
for k = 2:numel(idVars)
    key = strcat(key,'_',string(S.(idVars{k})));
end
[uKey,ia] = unique(key,'stable');
v = nan(numel(uKey),1);
for c = 1:numel(uKey)
    idx = find(key == uKey(c));
    if numel(idx) > 1; info.nDuplicateCells = info.nDuplicateCells + 1; end
    x = raw(idx);
    switch red
        case 'mean';  v(c) = mean(x,'omitnan');
        case 'any';   v(c) = double(any(x == 1));
        case 'all';   v(c) = double(all(x == 1));
        case 'first'; v(c) = x(1);
    end
end

cells = S(ia, idVars);
info.ok = true;
end
