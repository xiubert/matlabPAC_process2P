function out = sigCellCounts(T,family,varargin)
% sigCellCounts  How many cells respond significantly at each stimulus level.
%
%   out = sigCellCounts(T,family)
%   out = sigCellCounts(T,family,'rowMask',mask)
%
%   Family-agnostic tally of significant cells per level of the family's
%   stimulus axis: sound level for BPN, contrast for CGC. Also reports how
%   many cells are significant at EVERY level, which for CGC is the
%   "significant in both contrasts" criterion that gates the scatter and the
%   paired test.
%
%   Inputs:
%     T      - group table (or any per-animal table of the same family).
%     family - 'BPN' | 'CGC'.
%
%   Name/Value:
%     'rowMask' - restrict to a subset of rows before tallying.
%     'idVars'  - cell identity columns. Default from stimGroupSpec.
%
%   Output (struct):
%     .family     the family
%     .levels     1 x nLevels stimulus axis, ascending
%     .labels     1 x nLevels display labels ("50 dB", "Low contrast", ...)
%     .levelVar   name of the level column
%     .nTotal     cells in the (masked) table
%     .nSig       1 x nLevels count significant at each level
%     .pctSig     1 x nLevels percent of nTotal
%     .nSigAll    cells significant at EVERY level AND present at every level
%     .pctSigAll  percent of nTotal
%     .nSigAny    cells significant at one or more levels
%     .pctSigAny  percent of nTotal
%     .sigMatrix  nCells x nLevels logical
%     .present    nCells x nLevels logical (false where a level is missing)
%     .cells      table of cell identities, rows aligned with sigMatrix
%
%   nSigAll requires presence as well as significance, so a cell missing a
%   level cannot be counted as significant everywhere by default.
%
%   A cell contributing more than one row at a level counts as significant if
%   ANY of them is (gatherCellValues 'reduce','any'), matching how
%   plotCGCgroup builds its significance matrix.
%
%   See also gatherCellValues, groupN, stimGroupSpec, plotCGCgroup.

spec = stimGroupSpec(family);

p = inputParser;
addParameter(p,'rowMask',[],@(x) isempty(x)||islogical(x)||isnumeric(x));
addParameter(p,'idVars',spec.idVars,@(x) iscellstr(x)||isstring(x)); %#ok<ISCLSTR>
parse(p,varargin{:});
rowMask = p.Results.rowMask;
idv     = cellstr(p.Results.idVars);

if ~isempty(rowMask); T = T(rowMask,:); end

out.family   = spec.family;
out.levelVar = spec.levelVar;

if height(T) == 0
    out.levels = []; out.labels = strings(1,0);
    out.nTotal = 0; out.nSig = []; out.pctSig = [];
    out.nSigAll = 0; out.pctSigAll = NaN; out.nSigAny = 0; out.pctSigAny = NaN;
    out.sigMatrix = false(0,0); out.present = false(0,0);
    out.cells = T([],idv);
    return
end

levels = sort(unique(T.(spec.levelVar)))';
nL     = numel(levels);

allKeys = unique(cellKey(T,idv),'stable');
nCell   = numel(allKeys);
S       = false(nCell,nL);
present = false(nCell,nL);

for k = 1:nL
    mask = T.(spec.levelVar) == levels(k);
    [v,cellsK,gi] = gatherCellValues(T,mask,spec.sigVar,'idVars',idv,'reduce','any');
    if ~gi.ok; continue; end
    [tf,loc] = ismember(cellKey(cellsK,idv),allKeys);
    S(loc(tf),k)       = v(tf) == 1;
    present(loc(tf),k) = true;
end

out.levels    = levels;
out.labels    = levelLabels(spec,levels);
out.nTotal    = nCell;
out.nSig      = sum(S,1);
out.pctSig    = 100 * out.nSig / max(nCell,1);
out.nSigAll   = sum(all(S,2) & all(present,2));
out.pctSigAll = 100 * out.nSigAll / max(nCell,1);
out.nSigAny   = sum(any(S,2));
out.pctSigAny = 100 * out.nSigAny / max(nCell,1);
out.sigMatrix = S;
out.present   = present;
out.cells     = firstRowsFor(T,idv,allKeys);
end

%% ---- helpers ----
function k = cellKey(T,idVars)
k = string(T.(idVars{1}));
for i = 2:numel(idVars); k = strcat(k,'_',string(T.(idVars{i}))); end
end

function C = firstRowsFor(T,idVars,keys)
k = cellKey(T,idVars);
[~,loc] = ismember(keys,k);
C = T(loc,idVars);
end

function lab = levelLabels(spec,levels)
lab = strings(1,numel(levels));
switch spec.family
    case 'CGC'
        % dBdelta is the DRC's dB RANGE, so smaller = lower contrast.
        for i = 1:numel(levels); lab(i) = sprintf('%g dB range',levels(i)); end
        if numel(levels) == 2
            lab(1) = "Low contrast"; lab(2) = "High contrast";
        end
    otherwise
        for i = 1:numel(levels); lab(i) = sprintf('%g dB',levels(i)); end
end
end
