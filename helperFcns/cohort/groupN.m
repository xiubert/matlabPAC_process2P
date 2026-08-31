function N = groupN(T,varargin)
% groupN  Cell and animal counts for a group table, with per-animal breakdown.
%
%   N = groupN(T)
%   N = groupN(T,rowMask)
%   N = groupN(...,'idVars',{'animal','roiID'})
%
%   Every cohort figure needs to state how many cells and how many animals it
%   summarises. This computes that once, including the per-animal breakdown,
%   so a group dominated by one animal is visible rather than hidden behind a
%   total.
%
%   Inputs:
%     T       - group table (anmlROIbyStim-style) with animal + roiID columns.
%     rowMask - (optional) logical mask or row indices to restrict to. Use
%               this for counts of a filtered subset, e.g. cells that passed
%               a significance criterion.
%
%   Name/Value:
%     'idVars' - columns identifying a unique cell. Default {'animal','roiID'}.
%
%   Output (struct):
%     .nRows        rows in the (masked) table
%     .nCells       unique cells
%     .nAnimals     unique animals
%     .animals      string array of animal IDs
%     .perAnimal    table: animal, nCells, nRows  (sorted by animal)
%     .label        e.g. "n = 75 cells / 4 mice"
%     .singleAnimal true when nAnimals == 1
%     .empty        true when nRows == 0
%
%   Cell identity is always (animal, roiID) because roiID collides across
%   animals -- every animal numbers its ROIs from 1.
%
%   See also cohortMeanSEM, gatherCellTraces, annotateN.

p = inputParser;
addRequired(p,'T',@istable);
addOptional(p,'rowMask',[],@(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p,'idVars',{'animal','roiID'},@(x) iscellstr(x) || isstring(x)); %#ok<ISCLSTR>
parse(p,T,varargin{:});
T       = p.Results.T;
rowMask = p.Results.rowMask;
idVars  = cellstr(p.Results.idVars);

if ~isempty(rowMask); T = T(rowMask,:); end

missing = setdiff(idVars, T.Properties.VariableNames);
if ~isempty(missing)
    error('groupN:missingVars','Table is missing id column(s): %s', strjoin(missing,', '));
end

N.nRows = height(T);
if N.nRows == 0
    N.nCells = 0; N.nAnimals = 0;
    N.animals = strings(0,1);
    N.perAnimal = table(strings(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',{'animal','nCells','nRows'});
    N.label = 'n = 0 cells';
    N.singleAnimal = false;
    N.empty = true;
    return
end

anml = string(T.(idVars{1}));
key  = anml;
for k = 2:numel(idVars)
    key = strcat(key, '_', string(T.(idVars{k})));
end

N.animals  = unique(anml);
N.nAnimals = numel(N.animals);
N.nCells   = numel(unique(key));

nC = zeros(N.nAnimals,1); nR = zeros(N.nAnimals,1);
for a = 1:N.nAnimals
    m = anml == N.animals(a);
    nC(a) = numel(unique(key(m)));
    nR(a) = sum(m);
end
N.perAnimal = table(N.animals, nC, nR, 'VariableNames',{'animal','nCells','nRows'});

N.singleAnimal = N.nAnimals == 1;
N.empty = false;
N.label = sprintf('n = %d cell%s / %d %s', ...
    N.nCells, plural(N.nCells), N.nAnimals, mouseWord(N.nAnimals));
end

% ---- helpers ----
function s = plural(n)
if n == 1; s = ''; else; s = 's'; end
end

function s = mouseWord(n)
if n == 1; s = 'mouse'; else; s = 'mice'; end
end
