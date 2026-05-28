function out = cellRLF(cellT,varargin)
% CELLRLF  Build a response-level function for a single cell.
%
%   out = cellRLF(cellT)
%   out = cellRLF(cellT,'dBlist',dBlist,'nConsec',3)
%
%   Inputs:
%       cellT   - subtable of an anmlROIbyStim-style table containing all
%                 rows for ONE cell (one animal+roiID). Must have columns
%                 'BPNdBAmpl', 'sig', and 'pkResp'. sig and pkResp may be
%                 numeric or scalar-in-cell.
%
%   Name/Value:
%       'dBlist'  - dB levels to align the RLF to. Default: sort(unique(cellT.BPNdBAmpl)).
%                   Levels missing from cellT are filled with NaN.
%       'nConsec' - minimum number of consecutive sig==1 levels required for
%                   the cell to be included. Default: 3.
%
%   Output (struct):
%       .dBlist     1 x nDB vector of dB levels
%       .RLF        1 x nDB vector of pkResp at each dB (NaN if level absent)
%       .sigAtDB    1 x nDB logical vector of sig at each dB
%       .included   true if cell has >=nConsec consecutive sig==1 levels
%       .threshold  lowest dB with sig==1 anywhere on the dB axis (NOT
%                   constrained to fall inside the qualifying consecutive
%                   run). NaN if the cell is not included.
%
%   Note: if cellT contains multiple rows at the same dB (e.g. replicate
%   presentations or multiple stim conditions per dB), pkResp is averaged
%   across rows and sig is OR'd (sig==1 if any row at that dB is
%   significant). A warning is issued when this occurs.
%
%   See also tableRLF, plotRLF.

p = inputParser;
addRequired(p,'cellT',@istable);
addParameter(p,'dBlist',[],@(x) isnumeric(x) && (isempty(x) || isvector(x)));
addParameter(p,'nConsec',3,@(x) isnumeric(x) && isscalar(x) && x>=1);
parse(p,cellT,varargin{:});
cellT   = p.Results.cellT;
dBlist  = p.Results.dBlist;
nConsec = p.Results.nConsec;

if isempty(dBlist)
    dBlist = reshape(unique(cellT.BPNdBAmpl),1,[]);
else
    dBlist = reshape(sort(dBlist),1,[]);
end

sigVec    = unwrapNumCol(cellT.sig);
pkVec     = unwrapNumCol(cellT.pkResp);
dBcellVec = cellT.BPNdBAmpl(:);

nDB = numel(dBlist);
RLF     = nan(1,nDB);
sigAtDB = false(1,nDB);

dupDB = [];
for k = 1:nDB
    idx = find(dBcellVec==dBlist(k));
    if isempty(idx); continue; end
    if numel(idx)>1
        dupDB(end+1) = dBlist(k); %#ok<AGROW>
        RLF(k)     = mean(pkVec(idx),'omitnan');
        sigAtDB(k) = any(sigVec(idx)==1);
    else
        RLF(k)     = pkVec(idx);
        sigAtDB(k) = sigVec(idx)==1;
    end
end
if ~isempty(dupDB)
    warning('cellRLF:duplicateDB',...
        'Multiple rows at dB level(s) [%s]; pkResp averaged and sig OR''d across rows.',...
        num2str(dupDB));
end

included = maxRunLength(sigAtDB) >= nConsec;

% threshold: lowest sig==1 dB on the axis (NOT constrained to qualifying run).
% NaN for excluded cells.
if included
    threshold = dBlist(find(sigAtDB,1,'first'));
else
    threshold = NaN;
end

out.dBlist    = dBlist;
out.RLF       = RLF;
out.sigAtDB   = sigAtDB;
out.included  = included;
out.threshold = threshold;
end

% ---- helpers ----
function v = unwrapNumCol(c)
% Accept either a numeric column or a cell column of scalars.
if iscell(c)
    v = cellfun(@(x) double(x), c);
else
    v = double(c);
end
v = v(:);
end

function r = maxRunLength(b)
% Longest run of true values in logical vector b.
b = logical(b(:)');
if ~any(b); r = 0; return; end
d = diff([0 b 0]);
runStarts = find(d==1);
runEnds   = find(d==-1)-1;
r = max(runEnds-runStarts+1);
end
