function out = tableRLF(T,varargin)
% TABLERLF  Compute response-level functions for all cells in a stim table.
%
%   out = tableRLF(T)
%   out = tableRLF(T,'cellIDvars',{'animal','roiID'},'dBlist',[],'nConsec',3)
%
%   Inputs:
%       T   - anmlROIbyStim-style table with columns dBampl, sig, pkResp,
%             plus any columns identifying unique cells (default
%             {'animal','roiID'}).
%
%   Name/Value:
%       'cellIDvars' - cellstr of columns that uniquely identify a cell.
%                      Default {'animal','roiID'}.
%       'dBlist'     - dB levels to align to. Default: sort(unique(T.dBampl)).
%       'nConsec'    - min consecutive sig==1 levels for inclusion. Default 3.
%
%   Output (struct):
%       .dBlist        1 x nDB dB axis
%       .cellInfo      table: one row per cell, with cellID vars +
%                      threshold + included flag + nSig
%       .RLFall        nCellAll x nDB matrix of pkResp (all cells)
%       .sigAll        nCellAll x nDB logical matrix
%       .RLFincl       nCellIncl x nDB matrix (included cells only)
%       .meanRLF       1 x nDB mean across included cells
%       .semRLF        1 x nDB SEM across included cells
%       .nIncluded     scalar count of included cells
%       .nTotal        scalar count of cells in T
%
%   See also cellRLF, plotRLF.

p = inputParser;
addRequired(p,'T',@istable);
addParameter(p,'cellIDvars',{'animal','roiID'},@(x) iscellstr(x) || isstring(x));
addParameter(p,'dBlist',[],@(x) isnumeric(x) && (isempty(x) || isvector(x)));
addParameter(p,'nConsec',3,@(x) isnumeric(x) && isscalar(x) && x>=1);
parse(p,T,varargin{:});
T          = p.Results.T;
cellIDvars = cellstr(p.Results.cellIDvars);
dBlist     = p.Results.dBlist;
nConsec    = p.Results.nConsec;

if isempty(dBlist)
    dBlist = reshape(unique(T.dBampl),1,[]);
else
    dBlist = reshape(sort(dBlist),1,[]);
end
nDB = numel(dBlist);

[gID,gT] = findgroups(T(:,cellIDvars));
nCell = height(gT);

RLFall  = nan(nCell,nDB);
sigAll  = false(nCell,nDB);
included  = false(nCell,1);
threshold = nan(nCell,1);
nSig      = zeros(nCell,1);

for c = 1:nCell
    cellOut = cellRLF(T(gID==c,:),'dBlist',dBlist,'nConsec',nConsec);
    RLFall(c,:)  = cellOut.RLF;
    sigAll(c,:)  = cellOut.sigAtDB;
    included(c)  = cellOut.included;
    threshold(c) = cellOut.threshold;
    nSig(c)      = sum(cellOut.sigAtDB);
end

cellInfo = gT;
cellInfo.included  = included;
cellInfo.threshold = threshold;
cellInfo.nSig      = nSig;

RLFincl = RLFall(included,:);

out.dBlist    = dBlist;
out.cellInfo  = cellInfo;
out.RLFall    = RLFall;
out.sigAll    = sigAll;
out.RLFincl   = RLFincl;
out.meanRLF   = mean(RLFincl,1,'omitnan');
out.semRLF    = SEMcalc(RLFincl,1);
out.nIncluded = sum(included);
out.nTotal    = nCell;
end
