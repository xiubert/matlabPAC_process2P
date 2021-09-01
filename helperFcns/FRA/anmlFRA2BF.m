function BFuDB = anmlFRA2BF(FRA,varargin)
p = inputParser;
addRequired(p,'FRA',@isstruct)
addOptional(p,'F',FRA.CellSigPkLinDBfreq,...
    @(x) ismatrix(x) && all(size(x)>1));

parse(p,FRA,varargin{:});
FRA = p.Results.FRA;
F = p.Results.F;

nCell = size(F,1);
[~,BFfqIDX] = max(cell2mat(cellfun(@(c) nanmean(reshape(c,[numel(FRA.dBlist) numel(FRA.freqList)]),1),...
    mat2cell(F,repmat(1,[1 nCell]),numel(FRA.freqList)*numel(FRA.dBlist)),'uni',0)),[],2);

BFuDB = FRA.freqList(BFfqIDX);
end