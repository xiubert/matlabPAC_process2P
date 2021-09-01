function dPrime = anmlFRA2dPrime(FRA,varargin)
%if working from TtrfAnml table
% FRA = table2struct(TtrfAnml(anmlRow,:))

p = inputParser;
addRequired(p,'FRA',@isstruct)
addOptional(p,'nBoots',1000,@isnumeric);
addOptional(p,'F',FRA.CellSigPkLinDBfreq,...
    @(x) ismatrix(x) && all(size(x)>1));

parse(p,FRA,varargin{:});
FRA = p.Results.FRA;
F = p.Results.F;
boots = p.Results.nBoots;
nAmpl = length(FRA.dBlist);
nFreq = length(FRA.freqList);

[~,maxRespLinIDX] = max(F,[],2);
[dBidx,fqIDX] = ind2sub([nAmpl nFreq],maxRespLinIDX);

nCell = size(F,1);
%get adjacent freq/ampl indices
tmp = [dBidx+1 dBidx dBidx-1];
dBadjIDX = zero2nan(tmp.*ismember(tmp,[1:nAmpl]));
tmp = [fqIDX+1 fqIDX-1];
fqAdjIDX = zero2nan(tmp.*ismember(tmp,[1:nFreq]));

%back to linear IDX
adjRespIDX = [sub2ind([nAmpl nFreq],dBadjIDX,repmat(fqIDX,[1 3])),...
sub2ind([nAmpl nFreq],repmat(dBidx,[1 2]),fqAdjIDX)];

nResp = sum(~isnan(adjRespIDX),2);

adjRespIDX = num2cell(adjRespIDX,2);
adjRespIDX = cellfun(@(c) c(~isnan(c)),adjRespIDX,'uni',0);
adjRespMean = arrayfun(@(r) nanmean(F(r,adjRespIDX{r})),...
    [1:nCell],'uni',1)';

tmpDp = zeros(nCell,boots);
for b = 1:boots
    randRespIDX = arrayfun(@(r) randi(nAmpl*nFreq,[1 nResp(r)]),[1:nCell],'uni',0)';
    adjRandRespMean = arrayfun(@(r) nanmean(F(r,randRespIDX{r})),...
    [1:nCell],'uni',1)';
    tmpDp(:,b) = adjRespMean-adjRandRespMean;
end

dPrime = nanmean(tmpDp,2);



