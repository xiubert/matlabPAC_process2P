function dPrime = cellFRA2dPrime(FRA,varargin)
p = inputParser;
addRequired(p,'FRA',@ismatrix)
addOptional(p,'nAmpl',3,@(x) isnumeric(x) && isscalar(x));
addOptional(p,'nFreq',28,@(x) isnumeric(x) && isscalar(x));
addOptional(p,'nBoots',1000,@isnumeric);

parse(p,FRA,varargin{:});
FRA = p.Results.FRA;
boots = p.Results.nBoots;

%assumes more freq than Ampl if input is 2d
if ~isvector(FRA)
    nFreq = max(size(FRA));
    nAmpl = min(size(FRA));
else
    nFreq = p.Results.nFreq;
    nAmpl = p.Results.nAmpl;
end
[~,maxRespLinIDX] = max(FRA);
[dBidx,fqIDX] = ind2sub([nAmpl nFreq],maxRespLinIDX);

%get dBidx
tmp = [dBidx+1 dBidx dBidx-1];
dBadjIDX = zero2nan(tmp.*ismember(tmp,[1:nAmpl]));
%get fqIDX
tmp = [fqIDX+1 fqIDX-1];
fqAdjIDX = zero2nan(tmp.*ismember(tmp,[1:nFreq]));

adjRespIDX = [sub2ind([nAmpl nFreq],dBadjIDX,repmat(fqIDX,[1 3])),...
sub2ind([nAmpl nFreq],repmat(dBidx,[1 2]),fqAdjIDX)];

adjRespIDX = adjRespIDX(~isnan(adjRespIDX));
nResp = length(adjRespIDX);

adjRespMean = nanmean(FRA(adjRespIDX));

tmpDp = zeros(1,boots);
for b = 1:boots
    tmpDp(b) = adjRespMean-nanmean(FRA(randi(nFreq*nAmpl,[1 nResp])));
end
dPrime = nanmean(tmpDp);
    