function testConsensusROIsets()
% testConsensusROIsets  Unit test for voting across per-tif ROI sets.
% The headless path segments every tif separately and keeps the cells enough
% tifs agree on. This pins the vote semantics -- union at one extreme, strict
% inner join at the other -- and the per-cell detection count that the single
% session-mean pass cannot provide.
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
testConfig();

sz = [64 64];
%three ground-truth cells; A is in every set, B in most, C in one
A = disc(sz,[20 20],5);
B = disc(sz,[20 45],5);
C = disc(sz,[45 20],5);

mk = @(varargin) labelImg2moCorROI(stack(sz,varargin{:}),...
    'edgeMarginPx',0,'minAreaPx',5);

sets = { mk(A,B), mk(A,B), mk(A,B,C), mk(A), mk(A,B) };   % A:5 B:4 C:1

% CASE 1: a low threshold keeps everything seen twice; the singleton drops
[roi2,info2] = consensusROIsets(sets,'minVotes',2);
assert(numel(roi2)==2,'minVotes=2 should keep A and B, got %d',numel(roi2));
assert(info2.nClusters==3,'expected 3 clusters, got %d',info2.nClusters);
assert(all(ismember([roi2.votes],[4 5])),'votes not carried onto the ROIs');

% CASE 2: minVotes=1 is the UNION -- every cluster, singleton included
roi1 = consensusROIsets(sets,'minVotes',1);
assert(numel(roi1)==3,'minVotes=1 (union) should keep all 3, got %d',numel(roi1));

% CASE 3: minVotes = nSets is the strict inner join -- only the cell in every
% set survives. This is why the pipeline does not default to it: one bad tif
% governs the whole result.
roiAll = consensusROIsets(sets,'minVotes',numel(sets));
assert(numel(roiAll)==1,'strict join should keep only A, got %d',numel(roiAll));
assert(roiAll(1).votes==numel(sets));

% CASE 4: a fraction <= 1 is read as a fraction of the set count
roiFrac = consensusROIsets(sets,'minVotes',0.5);      % -> 3 of 5 votes
assert(numel(roiFrac)==2,'minVotes=0.5 should keep A and B, got %d',numel(roiFrac));
assert(isequal(sort([roiFrac.votes]),[4 5]));

% CASE 5: votes and voteFraction agree, and are ordered with the ROIs
for k = 1:numel(roi2)
    assert(abs(roi2(k).voteFraction - roi2(k).votes/numel(sets)) < 1e-12,...
        'voteFraction inconsistent with votes');
end

% CASE 6: cells that merely touch are not merged. Two discs 11 px apart
% overlap in neither set, so they must stay two clusters.
near1 = disc(sz,[32 28],5);
near2 = disc(sz,[32 39],5);
setsNear = { mk(near1,near2), mk(near1,near2) };
roiNear = consensusROIsets(setsNear,'minVotes',2);
assert(numel(roiNear)==2,'adjacent cells were merged (%d clusters)',numel(roiNear));

% CASE 7: asking for more votes than there are sets is an error, not a silent
% clamp to the strict join
threwHigh = false;
try
    consensusROIsets(sets,'minVotes',99);
catch ME
    threwHigh = strcmp(ME.identifier,'consensusROIsets:minVotesTooHigh');
end
assert(threwHigh,'minVotes above the set count must raise minVotesTooHigh');

% ...but a threshold that is reachable yet met by nothing warns and returns
% empty rather than erroring
setsThin = { mk(A), mk(B), mk(C) };      % no cell appears twice
w = warning('off','consensusROIsets:noneKept');
c = onCleanup(@() warning(w)); %#ok<NASGU>
roiNone = consensusROIsets(setsThin,'minVotes',2);
assert(isempty(roiNone),'unmet minVotes should give an empty set');

% CASE 8: sets of differing frame size are refused, not silently voted on
bad = { mk(A), labelImg2moCorROI(stack([32 32],disc([32 32],[16 16],4)),...
    'edgeMarginPx',0,'minAreaPx',5) };
threw = false;
try
    consensusROIsets(bad,'minVotes',1);
catch ME
    threw = strcmp(ME.identifier,'consensusROIsets:sizeMismatch');
end
assert(threw,'mismatched frame sizes must raise consensusROIsets:sizeMismatch');

disp('testConsensusROIsets PASS: union/join extremes, vote counts, no false merges');
end

% ------------------------------------------------------------------------
function m = disc(sz,c,r)
[X,Y] = meshgrid(1:sz(2),1:sz(1));
m = (Y-c(1)).^2 + (X-c(2)).^2 <= r^2;
end
% ------------------------------------------------------------------------
function L = stack(sz,varargin)
L = zeros(sz);
for k = 1:numel(varargin)
    L(varargin{k}) = k;
end
end
