function [roiOut,info] = consensusROIsets(roiSets,opts)
% consensusROIsets  Combine per-tif ROI sets into one voted consensus set.
%
%   [roiOut,info] = consensusROIsets(roiSets)
%   [roiOut,info] = consensusROIsets(roiSets,'minVotes',0.5,'iouThreshold',0.3)
%
%   Segmenting each tif's mean projection separately gives one ROI set per
%   tif, each a noisy sample of the same underlying cells. This clusters those
%   sets by spatial overlap and keeps the cells that enough tifs agree on,
%   producing a single moCorROI set with a detection count per cell.
%
%   Why not just segment the session mean? On this data Cellpose-SAM's
%   detection is strongly dependent on the granular noise texture of a single
%   tif's projection: averaging tifs removes that texture and detection
%   collapses non-monotonically (on TO0003: 1 tif -> 33 cells, 2 -> 19,
%   4 -> 0, all 29 -> 0), while each individual tif's detections agree well
%   with hand-drawn ROIs (median IoU 0.6-0.7). Voting across the per-tif
%   segmentations keeps the regime where the segmenter works and averages
%   away its instability, instead of averaging away the signal it needs.
%
%   The vote threshold spans the obvious combining rules: minVotes = 1 is the
%   union, minVotes = numel(roiSets) is the strict inner join, and the useful
%   settings are in between -- a strict join is dominated by the worst tif,
%   which may have found almost nothing.
%
%   Inputs
%     roiSets  cell array of moCorROI struct arrays (from labelImg2moCorROI),
%              all on the same frame size. Empty sets are allowed and count
%              towards the total.
%
%   Name-value
%     'minVotes'      how many sets must detect a cell to keep it. A value
%                     <= 1 is read as a FRACTION of the number of sets
%                     (default 0.5); > 1 as an absolute count.
%     'iouThreshold'  overlap needed to call two ROIs the same cell.
%                     Default 0.3.
%     'maskQuantile'  a consensus cell's mask is the pixels claimed by at
%                     least this fraction of the sets that detected it.
%                     Default 0.5 (per-pixel majority).
%     'minAreaPx'     drop consensus masks below this. Default 20.
%     'verbose'       Default false.
%
%   Outputs
%     roiOut  moCorROI struct array, IDs renumbered 1..N in the usual
%             row-major centroid order, each with an extra .votes field (how
%             many sets found it) and .voteFraction.
%     info    .nSets .nClusters .nKept .votes (per kept cell)
%             .voteHistogram .droppedArea .droppedVotes
%
%   See also labelImg2moCorROI, cellposeSegment, cellposeROIset,
%   compareROIsets

arguments
    roiSets            cell
    opts.minVotes     (1,1) double {mustBePositive} = 0.5
    opts.iouThreshold (1,1) double {mustBeInRange(opts.iouThreshold,0,1)} = 0.3
    opts.maskQuantile (1,1) double {mustBeInRange(opts.maskQuantile,0,1)} = 0.5
    opts.minAreaPx    (1,1) double {mustBeNonnegative} = 20
    opts.verbose      (1,1) logical = false
end

nSets = numel(roiSets);
if nSets == 0
    error('consensusROIsets:noSets','roiSets is empty.');
end

%--- frame size, from the first non-empty set -----------------------------
sz = [];
for s = 1:nSets
    if ~isempty(roiSets{s})
        sz = size(roiSets{s}(1).mask);
        break
    end
end
if isempty(sz)
    error('consensusROIsets:allEmpty','Every ROI set is empty.');
end
for s = 1:nSets
    if ~isempty(roiSets{s}) && ~isequal(size(roiSets{s}(1).mask),sz)
        error('consensusROIsets:sizeMismatch',...
            'Set %d has %s masks but set 1 has %s.',...
            s,mat2str(size(roiSets{s}(1).mask)),mat2str(sz));
    end
end

if opts.minVotes <= 1
    minVotes = max(1,round(opts.minVotes*nSets));
else
    minVotes = round(opts.minVotes);
end

%--- cluster ROIs across sets by overlap ----------------------------------
%process the richest sets first so clusters are seeded from the tifs the
%segmenter did best on, rather than from whichever happened to be first
[~,order] = sort(cellfun(@numel,roiSets),'descend');

cSum   = {};       % per-cluster pixel vote accumulator (uint16)
cMask  = {};       % current consensus mask, for matching
cVotes = [];
cBox   = {};

for s = order(:)'
    set_ = roiSets{s};
    claimed = false(1,numel(cVotes));
    for r = 1:numel(set_)
        m   = set_(r).mask;
        box = maskBox(m);
        a   = nnz(m);

        best = 0; bestIoU = opts.iouThreshold;
        for c = 1:numel(cVotes)
            if claimed(c) || ~boxesOverlap(box,cBox{c}), continue, end
            inter = nnz(m & cMask{c});
            if inter == 0, continue, end
            iou = inter/(a + nnz(cMask{c}) - inter);
            if iou >= bestIoU
                bestIoU = iou;
                best = c;
            end
        end

        if best > 0
            cSum{best}   = cSum{best} + uint16(m);
            cVotes(best) = cVotes(best) + 1; %#ok<AGROW>
            cMask{best}  = cSum{best} >= max(1,ceil(opts.maskQuantile*cVotes(best)));
            cBox{best}   = maskBox(cMask{best});
            claimed(best)= true;
        else
            cSum{end+1}   = uint16(m);   %#ok<AGROW>
            cVotes(end+1) = 1;           %#ok<AGROW>
            cMask{end+1}  = m;           %#ok<AGROW>
            cBox{end+1}   = box;         %#ok<AGROW>
            claimed(end+1)= true;        %#ok<AGROW>
        end
    end
end

nClusters = numel(cVotes);

%--- keep what enough sets agree on --------------------------------------
keep = false(1,nClusters);
finalMask = cell(1,nClusters);
droppedArea = 0;
for c = 1:nClusters
    if cVotes(c) < minVotes, continue, end
    m = cSum{c} >= max(1,ceil(opts.maskQuantile*cVotes(c)));
    %the per-pixel majority can fragment when member masks disagree at the
    %edges; keep the largest piece so mask and outline describe one object
    CC = bwconncomp(m,8);
    if CC.NumObjects > 1
        [~,big] = max(cellfun(@numel,CC.PixelIdxList));
        m = false(sz);
        m(CC.PixelIdxList{big}) = true;
    end
    if nnz(m) < opts.minAreaPx
        droppedArea = droppedArea + 1;
        continue
    end
    keep(c) = true;
    finalMask{c} = m;
end

keepIdx = find(keep);
info = struct('nSets',nSets,'nClusters',nClusters,'nKept',numel(keepIdx),...
    'minVotes',minVotes,'iouThreshold',opts.iouThreshold,...
    'votes',cVotes(keepIdx),'voteHistogram',histcounts(cVotes,0.5:1:(nSets+0.5)),...
    'droppedArea',droppedArea,'droppedVotes',nClusters-numel(keepIdx)-droppedArea);

if isempty(keepIdx)
    roiOut = labelImg2moCorROI(zeros(sz));
    warning('consensusROIsets:noneKept',...
        'No cluster reached %d of %d votes (best was %d).',...
        minVotes,nSets,max(cVotes));
    return
end

%--- rebuild as a label image so the moCorROI construction stays in one place
labelImg = zeros(sz);
for n = 1:numel(keepIdx)
    labelImg(finalMask{keepIdx(n)}) = n;
end
%overlapping consensus masks would have overwritten each other above; the
%later cell wins, and any cell fully overwritten simply disappears
[roiOut,~] = labelImg2moCorROI(labelImg,'minAreaPx',opts.minAreaPx,'edgeMarginPx',0);

%carry the votes onto the ROIs, matching by pixel identity
votesKept = cVotes(keepIdx);
for k = 1:numel(roiOut)
    lbl = mode(labelImg(roiOut(k).mask));
    roiOut(k).votes        = votesKept(lbl);
    roiOut(k).voteFraction = votesKept(lbl)/nSets;
end

if opts.verbose
    fprintf(['consensusROIsets: %d sets -> %d clusters -> %d kept '...
        '(>=%d votes; %d below threshold, %d too small)\n'],...
        nSets,nClusters,numel(roiOut),minVotes,info.droppedVotes,droppedArea);
end
end

% ------------------------------------------------------------------------
function b = maskBox(m)
rws = find(any(m,2));
cls = find(any(m,1));
if isempty(rws)
    b = [NaN NaN NaN NaN];
else
    b = [rws(1) rws(end) cls(1) cls(end)];
end
end
% ------------------------------------------------------------------------
function tf = boxesOverlap(a,b)
tf = ~(any(isnan(a)) || any(isnan(b))) && ...
     a(1) <= b(2) && b(1) <= a(2) && a(3) <= b(4) && b(3) <= a(4);
end
