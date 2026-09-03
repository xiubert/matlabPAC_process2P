function cmp = compareROIsets(roiA,roiB,opts)
% compareROIsets  Match two moCorROI sets by mask overlap and report agreement.
%
%   cmp = compareROIsets(roiA,roiB)
%   cmp = compareROIsets(roiA,roiB,'iouThreshold',0.3)
%
%   Pairs ROIs between two sets drawn on the SAME field (e.g. hand-drawn vs
%   segmented) by intersection-over-union, using an optimal assignment rather
%   than greedy nearest-neighbour so one large mask cannot hoover up several
%   small ones. Written for calibrating a segmenter against a hand-drawn
%   corpus, but it is just as useful for comparing two segmentation settings.
%
%   A note on how to read the result: an unmatched ROI in set B is NOT
%   automatically a false positive. If set A is hand-drawn it is very likely
%   to UNDER-count, so B-only ROIs are candidate cells the human skipped. The
%   informative comparison is .areaRatio and .medianIoU on the MATCHED pairs
%   (which is what a size parameter should be tuned on), plus whether the
%   B-only masks look like the matched ones (.areaOnlyB vs .areaB) or like
%   debris.
%
%   Inputs
%     roiA, roiB  moCorROI struct arrays with a logical .mask field, both on
%                 the same frame size.
%
%   Name-value
%     'iouThreshold'  minimum IoU for a pair to count as matched. Default 0.3.
%     'verbose'       print the summary. Default false.
%
%   Output (struct)
%     .pairs        nMatch x 2, [indexInA indexInB]
%     .iou          nMatch x 1 IoU of each pair
%     .medianIoU    median over matched pairs
%     .nA .nB .nMatched .onlyA .onlyB   (counts and index vectors)
%     .areaA .areaB                     areas of the matched masks
%     .areaRatio    median areaB/areaA over matched pairs -- the number a
%                   dilation setting is tuned to bring to 1
%     .areaOnlyB    areas of the unmatched B masks
%     .iouMatrix    nA x nB full IoU matrix
%
%   See also labelImg2moCorROI, cellposeROIset, intersectROIfiles

arguments
    roiA struct
    roiB struct
    opts.iouThreshold (1,1) double {mustBeInRange(opts.iouThreshold,0,1)} = 0.3
    opts.verbose      (1,1) logical = false
end

nA = numel(roiA);
nB = numel(roiB);
if nA == 0 || nB == 0
    error('compareROIsets:emptySet','Both ROI sets must be non-empty (%d vs %d).',nA,nB);
end
if ~isequal(size(roiA(1).mask),size(roiB(1).mask))
    error('compareROIsets:sizeMismatch',...
        'Masks are %s vs %s; the two sets are not on the same frame.',...
        mat2str(size(roiA(1).mask)),mat2str(size(roiB(1).mask)));
end

%--- IoU matrix, bounding boxes first so we only intersect plausible pairs -
boxA = arrayfun(@(r) maskBox(r.mask),roiA,'uni',0);
boxB = arrayfun(@(r) maskBox(r.mask),roiB,'uni',0);
areaA = arrayfun(@(r) nnz(r.mask),roiA);
areaB = arrayfun(@(r) nnz(r.mask),roiB);

iouMatrix = zeros(nA,nB);
for i = 1:nA
    for j = 1:nB
        if ~boxesOverlap(boxA{i},boxB{j}), continue, end
        inter = nnz(roiA(i).mask & roiB(j).mask);
        if inter == 0, continue, end
        iouMatrix(i,j) = inter/(areaA(i)+areaB(j)-inter);
    end
end

%--- optimal assignment ---------------------------------------------------
%matchpairs minimises cost; unmatched costs 1 - threshold so any pair below
%the threshold is never worth matching
cost = 1 - iouMatrix;
unmatchedCost = 1 - opts.iouThreshold;
M = matchpairs(cost,unmatchedCost);

if isempty(M)
    pairs = zeros(0,2);
    iou   = zeros(0,1);
else
    idx   = sub2ind([nA nB],M(:,1),M(:,2));
    good  = iouMatrix(idx) >= opts.iouThreshold;
    pairs = M(good,:);
    iou   = iouMatrix(idx(good));
end

onlyA = setdiff(1:nA,pairs(:,1));
onlyB = setdiff(1:nB,pairs(:,2));

cmp = struct();
cmp.pairs      = pairs;
cmp.iou        = iou(:);
cmp.medianIoU  = median(iou);
cmp.nA         = nA;
cmp.nB         = nB;
cmp.nMatched   = size(pairs,1);
cmp.onlyA      = onlyA(:)';
cmp.onlyB      = onlyB(:)';
cmp.areaA      = areaA(:)';
cmp.areaB      = areaB(:)';
cmp.areaOnlyB  = areaB(onlyB);
cmp.iouMatrix  = iouMatrix;
cmp.iouThreshold = opts.iouThreshold;
if isempty(pairs)
    cmp.areaRatio = NaN;
    cmp.medianIoU = NaN;
else
    cmp.areaRatio = median(areaB(pairs(:,2))./areaA(pairs(:,1)));
end

if opts.verbose
    fprintf(['compareROIsets: A=%d B=%d matched=%d (IoU>=%.2f) '...
        'medianIoU=%.3f areaRatio(B/A)=%.2f\n'],...
        nA,nB,cmp.nMatched,opts.iouThreshold,cmp.medianIoU,cmp.areaRatio);
    if ~isempty(onlyB)
        fprintf('  B-only: %d ROIs, median area %.0f px (matched B median %.0f px)\n',...
            numel(onlyB),median(cmp.areaOnlyB),median(areaB(pairs(:,2))));
    end
    if ~isempty(onlyA)
        fprintf('  A-only: %d ROIs (IDs %s)\n',numel(onlyA),...
            strjoin({roiA(onlyA).ID},','));
    end
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
