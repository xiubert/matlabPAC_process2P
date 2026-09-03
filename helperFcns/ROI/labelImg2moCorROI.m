function [moCorROI,info] = labelImg2moCorROI(labelImg,opts)
% labelImg2moCorROI  Turn a segmentation label image into a moCorROI struct.
%
%   moCorROI = labelImg2moCorROI(labelImg)
%   [moCorROI,info] = labelImg2moCorROI(labelImg,'dilatePx',1,...)
%
%   Builds the SAME struct array TIFcatROIgui writes when a human draws ROIs,
%   so a segmented ROI set drops into the rest of the pipeline unchanged
%   (moCorRawF2tifList, intersectROIfiles, remapROItoAcq, the FISSA driver).
%   labelImg is a 2-D integer-labelled image -- 0 = background, each positive
%   integer one cell -- e.g. Cellpose's <name>_cp_masks.tif.
%
%   Inputs
%     labelImg  H x W numeric/integer label image.
%
%   Name-value
%     'dilatePx'     grow every mask by this many pixels (disk strel) before
%                    anything else. Default 0. Segmenters trace the soma
%                    boundary tightly; hand-drawn ellipses are more generous,
%                    and mask size sets both F and neuropil contamination, so
%                    this is the knob that makes dF/F amplitudes comparable
%                    to a hand-drawn corpus. Calibrate it, do not guess --
%                    see etc/calibrateCellposeROI.m.
%     'preventOverlap'  after dilation, drop pixels claimed by more than one
%                    ROI so masks stay disjoint the way the label image was.
%                    Default true. Only matters when dilatePx > 0.
%     'minAreaPx'    drop masks smaller than this. Default 20.
%     'maxAreaPx'    drop masks larger than this. Default Inf. Catches the
%                    occasional merged blob.
%     'edgeMarginPx' drop masks whose bounding box comes within this many
%                    pixels of the frame edge. Default 0. Set it to the
%                    neuropil annulus width so FISSA's ring is not clipped
%                    at the border -- the same concern remapROItoAcq's
%                    neuropilMarginPx encodes.
%     'startID'      first numeric ID. Default 1.
%
%   Outputs
%     moCorROI  1 x nROI struct array with the fields TIFcatROIgui emits:
%               ID ('1','2',...), pos ([x y w h] bounding box), XYvertices
%               (nPts x 2 [x y] traced outline), frame, deleted, type
%               ('cellpose'), mask (logical H x W), ROIxyCoord (2 x nPts,
%               row1 = x), ROIcurveOrderedXY (2 x nPts, closed contour),
%               label (.String / .Position).
%               Ordered by centroid, row-major, so the same label image
%               always yields the same ID assignment.
%     info      struct: .nLabels .nKept .dropped (per-reason counts)
%               .areaPx .equivDiamPx (kept ROIs) .droppedIDs
%
%   Why the outline comes from bwboundaries and NOT orderEllipsePtOnCurve:
%   that helper angular-sorts the perimeter points, which is only valid for a
%   star-convex blob (its own header says so). A segmented mask is often not,
%   and the sort then folds the polygon into a self-intersecting shape. That
%   polygon is exactly what the FISSA driver reads out of ROIcurveOrderedXY to
%   build each cell's neuropil annulus, so a folded contour would silently
%   sample the wrong neuropil. bwboundaries traces a real closed contour.
%
%   See also cellposeSegment, cellposeROIset, TIFcatROIgui, mask2polyCoord,
%   remapROItoAcq, moCorRawF2tifList

arguments
    labelImg                        {mustBeNumeric,mustBeNonempty}
    opts.dilatePx       (1,1) double {mustBeNonnegative,mustBeInteger} = 0
    opts.preventOverlap (1,1) logical = true
    opts.minAreaPx      (1,1) double {mustBeNonnegative} = 20
    opts.maxAreaPx      (1,1) double {mustBePositive}    = Inf
    opts.edgeMarginPx   (1,1) double {mustBeNonnegative,mustBeInteger} = 0
    opts.startID        (1,1) double {mustBeInteger}     = 1
end

if ndims(labelImg) > 2 %#ok<ISMAT>
    error('labelImg2moCorROI:not2D',...
        'labelImg must be 2-D (got %s). Segment one plane at a time.',...
        mat2str(size(labelImg)));
end

labelImg = double(labelImg);
labelIDs = unique(labelImg(:));
labelIDs = labelIDs(labelIDs > 0);
nLabels  = numel(labelIDs);

info = struct('nLabels',nLabels,'nKept',0,...
    'dropped',struct('empty',0,'minArea',0,'maxArea',0,'edge',0),...
    'areaPx',[],'equivDiamPx',[],'droppedIDs',[]);

if nLabels == 0
    moCorROI = emptyROIstruct();
    warning('labelImg2moCorROI:noLabels','Label image contains no ROIs.');
    return
end

[H,W] = size(labelImg);

%--- masks, optionally grown ----------------------------------------------
masks = false(H,W,nLabels);
for k = 1:nLabels
    masks(:,:,k) = labelImg == labelIDs(k);
end

if opts.dilatePx > 0
    se = strel('disk',opts.dilatePx);
    for k = 1:nLabels
        masks(:,:,k) = imdilate(masks(:,:,k),se);
    end
    if opts.preventOverlap
        %a pixel two grown masks both claim belongs to neither: the label
        %image had disjoint cells and F extraction assumes it still does
        contested = sum(masks,3) > 1;
        if any(contested(:))
            for k = 1:nLabels
                mk = masks(:,:,k);
                mk(contested) = false;
                masks(:,:,k) = mk;
            end
        end
    end
end

%--- QC + geometry --------------------------------------------------------
keep       = false(nLabels,1);
centroids  = nan(nLabels,2);   % [x y]
areaPx     = zeros(nLabels,1);
droppedIDs = [];

for k = 1:nLabels
    m = masks(:,:,k);
    a = nnz(m);
    areaPx(k) = a;

    if a == 0
        info.dropped.empty = info.dropped.empty + 1;
        droppedIDs(end+1) = labelIDs(k); %#ok<AGROW>
        continue
    end
    if a < opts.minAreaPx
        info.dropped.minArea = info.dropped.minArea + 1;
        droppedIDs(end+1) = labelIDs(k); %#ok<AGROW>
        continue
    end
    if a > opts.maxAreaPx
        info.dropped.maxArea = info.dropped.maxArea + 1;
        droppedIDs(end+1) = labelIDs(k); %#ok<AGROW>
        continue
    end

    rws = find(any(m,2));
    cls = find(any(m,1));
    mg  = opts.edgeMarginPx;
    if rws(1) <= mg || rws(end) > H-mg || cls(1) <= mg || cls(end) > W-mg
        info.dropped.edge = info.dropped.edge + 1;
        droppedIDs(end+1) = labelIDs(k); %#ok<AGROW>
        continue
    end

    keep(k) = true;
    [r,c] = find(m);
    centroids(k,:) = [mean(c) mean(r)];
end

info.droppedIDs = droppedIDs(:)';
keepIdx = find(keep);
if isempty(keepIdx)
    moCorROI = emptyROIstruct();
    warning('labelImg2moCorROI:allDropped',...
        'All %d labels dropped by QC (minArea=%g maxArea=%g edgeMargin=%d).',...
        nLabels,opts.minAreaPx,opts.maxAreaPx,opts.edgeMarginPx);
    return
end

%--- deterministic ordering: centroid, row-major --------------------------
%the label image's own numbering is an implementation detail of the
%segmenter; sorting by position means the same image always gives the same
%ID assignment, which is what makes an ROI set comparable across re-runs.
[~,ord] = sortrows([round(centroids(keepIdx,2)) round(centroids(keepIdx,1))]);
keepIdx = keepIdx(ord);

%--- build the struct -----------------------------------------------------
moCorROI = emptyROIstruct();
for n = 1:numel(keepIdx)
    k = keepIdx(n);
    m = masks(:,:,k);

    B = bwboundaries(m,8,'noholes');
    if numel(B) > 1
        %a label that is not one connected blob: keep the largest piece so
        %the outline and the mask describe the same object
        [~,big] = max(cellfun(@(b) size(b,1),B));
        warning('labelImg2moCorROI:splitLabel',...
            'Label %d has %d connected components; keeping the largest outline.',...
            labelIDs(k),numel(B));
        B = B(big);
    end
    bnd = B{1};                       % [row col], closed (last == first)
    xy  = [bnd(:,2) bnd(:,1)];        % -> [x y]

    rws = find(any(m,2));
    cls = find(any(m,1));
    pos = [cls(1) rws(1) cls(end)-cls(1) rws(end)-rws(1)];   % imellipse-style

    ID = num2str(opts.startID + n - 1);

    moCorROI(n).ID                = ID;
    moCorROI(n).pos               = pos;
    moCorROI(n).XYvertices        = xy;
    moCorROI(n).frame             = 1;
    moCorROI(n).deleted           = false;
    moCorROI(n).type              = 'cellpose';
    moCorROI(n).mask              = m;
    moCorROI(n).ROIxyCoord        = mask2polyCoord(m);
    moCorROI(n).ROIcurveOrderedXY = closeContour(xy.');   % 2 x nPts, row1 = x
    moCorROI(n).label             = struct('String',ID,...
                                           'Position',[pos(1) pos(2)-5 0]);
end

info.nKept       = numel(moCorROI);
info.areaPx      = areaPx(keepIdx)';
info.equivDiamPx = 2*sqrt(info.areaPx/pi);
end

% ------------------------------------------------------------------------
function s = emptyROIstruct()
s = struct('ID',{},'pos',{},'XYvertices',{},'frame',{},'deleted',{},...
    'type',{},'mask',{},'ROIxyCoord',{},'ROIcurveOrderedXY',{},'label',{});
end
% ------------------------------------------------------------------------
function c = closeContour(c)
%same closing convention orderEllipsePtOnCurve uses, so FISSA sees the
%polygon format it already gets from the GUI path
if ~all(c(:,end) == c(:,1))
    c(:,end+1) = c(:,1);
end
end
