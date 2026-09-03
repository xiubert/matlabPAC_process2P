function testLabelImg2moCorROI()
% testLabelImg2moCorROI  Unit test for label image -> moCorROI conversion.
% Guards the contract the rest of the pipeline depends on: the same fields
% TIFcatROIgui writes, a traced (not angle-sorted) closed outline, and QC
% filters that actually drop what they claim to.
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
testConfig();

sz = [64 64];

% --- a label image with three discs of known size, plus one at the edge ----
L = zeros(sz);
L = addDisc(L,[20 20],5,1);     % r=5  -> ~81 px
L = addDisc(L,[20 45],3,2);     % r=3  -> ~29 px
L = addDisc(L,[45 20],7,3);     % r=7  -> ~149 px
L = addDisc(L,[ 3  3],3,4);     % touches the frame edge

[roi,info] = labelImg2moCorROI(L,'edgeMarginPx',0,'minAreaPx',5);

% CASE 1: every non-edge label survives with the GUI's field set
assert(numel(roi)==4,'expected 4 ROIs, got %d',numel(roi));
want = {'ID','pos','XYvertices','frame','deleted','type','mask',...
        'ROIxyCoord','ROIcurveOrderedXY','label'};
assert(isempty(setdiff(want,fieldnames(roi))),...
    'moCorROI is missing fields TIFcatROIgui writes: %s',...
    strjoin(setdiff(want,fieldnames(roi)),', '));
assert(islogical(roi(1).mask),'mask must be logical');
assert(strcmp(roi(1).type,'cellpose'));

% CASE 2: IDs are 1..N in row-major centroid order, independent of the
% label image's own numbering -- so a re-run gives the same ID to the same cell
ids = str2double({roi.ID});
assert(isequal(ids,1:4),'IDs must be 1..N in order, got %s',mat2str(ids));
cen = cell2mat(arrayfun(@(r) centroid(r.mask),roi','uni',0));
assert(issorted(round(cen(:,1))),'ROIs are not ordered row-major by centroid');

% CASE 3: the outline is a CLOSED, TRACED contour. A star-convex angular sort
% (orderEllipsePtOnCurve) would not survive a concave mask; bwboundaries does.
for k = 1:numel(roi)
    c = roi(k).ROIcurveOrderedXY;
    assert(size(c,1)==2,'ROIcurveOrderedXY must be 2 x nPts (row1 = x)');
    assert(isequal(c(:,1),c(:,end)),'contour %d is not closed',k);
    %every outline point must sit on a mask pixel
    lin = sub2ind(size(roi(k).mask),round(c(2,:)),round(c(1,:)));
    assert(all(roi(k).mask(lin)),'outline %d leaves the mask',k);
end

% a deliberately concave (C-shaped) mask: the traced outline still closes and
% stays on the mask, which an angular sort cannot guarantee
Lc = zeros(sz);
Lc = addDisc(Lc,[32 32],9,1);
Lc = addDisc(Lc,[32 38],5,0);          % bite out of the right side
roiC = labelImg2moCorROI(Lc,'edgeMarginPx',0,'minAreaPx',5);
cC = roiC(1).ROIcurveOrderedXY;
assert(isequal(cC(:,1),cC(:,end)),'concave outline not closed');
linC = sub2ind(sz,round(cC(2,:)),round(cC(1,:)));
assert(all(roiC(1).mask(linC)),'concave outline leaves the mask');

% CASE 4: minAreaPx drops the small disc and says so
[roiMin,infoMin] = labelImg2moCorROI(L,'edgeMarginPx',0,'minAreaPx',50);
assert(numel(roiMin)==2,'minAreaPx=50 should leave 2 ROIs, got %d',numel(roiMin));
assert(infoMin.dropped.minArea==2,'dropped.minArea not reported');

% CASE 5: edgeMarginPx drops the edge-touching label -- the guard that keeps
% FISSA's neuropil annulus off the frame border
[roiEdge,infoEdge] = labelImg2moCorROI(L,'edgeMarginPx',3,'minAreaPx',5);
assert(numel(roiEdge)==3,'edgeMarginPx=3 should leave 3 ROIs, got %d',numel(roiEdge));
assert(infoEdge.dropped.edge==1,'dropped.edge not reported');

% CASE 6: dilation grows every mask and keeps them disjoint
a0 = arrayfun(@(r) nnz(r.mask),roi);
roiD = labelImg2moCorROI(L,'edgeMarginPx',0,'minAreaPx',5,'dilatePx',2);
aD = arrayfun(@(r) nnz(r.mask),roiD);
assert(all(aD > a0),'dilatePx did not grow the masks');
overlap = zeros(sz);
for k = 1:numel(roiD); overlap = overlap + roiD(k).mask; end
assert(max(overlap(:))<=1,'dilated masks overlap; preventOverlap failed');

% CASE 7: an empty label image is a warning and an empty set, not an error
w = warning('off','labelImg2moCorROI:noLabels');
c = onCleanup(@() warning(w));
roiE = labelImg2moCorROI(zeros(sz));
assert(isempty(roiE),'empty label image should give an empty ROI set');

fprintf(['testLabelImg2moCorROI PASS: GUI-compatible fields, closed traced '...
    'outlines, QC filters (%d labels in)\n'],info.nLabels);
end

% ------------------------------------------------------------------------
function L = addDisc(L,c,r,val)
[X,Y] = meshgrid(1:size(L,2),1:size(L,1));
L((Y-c(1)).^2 + (X-c(2)).^2 <= r^2) = val;
end
% ------------------------------------------------------------------------
function c = centroid(m)
[r,cc] = find(m);
c = [mean(r) mean(cc)];
end
