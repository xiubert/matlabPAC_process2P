function roiOut = remapROItoAcq(roiIn, srcHeader, tgtHeader, neuropilMarginPx)
% remapROItoAcq  Reuse 256x256 ROIs on a centered-crop 256x128 acquisition.
%
%   roiOut = remapROItoAcq(roiIn, srcHeader, tgtHeader)
%   roiOut = remapROItoAcq(roiIn, srcHeader, tgtHeader, neuropilMarginPx)
%
%   Maps an existing moCorROI struct array (drawn on a SOURCE acquisition,
%   e.g. 256x256 @ 5 Hz) onto a TARGET acquisition (e.g. 256x128 @ 10 Hz),
%   keeping only ROIs whose footprint lies fully within the target frame and
%   cropping their masks/polygons into target pixel coordinates.
%
%   This assumes the intended 10 Hz protocol: SAME zoom, full-width fast axis,
%   reduced slow-axis multiplier, no shift/rotation, and matched pixel pitch
%   -- i.e. a CENTERED VERTICAL CROP at 1:1 magnification. Under that
%   assumption the map is exact and integer:
%       offsetTop  = (srcRows - tgtRows)/2          % e.g. (256-128)/2 = 64
%       row_tgt = row_src - offsetTop ,  col_tgt = col_src
%   so a contained ROI's cropped mask is bit-identical to the source pixels.
%
%   If the two headers do NOT describe that clean crop (zoom differs, slow
%   axis not centered, pitch mismatch -> magnification, non-zero rotation),
%   the function ERRORS rather than silently sampling the wrong cells. This
%   is deliberate: an off-protocol acquisition (e.g. an erroneous zoom=2)
%   would otherwise yield correct-looking [tgtRows tgtCols] masks over the
%   wrong field of view.
%
%   Inputs
%     roiIn             moCorROI struct array (fields: ID, pos, XYvertices,
%                       frame, deleted, type, mask, ROIxyCoord,
%                       ROIcurveOrderedXY, label). Deleted ROIs are skipped.
%     srcHeader         ScanImage header for the source frame (from
%                       readSCIMtif(srcTif,'metaOnly')); needs hRoiManager.
%     tgtHeader         ScanImage header for the target frame.
%     neuropilMarginPx  optional; ROIs must sit at least this many pixels
%                       inside the crop on every side (so FISSA's neuropil
%                       annulus is not clipped at the frame edge). Default 0.
%
%   Output
%     roiOut    moCorROI struct array (same fields/format as roiIn) containing
%               only the contained ROIs, with mask, XYvertices, pos, label,
%               ROIxyCoord and ROIcurveOrderedXY all in target coordinates.
%               IDs are preserved so intersectROIfiles can match conditions.
%
%   See also remapROIfile, TIFcatROIgui, mask2polyCoord, orderEllipsePtOnCurve

arguments
    roiIn            struct
    srcHeader        struct
    tgtHeader        struct
    neuropilMarginPx (1,1) double {mustBeNonnegative} = 0
end

sg = getGeom(srcHeader, 'source');
tg = getGeom(tgtHeader, 'target');

% --- refuse anything that is not the expected clean centered crop ---------
validateCleanCrop(sg, tg);

offsetTop  = (sg.rows - tg.rows)/2;
offsetLeft = (sg.cols - tg.cols)/2;
bandRows = (offsetTop  + 1):(offsetTop  + tg.rows);
bandCols = (offsetLeft + 1):(offsetLeft + tg.cols);
mg = neuropilMarginPx;

% --- live ROIs only -------------------------------------------------------
if isfield(roiIn, 'deleted')
    live = roiIn(~[roiIn.deleted]);
else
    live = roiIn;
end

keep = false(1, numel(live));
out  = live;   % same struct fields; trimmed to `keep` at the end
for k = 1:numel(live)
    m = live(k).mask;
    rws = find(any(m, 2));
    cls = find(any(m, 1));
    if isempty(rws)
        continue
    end

    % fully contained (with optional neuropil margin)
    contained = rws(1)   >= bandRows(1)   + mg && ...
                rws(end) <= bandRows(end) - mg && ...
                cls(1)   >= bandCols(1)   + mg && ...
                cls(end) <= bandCols(end) - mg;
    if ~contained
        continue
    end
    keep(k) = true;

    r = live(k);
    r.mask = m(bandRows, bandCols);                 % bit-exact crop
    if isfield(r, 'XYvertices') && ~isempty(r.XYvertices)
        r.XYvertices = [r.XYvertices(:,1) - offsetLeft, r.XYvertices(:,2) - offsetTop];
    end
    if isfield(r, 'pos') && numel(r.pos) >= 2
        r.pos(1) = r.pos(1) - offsetLeft;
        r.pos(2) = r.pos(2) - offsetTop;
    end
    if isfield(r, 'label') && isstruct(r.label) && isfield(r.label, 'Position')
        r.label.Position(1) = r.label.Position(1) - offsetLeft;
        r.label.Position(2) = r.label.Position(2) - offsetTop;
    end
    % regenerate polygon coords from the cropped mask, exactly as
    % TIFcatROIgui does, so FISSA receives well-formed ROIcurveOrderedXY
    r.ROIxyCoord        = mask2polyCoord(r.mask);
    r.ROIcurveOrderedXY = orderEllipsePtOnCurve(r.ROIxyCoord);
    r.deleted = false;
    out(k) = r;
end

roiOut = out(keep);
end

% ------------------------------------------------------------------------
function validateCleanCrop(sg, tg)
% error unless (sg -> tg) is a centered, zoom-matched, unrotated, 1:1-pitch
% vertical crop.
tol = 1e-9;

if abs(sg.zoom - tg.zoom) > tol
    error('remapROItoAcq:zoomMismatch', ...
        ['Source and target scanZoomFactor differ (%g vs %g). The 10 Hz ' ...
         'protocol must use the same zoom as the source; a magnified frame ' ...
         'would sample different cells. Re-acquire at matching zoom.'], ...
        sg.zoom, tg.zoom);
end
if abs(sg.mFast - tg.mFast) > tol
    error('remapROItoAcq:fastAxisMismatch', ...
        ['scanAngleMultiplierFast differs (%g vs %g); the fast (X) axis ' ...
         'field of view must match for a pure vertical crop.'], sg.mFast, tg.mFast);
end
if any(abs([sg.sFast sg.sSlow tg.sFast tg.sSlow]) > tol)
    error('remapROItoAcq:offCenter', ...
        ['Non-zero scanAngleShift detected; remap assumes a centered crop ' ...
         '(all shifts = 0).']);
end
if any(abs([sg.rot tg.rot]) > tol)
    error('remapROItoAcq:rotated', ...
        'Non-zero scanRotation detected; remap assumes axis-aligned fields.');
end
% pixel pitch must match (no magnification): tgt/src pixel ratio == FOV ratio
if abs(tg.rows/sg.rows - tg.mSlow/sg.mSlow) > tol
    error('remapROItoAcq:pitchMismatchSlow', ...
        ['Slow-axis pixel pitch mismatch: linesPerFrame ratio %g != ' ...
         'multiplierSlow ratio %g (would imply vertical magnification).'], ...
        tg.rows/sg.rows, tg.mSlow/sg.mSlow);
end
if abs(tg.cols/sg.cols - tg.mFast/sg.mFast) > tol
    error('remapROItoAcq:pitchMismatchFast', ...
        ['Fast-axis pixel pitch mismatch: pixelsPerLine ratio %g != ' ...
         'multiplierFast ratio %g.'], tg.cols/sg.cols, tg.mFast/sg.mFast);
end
% centered crop must land on integer pixel offsets
if mod(sg.rows - tg.rows, 2) ~= 0 || (sg.rows - tg.rows) < 0
    error('remapROItoAcq:badRowCrop', ...
        'Cannot center-crop %d source rows to %d target rows.', sg.rows, tg.rows);
end
if mod(sg.cols - tg.cols, 2) ~= 0 || (sg.cols - tg.cols) < 0
    error('remapROItoAcq:badColCrop', ...
        'Cannot center-crop %d source cols to %d target cols.', sg.cols, tg.cols);
end
end

% ------------------------------------------------------------------------
function g = getGeom(h, which)
if ~isfield(h, 'hRoiManager')
    error('remapROItoAcq:noRoiManager', '%s header lacks hRoiManager.', which);
end
rm = h.hRoiManager;
g.cols  = getf(rm, 'pixelsPerLine');
g.rows  = getf(rm, 'linesPerFrame');
g.mFast = getf(rm, 'scanAngleMultiplierFast');
g.mSlow = getf(rm, 'scanAngleMultiplierSlow');
g.zoom  = getf(rm, 'scanZoomFactor');
g.sFast = getf(rm, 'scanAngleShiftFast', 0);
g.sSlow = getf(rm, 'scanAngleShiftSlow', 0);
g.rot   = getf(rm, 'scanRotation', 0);
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = double(s.(name));
elseif nargin == 3
    v = default;
else
    error('remapROItoAcq:missingField', 'hRoiManager.%s is required.', name);
end
end
