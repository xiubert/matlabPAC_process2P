function testRemapROItoAcq_centeredCrop()
% testRemapROItoAcq_centeredCrop  Spec + fixture test for remapping 256x256
%   ROIs onto a centered-crop 256x128 (10 Hz) acquisition.
%
%   WHY THIS TEST EXISTS
%   --------------------
%   The intended 10 Hz protocol is acquired with scanZoomFactor = 1 (same as
%   the 256x256 reference), scanAngleMultiplierSlow = 0.5, and
%   linesPerFrame = 128. With shift = 0 (centered) this is a pure CENTERED
%   VERTICAL CROP at 1:1 magnification:
%       - X (fast axis): full width, col_128 = col_256
%       - Y (slow axis): keep central 128 rows -> drop top 64 / bottom 64
%                        row_128 = row_256 - 64        (offsetTop = 64)
%   No magnification, no re-rasterization: a 256^2 ROI's pixels either lie
%   wholly inside the central band [65,192] (kept) or not (dropped).
%
%   NOTE: the on-disk AA0071 test file was acquired at zoom = 2 IN ERROR
%   (a magnified sub-field). This test does NOT use it. Instead it FABRICATES
%   the correct zoom = 1 target by cropping the central 128 rows of a real
%   256^2 TO0003 stack -- which is bit-for-bit what a zoom=1 / multSlow=0.5 /
%   128-line acquisition would have produced. That gives an exact ground
%   truth, so the test runs today without waiting for new acquisition.
%
%   WHAT IT CHECKS
%   --------------
%   Part A (oracle, no helper required -- always runs):
%       For every non-deleted ROI fully inside the central band, the mean F
%       over the CROPPED mask on the CROPPED stack equals, frame-for-frame
%       and BIT-EXACTLY, the mean F over the ORIGINAL mask on the FULL stack
%       (they are literally the same pixels). This validates the fixture and
%       the correctness oracle the helper must satisfy.
%
%   Part B (contract for remapROItoAcq -- runs once the helper exists):
%       Proposed signature:
%           moCorROIout = remapROItoAcq(moCorROIin, srcHeader, tgtHeader)
%       Assertions:
%           1. kept IDs == the contained-ROI IDs computed in Part A
%           2. every returned mask is [128 256] logical, non-empty
%           3. bit-exact F equivalence using the helper's masks
%           4. idempotency: src==tgt header returns all ROIs, masks unchanged
%           5. zoom mismatch (tgt zoom~=src zoom) raises a warning
%
%   Run:  matlab -batch "addpath('tests'); runTests('Filter','testRemapROItoAcq_centeredCrop')"

%% ---- config -------------------------------------------------------------
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();

dataDir = requireFixture(cfg.exampleAnimalDir,'TO0003 animal folder');
roiFile = fullfile(dataDir,'TO0003_moCorrROI_all.mat');
chanID  = 2;                 % green channel, as in the pipeline

% intended-protocol geometry (zoom=1, multSlow=0.5, centered)
srcRows  = 256;
srcCols  = 256;
tgtRows  = 128;
offsetTop = (srcRows - tgtRows)/2;   % = 64 ; central band is rows 65:192
bandRows = (offsetTop+1):(offsetTop+tgtRows);

fprintf('=== testRemapROItoAcq_centeredCrop ===\n');

%% ---- data guards --------------------------------------------------------
if exist(roiFile,'file')~=2
    fprintf('SKIP: ROI file not found: %s\n', roiFile);
    return
end
tifDir = dir(fullfile(dataDir,'TO0003AAAA_*_00001.tif'));
if isempty(tifDir)
    fprintf('SKIP: no TO0003 tif found in %s\n', dataDir);
    return
end
tifPath = fullfile(tifDir(1).folder, tifDir(1).name);

%% ---- load real ROIs + a real 256^2 stack --------------------------------
S = load(roiFile,'moCorROI');
moCorROI = S.moCorROI;
fullStack = readSCIMtif(tifPath, chanID);          % [256 256 nFrames]
assert(isequal(size(fullStack,1,2),[srcRows srcCols]), ...
    'Expected a %dx%d source stack; got %s', srcRows, srcCols, mat2str(size(fullStack)));
nFrames = size(fullStack,3);
fprintf('Loaded %d ROIs and stack %s\n', numel(moCorROI), mat2str(size(fullStack)));

%% ---- fabricate the EXPECTED zoom=1 256x128 target -----------------------
croppedStack = fullStack(bandRows, :, :);          % [128 256 nFrames]
assert(isequal(size(croppedStack,1,2),[tgtRows srcCols]));

tgtHeader = makeHeader(srcCols, tgtRows, 128, 0.5, 1);   % multSlow=0.5, zoom=1
srcHeader = makeHeader(srcCols, srcRows, 256, 1.0, 1);   % reference 256^2

% optionally persist a lightweight fixture for inspection (not the full stack)
% save(fullfile(cfg.testsRoot,'fixture_expected128.mat'),'tgtHeader','srcHeader', ...
%      'offsetTop','bandRows','-v7');

%% ---- determine the expected contained set (ground truth) ----------------
% A 256^2 ROI is contained iff all its mask pixels lie within the band.
isLive = ~[moCorROI.deleted];
contained = false(1,numel(moCorROI));
for k = 1:numel(moCorROI)
    if ~isLive(k); continue; end
    rowsWithPixels = find(any(moCorROI(k).mask,2));
    contained(k) = ~isempty(rowsWithPixels) && ...
        rowsWithPixels(1) >= bandRows(1) && rowsWithPixels(end) <= bandRows(end);
end
keptIDs    = string({moCorROI(contained).ID});
droppedIDs = string({moCorROI(isLive & ~contained).ID});
fprintf('Contained ROIs (%d): %s\n', numel(keptIDs), strjoin(keptIDs,', '));
fprintf('Dropped ROIs   (%d): %s\n', numel(droppedIDs), strjoin(droppedIDs,', '));
assert(~isempty(keptIDs), 'No ROIs fall within the central band -- cannot validate.');

%% ===== PART A: oracle (runs today, no helper) ============================
fprintf('\n-- Part A: ground-truth crop equivalence --\n');
nChecked = 0;
for k = find(contained)
    refMask = moCorROI(k).mask;            % 256x256
    cropMask = refMask(bandRows, :);       % 128x256  (== row_256-64 shift)

    Ffull = roiTrace(fullStack,  refMask);
    Fcrop = roiTrace(croppedStack, cropMask);

    assert(isequaln(Ffull, Fcrop), ...
        'ROI %s: cropped trace not bit-exact to full trace (max|d|=%g)', ...
        moCorROI(k).ID, max(abs(Ffull - Fcrop)));
    nChecked = nChecked + 1;
end
fprintf('PASS: %d contained ROIs give bit-exact F before/after centered crop.\n', nChecked);

%% ===== PART B: remapROItoAcq contract ====================================
fprintf('\n-- Part B: remapROItoAcq helper --\n');
if exist('remapROItoAcq','file')~=2
    fprintf(['PENDING: remapROItoAcq not on path. Part A defines its oracle; ' ...
             'implement the helper, then this section will validate it.\n']);
    fprintf('=== done (Part A passed; Part B pending) ===\n');
    return
end

out = remapROItoAcq(moCorROI, srcHeader, tgtHeader);

% 1. kept IDs match
outIDs = string({out.ID});
assert(isequal(sort(outIDs), sort(keptIDs)), ...
    'Helper kept IDs %s differ from expected %s', ...
    strjoin(sort(outIDs),','), strjoin(sort(keptIDs),','));

% 2. mask sizes + non-empty
for k = 1:numel(out)
    assert(islogical(out(k).mask) && isequal(size(out(k).mask),[tgtRows srcCols]), ...
        'ROI %s: mask is not %dx%d logical', out(k).ID, tgtRows, srcCols);
    assert(any(out(k).mask(:)), 'ROI %s: empty mask after remap', out(k).ID);
end

% 3. bit-exact F equivalence via helper masks
for k = 1:numel(out)
    src = moCorROI(strcmp({moCorROI.ID}, out(k).ID));
    Ffull = roiTrace(fullStack, src.mask);
    Fcrop = roiTrace(croppedStack, out(k).mask);
    assert(isequaln(Ffull, Fcrop), ...
        'ROI %s: helper mask F not bit-exact to source F (max|d|=%g)', ...
        out(k).ID, max(abs(Ffull - Fcrop)));
end

% 4. idempotency: src==tgt -> all ROIs, masks unchanged
same = remapROItoAcq(moCorROI, srcHeader, srcHeader);
assert(numel(same)==sum(isLive), 'Idempotent remap dropped ROIs');
for k = 1:numel(same)
    src = moCorROI(strcmp({moCorROI.ID}, same(k).ID));
    assert(isequal(same(k).mask, src.mask), 'ROI %s: mask changed under identity remap', same(k).ID);
end

% 5. zoom mismatch ERRORS (refuses to silently mis-sample a magnified frame)
badHeader = tgtHeader; badHeader.hRoiManager.scanZoomFactor = 2;
errId = '';
try
    remapROItoAcq(moCorROI, srcHeader, badHeader);
catch ME
    errId = ME.identifier;
end
assert(~isempty(errId), 'Expected an error when source/target zoom differ');

fprintf('PASS: remapROItoAcq satisfies the contract.\n');
fprintf('=== all checks passed ===\n');
end

%% ------------------------------------------------------------------------
function F = roiTrace(stack, mask)
% mean pixel intensity within mask, per frame (vectorized; matches
% moCorRawF2tifList's mean(f(mask)) exactly).
sz = size(stack);
flat = reshape(stack, sz(1)*sz(2), sz(3));
F = mean(flat(mask(:), :), 1);
end

function h = makeHeader(pixelsPerLine, linesPerFrame, lpf, multSlow, zoom)
% minimal ScanImage-like header carrying the fields the remap geometry needs.
h.imWidth  = pixelsPerLine;
h.imHeight = linesPerFrame;
h.nFrames  = NaN;
h.hRoiManager.pixelsPerLine            = pixelsPerLine;
h.hRoiManager.linesPerFrame            = lpf;
h.hRoiManager.scanAngleMultiplierFast  = 1;
h.hRoiManager.scanAngleMultiplierSlow  = multSlow;
h.hRoiManager.scanAngleShiftFast       = 0;
h.hRoiManager.scanAngleShiftSlow       = 0;
h.hRoiManager.scanZoomFactor           = zoom;
h.hRoiManager.scanRotation             = 0;
end
