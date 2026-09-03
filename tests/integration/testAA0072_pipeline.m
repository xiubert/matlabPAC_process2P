function testAA0072_pipeline()
% testAA0072_pipeline  Real-data integration test of the 256x128 reuse path.
%
%   Dataset: /media/DATA/Ophys/Jinbo/Testdata_10Hz/AA0072 (mixed 256x256 @5Hz
%   + 256x128 @10Hz, no treatments, ROIs in a fluo2p roiOutput on tif 00047).
%
%   Exercises the automated path non-interactively, as far as is possible
%   without the interactive GUIs and the external Python FISSA step:
%     1. §2 frame-size auto-split  -> 256^2 condition + _128x256 condition
%     2. roiOutput(fluo2p) -> moCorROI source ROI set
%     3. §5b remapROItoAcq          -> contained ROIs cropped to 128 frame
%     4. §3 NoRMCorre on the 128 condition (size-driven; runs at 128 rows)
%     5. §7 moCorRawF2tifList        -> 10 Hz traces from the remapped ROIs
%
%   FISSA (§8) and the stimParam2ROI spont table (§10) require the external
%   FISSA run and are validated separately (resolveROIset has its own test).

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
dp = requireFixture(cfg.testdata10HzDir,'AA0072 10 Hz test dataset');
if exist('normcorre_batch','file')~=2
    addpath(genpath(requireFixture(cfg.normcorreDir,'NoRMCorre checkout')));
end
fprintf('=== testAA0072_pipeline ===\n');
assert(isfolder(dp),'dataset not found: %s',dp);

%% 1. inventory + split by full scan GEOMETRY (dims + zoom + multipliers)
% NB: pixel-dims alone are insufficient here -- tif 00040 is 128x256 like the
% other spont tifs but was acquired at zoom=2 (a different, magnified FOV), so
% it must NOT be grouped with the zoom=1 spont tifs. Splitting on the full
% geometry signature isolates it.
tifFiles = dir(fullfile(dp,'AA0072AAAA_*_00001.tif'));
geom = zeros(numel(tifFiles),5);   % [H W zoom mFast mSlow]
for i = 1:numel(tifFiles)
    [~,hi] = readSCIMtif(fullfile(tifFiles(i).folder,tifFiles(i).name),'metaOnly');
    rm = hi.hRoiManager;
    geom(i,:) = [hi.imHeight hi.imWidth rm.scanZoomFactor ...
                 rm.scanAngleMultiplierFast rm.scanAngleMultiplierSlow];
end
[uG,~,ig] = unique(geom,'rows','stable');
fprintf('1. geometry groups:\n');
for k = 1:size(uG,1)
    fprintf('   %dx%d zoom=%g mFast=%g mSlow=%g : %d tifs\n', ...
        uG(k,2),uG(k,1),uG(k,3),uG(k,4),uG(k,5), sum(ig==k));
end
src256IDX  = geom(:,1)==256 & geom(:,3)==1;                 % 256^2 source
spont128IDX = geom(:,1)==128 & geom(:,3)==1 & geom(:,5)==0.5; % valid 10Hz crop
stray128IDX = geom(:,1)==128 & geom(:,3)~=1;                 % misconfigured crop
fprintf('   -> 256^2 source: %d | valid 128 spont: %d | stray (zoom~=1) 128: %d\n', ...
    sum(src256IDX), sum(spont128IDX), sum(stray128IDX));
assert(sum(spont128IDX)>=2,'need >=2 correctly-configured 128 tifs');
if any(stray128IDX)
    fprintf('   WARNING: excluding misconfigured 128 tif(s): %s\n', ...
        strjoin({tifFiles(stray128IDX).name},', '));
end
src256   = tifFiles(src256IDX);
files128 = tifFiles(spont128IDX);

%% 2. roiOutput(fluo2p) -> moCorROI-format source
S = load(fullfile(dp,'AA0072AAAA_00047_00001_roiOutput.mat'),'fluo2p');
src = S.fluo2p.roi;
for k = 1:numel(src)
    src(k).ROIxyCoord = mask2polyCoord(src(k).mask);
    src(k).ROIcurveOrderedXY = orderEllipsePtOnCurve(src(k).ROIxyCoord);
end
fprintf('2. source ROIs from roiOutput: %d (mask %s)\n', numel(src), mat2str(size(src(1).mask)));

%% 3. §5b remap 256^2 -> 128
[~,srcH] = readSCIMtif(fullfile(src256(1).folder,src256(1).name),'metaOnly');
[~,tgtH] = readSCIMtif(fullfile(files128(1).folder,files128(1).name),'metaOnly');
tgtROI = remapROItoAcq(src, srcH, tgtH);
fprintf('3. remap: %d of %d ROIs contained in 128 crop; mask %s\n', ...
    numel(tgtROI), numel(src), mat2str(size(tgtROI(1).mask)));
assert(all(arrayfun(@(r) isequal(size(r.mask),[128 256]), tgtROI)),'remapped masks not 128x256');
assert(~isempty(tgtROI),'no ROIs survived the crop');

%% 4. §3 NoRMCorre on the valid 128 condition
fprintf('4. motion-correcting %d x 128 tifs ...\n', numel(files128));
[Ycon,~] = concatenate_files(files128);
raw = single(Ycon); raw = raw - min(raw(:));
nF = size(raw,3); ib = min(200,nF);
opts = NoRMCorreSetParms('d1',size(raw,1),'d2',size(raw,2),'grid_size',[32,32],...
    'mot_uf',4,'bin_width',ib,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',ib);
M = normcorre_batch(raw,opts);
fprintf('   corrected stack %s  OK\n', mat2str(size(M)));
assert(size(M,1)==128 && size(M,2)==256,'corrected stack not 128x256');

%% 5. §7 extract 10 Hz traces from the remapped ROIs
files128 = moCorRawF2tifList(files128, M, tgtROI, raw);
F = files128(1).moCorRawFroi;
fprintf('5. extracted traces: %s, frameRate=%.4f Hz\n', mat2str(size(F)), files128(1).frameRate);
assert(size(F,1)==numel(tgtROI),'trace row count != #remapped ROIs');
assert(all(isfinite(F(:))),'non-finite values in extracted traces');
assert(round(files128(1).frameRate)==10,'expected ~10 Hz');

fprintf('=== PASS: size-split -> remap -> motion-correct -> 10 Hz extraction ===\n');
fprintf('NOTE: FISSA (§8) + stimParam2ROI spont table (§10) require the external FISSA run.\n');
end
