function testFISSAgrouping()
% testFISSAgrouping  Self-contained regression test for the per-ROI-count
%   FISSA merge (mergeFISSAgroups) and the processAnimal2P section-9 grouped
%   branch. Builds synthetic per-group FISSA outputs + a manifest in the
%   layout the grouped FISSA driver (FISSAviaMatlab_prePostTreatment.py)
%   produces, so it needs neither FISSA nor Python.
%
%   Validates that a session mixing 256x256 (full ROI set) and 256x128
%   (reduced set) yields a tifFileList whose per-tif SCALEDfissaFroi row count
%   follows each tif's own group -- the contract stimParam2ROI's resolveROIset
%   relies on.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
testConfig();
fprintf('=== testFISSAgrouping ===\n');

nFr = 40; scale = 0.8;
[fdir,fdirCleanup] = testSandbox('FISSAgrouping'); %#ok<ASGLU>
mkGroup(fullfile(fdir,'g0'), 10, 2, nFr);   % 256^2 group: 10 ROIs, 2 tifs
mkGroup(fullfile(fdir,'g1'),  8, 2, nFr);   % 128   group:  8 ROIs, 2 tifs

manifest = [ ...
    struct('matfile','g0/matlab.mat','nROI',10, ...
           'tifNames',{{'A_00047_00001_NoRMCorre.tif','A_00048_00001_NoRMCorre.tif'}}); ...
    struct('matfile','g1/matlab.mat','nROI',8, ...
           'tifNames',{{'A_00041_00001_NoRMCorre.tif','A_00042_00001_NoRMCorre.tif'}}) ];
fid = fopen(fullfile(fdir,'groups.json'),'w'); fwrite(fid, jsonencode(manifest)); fclose(fid);

% synthesize tifList as section 7 leaves it (per-condition moCorRawFroi)
mk = @(nm,nroi,fr) struct('name',nm,'treatment','none','moCorRawFroi',rand(nroi,nFr),'frameRate',fr);
tifList.full(1) = mk('A_00047_00001.tif',10,5);
tifList.full(2) = mk('A_00048_00001.tif',10,5);
tifList.crop(1) = mk('A_00041_00001.tif',8,9.99);
tifList.crop(2) = mk('A_00042_00001.tif',8,9.99);

% ---- replicate the section-9 grouped branch ----
C = struct2cell(tifList); allcond = vertcat(C{:});
isMap = contains({allcond.treatment}','map');
tifFileList = struct();
tifFileList.stim = allcond(~isMap);
if any(isMap); tifFileList.map = allcond(isMap); end
tifFileList = mergeFISSAgroups(tifFileList, fdir, scale);

% ---- assertions ----
nm = {tifFileList.stim.name};
r256 = find(contains(nm,{'00047','00048'}));
r128 = find(contains(nm,{'00041','00042'}));
assert(numel(tifFileList.stim)==4,'expected 4 stim tifs');
for n = r256
    assert(isequal(size(tifFileList.stim(n).SCALEDfissaFroi),[10 nFr]),'256^2 SCALED size');
    assert(isequal(size(tifFileList.stim(n).fissaFroi),[10 nFr]),'256^2 fissaFroi size');
end
for n = r128
    assert(isequal(size(tifFileList.stim(n).SCALEDfissaFroi),[8 nFr]),'128 SCALED size');
    assert(tifFileList.stim(n).frameRate>9,'128 should be ~10Hz');
end
assert(all(cellfun(@(x) all(isfinite(x(:))),{tifFileList.stim.SCALEDfissaFroi})),'finite SCALED');

% resolveROIset contract: row count selects the matching ROI set
roiSets = {repmat(struct('ID','x'),10,1), repmat(struct('ID','y'),8,1)};
roiCounts = cellfun(@numel,roiSets);
for n = 1:numel(tifFileList.stim)
    c = size(tifFileList.stim(n).SCALEDfissaFroi,1);
    assert(roiCounts(find(roiCounts==c,1))==c,'resolveROIset would match row count');
end

fprintf('PASS: grouped merge -> 256^2 tifs=10 ROIs@5Hz, 128 tifs=8 ROIs@10Hz; row-count contract holds\n');
fprintf('=== done ===\n');
end

% ------------------------------------------------------------------------
function mkGroup(gdir, nCell, nTrial, nFr)
% write a synthetic FISSA save_to_matlab output: result.cell{c}.trial{t}
% = (nNeuropil+1) x nFrames, row 1 = ROI signal.
mkdir(gdir);
result = struct();
for c = 0:nCell-1
    cs = struct();
    for t = 0:nTrial-1
        cs.(sprintf('trial%d',t)) = rand(3,nFr);   % 1 ROI row + 2 neuropil rows
    end
    result.(sprintf('cell%d',c)) = cs;
end
save(fullfile(gdir,'matlab.mat'),'result'); %#ok<NASGU>
end
