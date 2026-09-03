function testResolveROIset()
% testResolveROIset  Unit test for the membership-first ROI-set resolver.
% Guards against silently labelling a stim group's traces with the wrong
% condition's ROI ID order when two conditions share an ROI count.
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
testConfig();
mkset=@(tag,n) arrayfun(@(i) struct('ID',num2str(i),'cond',tag),1:n);

% CASE 1: ambiguous count (both 8) -> membership must disambiguate
roiSets={mkset('full256',8),mkset('crop128',8)}; roiCounts=[8 8];
roiTifs={{'A_00001.tif','A_00002.tif'},{'A_00041.tif','A_00042.tif'}};
grp=struct('name',{'A_00041.tif','A_00042.tif'},'moCorRawFroi',{zeros(8,40),zeros(8,40)});
mc=resolveROIset(grp,roiSets,roiTifs,roiCounts);
assert(strcmp(mc(1).cond,'crop128'),'count-first would wrongly pick full256');

% CASE 2: subset membership
mc2=resolveROIset(struct('name',{'A_00001.tif'},'moCorRawFroi',{zeros(8,40)}),roiSets,roiTifs,roiCounts);
assert(strcmp(mc2(1).cond,'full256'));

% CASE 3: legacy (no tif lists) -> count fallback
mc3=resolveROIset(struct('name','X.tif','moCorRawFroi',zeros(8,40)),{mkset('a',10),mkset('b',8)},{{},{}},[10 8]);
assert(strcmp(mc3(1).cond,'b'));

% CASE 4: realistic AA0072 (256^2=10, 128=8) spont -> crop
roiSetsR={mkset('full',10),mkset('crop',8)};
roiTifsR={{'AA_00047_00001.tif','AA_00048_00001.tif'},{'AA_00041_00001.tif','AA_00042_00001.tif'}};
spont=struct('name',{'AA_00041_00001.tif','AA_00042_00001.tif'},'moCorRawFroi',{zeros(8,40),zeros(8,40)});
mc4=resolveROIset(spont,roiSetsR,roiTifsR,[10 8]);
assert(strcmp(mc4(1).cond,'crop'));

disp('testResolveROIset PASS: membership-first prevents trace/ID misattribution');
end
