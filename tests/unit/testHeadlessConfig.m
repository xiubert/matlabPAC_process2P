function testHeadlessConfig()
% testHeadlessConfig  Unit test for the headless run's config resolution.
% Every dialog processAnimal2P shows becomes a config field, so the risk moves
% from "the user mis-clicked" to "the selector silently matched the wrong
% tifs". This pins each selector form and checks that a selector matching
% nothing is an ERROR rather than an empty group.
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
testConfig();

[tmp,cleanup] = testSandbox('headlessConfig'); %#ok<ASGLU>

% a synthetic animal folder: 6 tifs, the last two deliberately large so the
% byte heuristic has something to find
names = {'TO9999AAAA_00001_00001.tif','TO9999AAAA_00002_00001.tif',...
         'TO9999AAAA_00003_00001.tif','TO9999AAAA_00004_00001.tif',...
         'TO9999AAAA_00041_00001.tif','TO9999AAAA_00042_00001.tif'};
bytes = [1e3 1e3 1e3 1e3 20e6 20e6];
for k = 1:numel(names)
    fid = fopen(fullfile(tmp,names{k}),'w');
    fwrite(fid,zeros(bytes(k),1,'uint8'));
    fclose(fid);
end

% CASE 1: defaults -- animal from the path, map tifs by the byte heuristic
cfg = headlessConfig(tmp,'animal','TO9999','verbose',false);
assert(strcmp(cfg.animal,'TO9999'));
assert(cfg.nTifsFound==6,'expected 6 tifs, got %d',cfg.nTifsFound);
assert(nnz(cfg.mapIDX)==2,'byte heuristic should flag 2 map tifs, got %d',nnz(cfg.mapIDX));
assert(~any(cfg.preIDX),'no treatment set -> no pre tifs');
assert(isequal(cfg.stages,1:11));

% CASE 1b: runs are ISOLATED BY DEFAULT -- a derived run name, and artifacts
% under analysis/<run>/ rather than loose in the animal folder
assert(~isempty(cfg.run),'a run name should be derived when none is given');
assert(~cfg.paths.isFlat,'runs must be isolated by default');
assert(strcmp(cfg.paths.artifacts,fullfile(tmp,'analysis',cfg.run)));
assert(strcmp(cfg.runLabel,cfg.run),'the stamp must follow the folder name');
assert(contains(cfg.run,'savedROI') || contains(cfg.run,'cellpose'), ...
    'the derived name should say which ROI source it used, got %s',cfg.run);

% ...and the flat layout is still reachable for reproducing an old run in place
flat = headlessConfig(tmp,'animal','TO9999','flatLayout',true,'verbose',false);
assert(flat.paths.isFlat && strcmp(flat.paths.artifacts,tmp));

% the calibrated ROI defaults must be what the headless run actually uses
assert(strcmp(cfg.roi.mode,'consensus'));
assert(cfg.roi.minVotes==2 && cfg.roi.dilatePx==2,...
    'ROI defaults drifted from the TO0003 calibration');

% CASE 2: numeric indices
cfg = headlessConfig(tmp,'animal','TO9999','treatmentName','ZX1',...
    'preTifs',1:3,'verbose',false);
assert(isequal(find(cfg.preIDX)',1:3));

% CASE 3: explicit file names
cfg = headlessConfig(tmp,'animal','TO9999','treatmentName','ZX1',...
    'preTifs',names(1:2),'verbose',false);
assert(nnz(cfg.preIDX)==2 && all(cfg.preIDX(1:2)));

% CASE 4: a regular expression on the names
cfg = headlessConfig(tmp,'animal','TO9999','treatmentName','ZX1',...
    'preTifs','_0000[12]_','verbose',false);
assert(nnz(cfg.preIDX)==2 && all(cfg.preIDX(1:2)),...
    'regexp selector matched %s',mat2str(find(cfg.preIDX)'));

% CASE 5: a logical mask
mask = false(6,1); mask([2 4]) = true;
cfg = headlessConfig(tmp,'animal','TO9999','treatmentName','ZX1',...
    'preTifs',mask,'verbose',false);
assert(isequal(cfg.preIDX,mask));

% CASE 6: mapTifs off entirely
cfg = headlessConfig(tmp,'animal','TO9999','mapTifs',[],'verbose',false);
assert(~any(cfg.mapIDX));

% CASE 7: [from to] expands to a stage range
cfg = headlessConfig(tmp,'animal','TO9999','stages',[4 11],'verbose',false);
assert(isequal(cfg.stages,4:11),'stage range not expanded: %s',mat2str(cfg.stages));

% CASE 8: selectors that match nothing must ERROR. This is the whole point:
% a typo'd regexp that silently selects no tifs would relabel the entire
% session as post-treatment and nobody would see it.
assertThrows(@() headlessConfig(tmp,'animal','TO9999','preTifs',[1 99],...
    'treatmentName','ZX1','verbose',false),'headlessConfig:badSelector',...
    'out-of-range index');
assertThrows(@() headlessConfig(tmp,'animal','TO9999','mapTifs','NOSUCH',...
    'verbose',false),'headlessConfig:badSelector','unmatched regexp');
assertThrows(@() headlessConfig(tmp,'animal','TO9999','preTifs',{'nope.tif'},...
    'treatmentName','ZX1','verbose',false),'headlessConfig:badSelector',...
    'unknown file name');
assertThrows(@() headlessConfig(tmp,'animal','TO9999','stages',[0 5],...
    'verbose',false),'headlessConfig:badStages','out-of-range stage');

% CASE 9: a struct is accepted in place of name-value pairs, so a saved config
% round-trips
s = struct('animal','TO9999','treatmentName','ZX1','preTifs',1:3,'verbose',false);
cfgS = headlessConfig(tmp,s);
assert(isequal(find(cfgS.preIDX)',1:3),'struct form did not resolve selectors');

disp('testHeadlessConfig PASS: all selector forms resolve, bad ones error');
end

% ------------------------------------------------------------------------
function assertThrows(fh,id,what)
threw = false;
try
    fh();
catch ME
    threw = strcmp(ME.identifier,id);
    if ~threw
        error('expected %s for %s, got %s',id,what,ME.identifier);
    end
end
assert(threw,'%s should have raised %s',what,id);
end
