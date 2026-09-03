function testAnimalPaths()
% testAnimalPaths  Unit test for the artifact layout contract.
% Pins which artifacts are SHARED across analyses of an animal and which
% belong to one run, and that a named run never writes where another run
% would -- the collision the layout exists to prevent.
addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
testConfig();

[tmp,cleanup] = testSandbox('animalPaths'); %#ok<ASGLU>
dp = fullfile(tmp,'TO9999');
mkdir(fullfile(dp,'NoRMCorred'));

% CASE 1: no run -> the flat legacy layout, unchanged from before
P = animalPaths(dp);
assert(P.isFlat,'no run should resolve flat');
assert(strcmp(P.artifacts,dp));
assert(strcmp(P.fissaDir,fullfile(dp,'NoRMCorred','FISSAoutput')), ...
    'a legacy folder must still find FISSA output under NoRMCorred/');
assert(strcmp(P.tifFileList,fullfile(dp,'TO9999_tifFileList.mat')));

% CASE 2: a named run puts per-run artifacts under analysis/<run>/
R = animalPaths(dp,'run','cellpose_20260903');
assert(~R.isFlat);
assert(strcmp(R.artifacts,fullfile(dp,'analysis','cellpose_20260903')));
assert(strcmp(R.tifFileList,fullfile(R.artifacts,'TO9999_tifFileList.mat')));

% CASE 3: SHARED artifacts do not move. The tif inventory, the condition
% split and NoRMCorred/ describe the acquisition and the motion correction,
% are identical for every analysis, and are expensive -- duplicating them per
% run would mean re-running NoRMCorre for each.
assert(strcmp(R.legend,P.legend),   'tif legend must stay shared');
assert(strcmp(R.condSplit,P.condSplit),'condition split must stay shared');
assert(strcmp(R.moCorrDir,P.moCorrDir),'NoRMCorred must stay shared');
assert(strcmp(R.ncParams,P.ncParams),  'NoRMCorre params must stay shared');

% CASE 4: a named run OWNS its FISSA output. Falling back to the shared
% NoRMCorred/FISSAoutput when the per-run one does not exist yet would send
% every run's output to the same place -- the exact collision being prevented.
mkdir(fullfile(dp,'NoRMCorred','FISSAoutput'));
R2 = animalPaths(dp,'run','cellpose_20260903');
assert(strcmp(R2.fissaDir,fullfile(R2.artifacts,'FISSAoutput')), ...
    'a named run must not inherit the shared FISSAoutput');

% CASE 5: two runs collide nowhere
A = animalPaths(dp,'run','handdrawn_20260903');
B = animalPaths(dp,'run','cellpose_20260903');
for f = {'artifacts','tifFileList','moCorrTifs','fissaDir'}
    assert(~strcmp(A.(f{1}),B.(f{1})), ...
        'two runs share %s -- they would overwrite each other',f{1});
end

% CASE 7: 'create' makes the folders, and not creating is the default
assert(~isfolder(R.artifacts),'animalPaths must not create folders by default');
C = animalPaths(dp,'run','made','create',true);
assert(isfolder(C.artifacts) && isfolder(C.fissaDir));

% CASE 8: a run name that could escape the folder or break a glob is refused
for bad = {'../escape','a/b','a:b','a*b'}
    threw = false;
    try
        animalPaths(dp,'run',bad{1});
    catch ME
        threw = strcmp(ME.identifier,'animalPaths:badRunName');
    end
    assert(threw,'run name "%s" should have been refused',bad{1});
end

disp('testAnimalPaths PASS: shared stays shared, runs cannot collide');
end
