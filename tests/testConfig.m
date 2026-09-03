function cfg = testConfig()
%testConfig  Resolve repo and fixture paths for the matlabPAC_process2P tests.
%
%   cfg = testConfig()
%
%   Every test starts with this. It puts the repo on the MATLAB path and
%   returns the locations of the recorded datasets the integration tests read.
%   Nothing here asserts that a fixture EXISTS -- use requireFixture for that,
%   so a machine without the data drive reports SKIP instead of FAIL.
%
%   Fixture roots are overridable by environment variable, because the data
%   lives outside the repo on a lab drive and its mount point is not the same
%   on every machine:
%
%     MATLABPAC_TESTDATA   root holding the Ophys datasets
%                          (default /media/DATA/Ophys/Jinbo)
%     MATLABPAC_NORMCORRE  NoRMCorre checkout, needed only by the motion
%                          correction test
%                          (default ~/Documents/MATLAB/NoRMCorre)
%
%   Output (struct):
%     .repoRoot          repo top level
%     .testsRoot         this folder
%     .dataRoot          fixture root (see MATLABPAC_TESTDATA)
%     .aggregateDir      TOMT cohort group files (<Family>_Group<g>.mat)
%     .animalsRoot       TOMT per-animal folders (TO0001 ... TO0014)
%     .exampleAnimalDir  the TO0003 working copy CLAUDE.md points tests at
%     .testdata10HzDir   AA0072, the mixed 5 Hz / 10 Hz acquisition
%     .normcorreDir      NoRMCorre checkout
%
%   NOTE: .exampleAnimalDir and fullfile(.animalsRoot,'TO0003') are two
%   near-identical copies of the same animal on disk. The BPN and remap tests
%   were written against the first, the FRA tests against the second; both are
%   kept so neither set silently starts reading the other's artifacts.
%
%   See also requireFixture, testSandbox, runTests.

testsRoot = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(testsRoot);

cfg.repoRoot  = repoRoot;
cfg.testsRoot = testsRoot;

dataRoot = getenv('MATLABPAC_TESTDATA');
if isempty(dataRoot); dataRoot = fullfile('/media','DATA','Ophys','Jinbo'); end
cfg.dataRoot = dataRoot;

cfg.aggregateDir     = fullfile(dataRoot,'TOMT','aggregate data');
cfg.animalsRoot      = fullfile(dataRoot,'TOMT','animals');
cfg.exampleAnimalDir = fullfile(dataRoot,'TO0003');
cfg.testdata10HzDir  = fullfile(dataRoot,'Testdata_10Hz','AA0072');

normcorre = getenv('MATLABPAC_NORMCORRE');
if isempty(normcorre)
    normcorre = fullfile(getenv('HOME'),'Documents','MATLAB','NoRMCorre');
end
cfg.normcorreDir = normcorre;

% Put the code under test on the path. genpath for helperFcns because it is
% several levels deep; stimulusSpecific holds the process*/plot* scripts.
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot,'helperFcns')));
addpath(fullfile(repoRoot,'stimulusSpecific'));
addpath(testsRoot);
end
