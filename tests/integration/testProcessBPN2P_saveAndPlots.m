% testProcessBPN2P_saveAndPlots
%
% Drive processBPN2P end-to-end against a COPY of the TO0003 _raw table, so
% the recorded animal folder is never written to.
%
% Verifies:
%   - script completes without error
%   - saved bundle contains the expected variables
%   - anmlROIbyStim has the expected derived columns
%   - dimensions in the saved table match the e2e expectation
%     (108 rows for TO0003, dFF [10 15] per row, etc.)
%   - the _raw input is left untouched (the two-stage convention)

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();

dp = requireFixture(cfg.exampleAnimalDir,'TO0003 animal folder');
animal = 'TO0003';
rawFile = fullfile(dp, [animal '_anmlROI_BPNstimTable_raw.mat']);
assert(isfile(rawFile), 'Missing _raw.mat; run stimParam2ROI first.');

% processBPN2P resolves its own input: it LOADS <animal>_..._raw.mat out of
% dataPath rather than reading an anmlROIbyStim left in the workspace. So the
% sandbox has to be a real animal folder with the _raw table copied into it --
% pre-setting the table variables here would be silently ignored.
[sandbox,sandboxCleanup] = testSandbox('processBPN2P_saveAndPlots'); %#ok<ASGLU>
sandboxDir = fullfile(sandbox,animal);
mkdir(sandboxDir);
sandboxRaw = fullfile(sandboxDir,[animal '_anmlROI_BPNstimTable_raw.mat']);
copyfile(rawFile, sandboxRaw);
rawBytes0 = dir(sandboxRaw).bytes;
dataPath = sandboxDir;

% Fingerprint the real animal folder's processed file so the run can be proven
% not to have touched it. This test is one clobbered variable away from
% overwriting the researcher's data, so the guard is worth the four lines.
realProc = fullfile(dp,[animal '_anmlROI_BPNstimTable.mat']);
realProcStamp0 = '';
if isfile(realProc); realProcStamp0 = dir(realProc).date; end

% NB: the _raw bundle stores a dataPath of its own, and processBPN2P loads it
% with a bare `load(rawFile)` -- so `dataPath` is clobbered back to the
% recording folder partway through the script. The save target is safe because
% procFile is resolved BEFORE that load, but every assertion below has to use
% sandboxDir, not dataPath, or it would inspect the real animal folder.

% Figures are created but never shown; the plotting blocks still have to run.
set(0, 'DefaultFigureVisible', 'off');

run(fullfile(cfg.repoRoot,'stimulusSpecific','processBPN2P.m'));

% Read back the saved file
outFile = fullfile(sandboxDir, [animal '_anmlROI_BPNstimTable.mat']);
assert(isfile(outFile), 'processBPN2P did not write %s', outFile);
Sout = load(outFile);

expectVars = {'anmlROIbyStim','stimTable','tifStimParamTable','dataPath'};
for v = expectVars
    assert(isfield(Sout, v{1}), 'Saved file missing variable %s', v{1});
end

T = Sout.anmlROIbyStim;
expectCols = {'animal','roiID','stimID','Pulse','BPNsOnset','BPNdBAmpl', ...
              'BPNmsPulseLen','BPNmsStimLen','BPNkHzFreqCutoff','treatment', ...
              'frameRate','rawFroi','moCorRawFroi','fissaFroi','SCALEDfissaFroi', ...
              't_total','t_dFF','dFF','dFF_avg','sigPeak','sig','pkResp'};
missing = setdiff(expectCols, T.Properties.VariableNames);
assert(isempty(missing), 'Saved anmlROIbyStim missing columns: %s', strjoin(missing,', '));

assert(height(T) == 108, 'Expected 108 rows in saved table, got %d', height(T));
sz = size(T.dFF{1});
assert(isequal(sz, [10 15]), 'Expected dFF{1} = [10 15], got %s', mat2str(sz));
sz = size(T.SCALEDfissaFroi{1});
assert(isequal(sz, [10 15]), 'Expected SCALEDfissaFroi{1} = [10 15], got %s', mat2str(sz));

% The _raw input is the one artifact processBPN2P must never touch.
assert(isfile(sandboxRaw), 'processBPN2P deleted its _raw input');
assert(dir(sandboxRaw).bytes == rawBytes0, 'processBPN2P rewrote its _raw input');

% ... and nothing was written into the recorded animal folder.
realProcStamp1 = '';
if isfile(realProc); realProcStamp1 = dir(realProc).date; end
assert(isequal(realProcStamp0,realProcStamp1), ...
    'processBPN2P wrote into the recorded animal folder: %s', realProc);

fprintf('Saved file OK: %s (%d rows, all expected columns present).\n', outFile, height(T));
fprintf('_raw input untouched.\n');
fprintf('All save-path assertions passed.\n');
