% traceBPNcombineDiffOnset   [INVESTIGATION SCRIPT -- no assertions]
%
% Not part of the test suite: this prints findings, it does not check them.
% Kept because it is the record of how the combineDiffOnset stacking question
% was answered. Run it by hand when revisiting that behaviour.
%
% Trace the BPN analysis pipeline on TO0003 example data, starting at
% stimParams2TifTable (the same call processAnimal2P makes), then
% anmlROIbyStimTable, then combineDiffOnset.
%
% Purpose:
%   - Inspect the raw anmlROIbyStim table as it appears BEFORE
%     processBPN2P.m adds t_total/t_dFF/dFF columns.
%   - Verify what combineDiffOnset.m does on this table:
%       * does it find the *sOnset column (case-sensitivity check)?
%       * does the vertcat of SCALEDfissaFroi cells produce a matrix or
%         a cell array?
%       * does anything actually merge?
%   - Produce a small set of facts to inform a revision plan for
%     combineDiffOnset.m + processBPN2P.m.
%
% Inputs (all already on disk under /media/DATA/Ophys/Jinbo/TO0003):
%   TO0003_tifFileList.mat      (loads tifFileList; .stim has SCALEDfissaFroi)
%   TO0003_stimGroupIDX.mat     (loads stimGroupIDX; .BPNStimIDX.tifFileList)
%   TO0003_moCorrROI_all.mat    (loads moCorROI)
%   *_Pulses.mat alongside the tifs (consumed by stimParams2TifTable)
%
% Side effect: saves the rebuilt raw table to
%   <tempdir>/TO0003_raw_anmlROIbyStim.mat
% for follow-up multi-onset synthesis tests.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();
dp = requireFixture(cfg.exampleAnimalDir,'TO0003 animal folder');
animal = 'TO0003';

S1 = load(fullfile(dp,[animal '_tifFileList.mat']),'tifFileList');
S2 = load(fullfile(dp,[animal '_stimGroupIDX.mat']),'stimGroupIDX');
S3 = load(fullfile(dp,[animal '_moCorrROI_all.mat']),'moCorROI');

tifFileListStimBPN = S1.tifFileList.stim(S2.stimGroupIDX.BPNStimIDX.tifFileList);
moCorROI = S3.moCorROI;
fprintf('Inputs: %d BPN tifs, %d ROIs\n', numel(tifFileListStimBPN), numel(moCorROI));

%% Step 1: stimParams2TifTable (per-tif stim params from *_Pulses.mat)
fprintf('\n========== Step 1: stimParams2TifTable ==========\n');
tifStimParamTable = stimParams2TifTable(tifFileListStimBPN, dp);
fprintf('size: %s\n', mat2str(size(tifStimParamTable)));
fprintf('vars: %s\n', strjoin(tifStimParamTable.Properties.VariableNames,', '));
if ismember('totalPulses', tifStimParamTable.Properties.VariableNames)
    fprintf('totalPulses per tif: %s\n', mat2str(tifStimParamTable.totalPulses'));
end

%% Step 2: anmlROIbyStimTable (raw anmlROIbyStim before processBPN2P)
fprintf('\n========== Step 2: anmlROIbyStimTable ==========\n');
[T_raw, stimTable] = anmlROIbyStimTable(animal, tifFileListStimBPN, moCorROI, tifStimParamTable);
fprintf('T_raw size: %s\n', mat2str(size(T_raw)));
fprintf('T_raw vars: %s\n', strjoin(T_raw.Properties.VariableNames,', '));
fprintf('stimTable size: %s\n', mat2str(size(stimTable)));
fprintf('stimTable vars: %s\n', strjoin(stimTable.Properties.VariableNames,', '));

% Detect sOnset column (case-insensitively for trace purposes)
onsetVar = T_raw.Properties.VariableNames(endsWith(T_raw.Properties.VariableNames,'sOnset','IgnoreCase',true));
fprintf('sOnset-like cols (case-insensitive): %s\n', strjoin(onsetVar,','));
ov = onsetVar{1};
fprintf('Unique %s in T_raw: %s\n', ov, mat2str(unique(T_raw.(ov))'));
if ismember('BPNdBAmpl', T_raw.Properties.VariableNames)
    fprintf('Unique BPNdBAmpl: %s\n', mat2str(unique(T_raw.BPNdBAmpl)'));
end
sz = cell2mat(cellfun(@size, T_raw.SCALEDfissaFroi, 'uni', 0));
fprintf('SCALEDfissaFroi rows range: [%d %d], cols range: [%d %d]\n', ...
    min(sz(:,1)), max(sz(:,1)), min(sz(:,2)), max(sz(:,2)));

%% Step 3: how many groups would combineDiffOnset see?
nonCell = T_raw.Properties.VariableNames(~varfun(@iscell,T_raw,'OutputFormat','uniform'));
groupVars = setdiff(nonCell, {ov,'Pulse','stimID'}, 'stable');
fprintf('\ncombineDiffOnset groupVars: %s\n', strjoin(groupVars,', '));
[G,~] = findgroups(T_raw(:, groupVars));
onsetsPerGroup = splitapply(@(s){unique(s)}, T_raw.(ov), G);
nMult = sum(cellfun(@(x) numel(x)>1, onsetsPerGroup));
fprintf('Groups total=%d, with multiple onsets=%d\n', numel(onsetsPerGroup), nMult);

%% Step 4: case-sensitivity check of endsWith (used by combineDiffOnset to find onset col)
fprintf('\n========== case-sensitivity of endsWith ==========\n');
fprintf('endsWith(''sonset'',''sOnset'') = %d   (case-sensitive)\n', endsWith('sonset','sOnset'));
fprintf('endsWith(''BPNsOnset'',''sOnset'') = %d\n', endsWith('BPNsOnset','sOnset'));
fprintf('endsWith(''sonset'',''sOnset'',''IgnoreCase'',true) = %d\n', endsWith('sonset','sOnset','IgnoreCase',true));

%% Step 5: run combineDiffOnset on raw table (no t_dFF / dFF columns yet)
fprintf('\n========== Step 5: combineDiffOnset on raw T_raw ==========\n');
try
    out = combineDiffOnset(T_raw);
    fprintf('OK; out height: %d\n', height(out));
    fprintf('out SCALEDfissaFroi class: %s, example dims = %s\n', ...
        class(out.SCALEDfissaFroi), mat2str(size(out.SCALEDfissaFroi{1})));
catch ME
    fprintf('ERROR: %s\n  -> %s\n', ME.identifier, ME.message);
end

% Save the rebuilt raw table for follow-up multi-onset synthesis tests
outFile = fullfile(tempdir,'TO0003_raw_anmlROIbyStim.mat');
save(outFile,'T_raw','stimTable','tifStimParamTable','-v7.3');
fprintf('\nSaved raw outputs to %s\n', outFile);
