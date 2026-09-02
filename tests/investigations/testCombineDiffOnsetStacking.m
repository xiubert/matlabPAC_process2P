% testCombineDiffOnsetStacking   [INVESTIGATION SCRIPT -- no assertions]
%
% Not part of the test suite: it prints what combineDiffOnset does to the cell
% columns, it does not assert. The behaviour it documents is now covered by
% integration/testProcessBPN2P_e2e.m. Kept for provenance.
%
% Verify how combineDiffOnset.m stacks the cell columns when actual
% merging happens (TO0003 has BPNsOnset = {1, 2}, every group has both,
% so combineDiffOnset merges every pair).
%
% Specific concern: when subTable has multiple rows for a group and
% combineDiffOnset does vertcat(subTable.(varName){:}), does the result
% end up as a single matrix inside the cell, or as a cell-of-cells?
%
% Also checks:
%   - The class and dims of every cell column after merge.
%   - Whether row counts double as expected (5 reps x 2 onsets -> 10).
%   - Whether SCALEDfissaFroi rows from the two onsets are time-aligned
%     (they should NOT be, because anmlROIbyStimTable slices each pulse
%     window per-tif without onset shifting).
%   - Whether dFF rows ARE time-aligned (they should be, because
%     processBPN2P's per-row dFoFcalc strips off "pre-(onset-1s)" first).
%
% Reads the TO0003 BPN _raw table straight off disk. It used to depend on a
% /tmp artifact left behind by traceBPNcombineDiffOnset.m, which made it fail
% unless that script had been run first in the same session -- the _raw file
% stimParam2ROI writes holds the same table, so read that instead.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % tests/ root
cfg = testConfig();

dp = requireFixture(cfg.exampleAnimalDir,'TO0003 animal folder');
raw = requireFixture(fullfile(dp,'TO0003_anmlROI_BPNstimTable_raw.mat'), ...
    'TO0003 BPN _raw table');
T = getfield(load(raw,'anmlROIbyStim'),'anmlROIbyStim'); %#ok<GFLD>

%% Step 1: replicate processBPN2P's t_total/t_dFF/dFF prep on this raw table
fprintf('\n========== Step 1: add t_total/t_dFF/dFF (processBPN2P prep) ==========\n');
T.t_total = rowfun(@(F,fr) {(1:size(F,2))/fr}, ...
    T,'InputVariables',{'SCALEDfissaFroi','frameRate'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');
T.t_dFF = rowfun(@(t,onset) {t(find((t>onset-1),1,'first'):end)}, ...
    T,'InputVariables',{'t_total','BPNsOnset'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');
T.dFF = rowfun(@(F,t,onset) ...
    {dFoFcalc(F,[find((t>onset-1),1,'first') find((t<=onset),1,'last')],1)}, ...
    T,'InputVariables',{'SCALEDfissaFroi','t_total','BPNsOnset'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');

fprintf('Per-row sizes BEFORE combineDiffOnset:\n');
for r = 1:2
    fprintf('  row %d (BPNsOnset=%g): SCALEDfissaFroi=%s, t_dFF=%s, dFF=%s\n', ...
        r, T.BPNsOnset(r), mat2str(size(T.SCALEDfissaFroi{r})), ...
        mat2str(size(T.t_dFF{r})), mat2str(size(T.dFF{r})));
end

%% Step 2: run combineDiffOnset
fprintf('\n========== Step 2: combineDiffOnset ==========\n');
out = combineDiffOnset(T);
fprintf('input rows: %d, output rows: %d (expected %d = 108)\n', ...
    height(T), height(out), height(T)/2);

%% Step 3: shape audit on output cell columns
fprintf('\n========== Step 3: shape audit on output cell columns ==========\n');
cellCols = out.Properties.VariableNames(...
    varfun(@iscell, out, 'OutputFormat','uniform'));
for c = 1:numel(cellCols)
    v = out.(cellCols{c});
    sample = v{1};
    fprintf('  %-18s | cell-content class=%-10s | row1 dims=%-12s | iscell=%d\n', ...
        cellCols{c}, class(sample), mat2str(size(sample)), iscell(sample));
end

%% Step 4: confirm dFF rows are time-aligned across onsets
fprintf('\n========== Step 4: dFF time-alignment across vertcat''d onsets ==========\n');
fr = out.frameRate(1);
sonsetSaved = out.BPNsOnset(1);
fprintf('  saved BPNsOnset (smallest of the merged): %g\n', sonsetSaved);
% In the original T, each per-row dFF starts at "first frame where t>onset-1".
% After processBPN2P trims t_dFF/dFF to minLength, frame 1 of every trial =
% "1s before each trial''s own onset". Stim onset is therefore at index ~fr.
fprintf('  expected stim-onset index in dFF (any onset): %g\n', fr);
% Locate the row that holds vertcat''d dFF for some ROI-dB combo
r = 1;
M = out.dFF{r};
fprintf('  out.dFF{%d} dims: %s\n', r, mat2str(size(M)));
% mean trace, find peak frame
mu = mean(M, 1, 'omitnan');
[~, peakFrame] = max(mu);
fprintf('  peak frame of mean dFF: %d (closer to %g = onset-aligned, OR closer to %g = sonset=2-aligned)\n', ...
    peakFrame, fr, fr*2);

% Per-onset peaks (rows 1:5 originally from sonset=1, rows 6:10 from sonset=2)
mu1 = mean(M(1:5,:), 1, 'omitnan');
mu2 = mean(M(6:10,:), 1, 'omitnan');
[~, p1] = max(mu1);
[~, p2] = max(mu2);
fprintf('  peak frame of sonset=1 sub-rows (1:5): %d\n', p1);
fprintf('  peak frame of sonset=2 sub-rows (6:10): %d\n', p2);

%% Step 5: confirm SCALEDfissaFroi rows are NOT time-aligned across onsets
fprintf('\n========== Step 5: SCALEDfissaFroi alignment check ==========\n');
M2 = out.SCALEDfissaFroi{r};
fprintf('  out.SCALEDfissaFroi{%d} dims: %s\n', r, mat2str(size(M2)));
mu1raw = mean(M2(1:5,:), 1, 'omitnan');
mu2raw = mean(M2(6:10,:), 1, 'omitnan');
[~, q1] = max(mu1raw);
[~, q2] = max(mu2raw);
fprintf('  peak frame of sonset=1 sub-rows raw F: %d (expect near %g = fr*1)\n', q1, fr*1);
fprintf('  peak frame of sonset=2 sub-rows raw F: %d (expect near %g = fr*2)\n', q2, fr*2);
fprintf('  -> if these differ by ~fr frames, SCALEDfissaFroi is NOT onset-aligned post-vertcat.\n');

%% Step 6: confirm column counts on cells before & after combineDiffOnset
fprintf('\n========== Step 6: column-count audit ==========\n');
% Pick a group manually and inspect the input rows that got merged
% Find one (roiID, BPNdBAmpl) group in T
idx = find(T.roiID==T.roiID(1) & T.BPNdBAmpl==T.BPNdBAmpl(1));
fprintf('  matched %d input rows for (%s, %gdB)\n', numel(idx), T.roiID(1), T.BPNdBAmpl(1));
for k = 1:numel(idx)
    fprintf('    row %d (BPNsOnset=%g): SCALEDfissaFroi cols=%d, dFF cols=%d\n', ...
        idx(k), T.BPNsOnset(idx(k)), size(T.SCALEDfissaFroi{idx(k)},2), size(T.dFF{idx(k)},2));
end
