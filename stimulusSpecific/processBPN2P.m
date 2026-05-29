% processBPN2P  Compute trial-aligned dF/F and peak responses for one
%               animal's BPN (band-pass noise) session.
%
% Workflow:
%   1. Load <animal>_anmlROI_BPNstimTable_raw.mat (produced by
%      stimParam2ROI). Resolves dataPath via uigetdir if not in workspace.
%   2. Compute per-row t_total, t_dFF, dFF using a configurable
%      pre-onset baseline (baselineSec). Each row's dFF starts at
%      "baselineSec before that row's own BPNsOnset", so the stim onset
%      lives at frame ~ frameRate*baselineSec in dFF coords for every
%      onset value.
%   3. Merge same-stim, different-onset rows via
%      combineDiffOnset(baselineSec). The merge onset-aligns the raw F
%      cells too, so re-running the script on the saved output stays
%      safe.
%   4. Trial-average dF/F per row -> dFF_avg.
%   5. pkFcalc on dFF_avg, with frameStart derived from t_dFF + BPNsOnset
%      (independent of baselineSec) and nFrameWindow rounded up so
%      sub-frame pulse lengths don't truncate the search window.
%   6. Save the processed bundle to <animal>_anmlROI_BPNstimTable.mat
%      (without the _raw suffix). Re-running this script will overwrite
%      that file but never the _raw input.
%
% Output table adds columns: t_total, t_dFF, dFF, dFF_avg, sigPeak,
% sig, pkResp.

if ~exist('dataPath','var')
    dataPath = uigetdir(pwd,'Select animal data folder');
    if isequal(dataPath,0)
        error('No data folder selected.');
    end
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
    load(fullfile(dataPath,[animal '_anmlROI_BPNstimTable_raw.mat']))
end
if ~exist('animal','var')
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
end

%% Parameters (EDIT IF NEEDED)
baselineSec      = 1;   % pre-onset baseline (s) for dFF + onset alignment
pkBPNsigSD       = 2;
nFramesPostPulse = 2;

%% Per-row time vectors and onset-normalized dF/F
anmlROIbyStim.t_total = rowfun(@(F,fr) {(1:size(F,2))/fr}, ...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','frameRate'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.t_dFF = rowfun(@(t,onset) ...
    {t(find((t>onset-baselineSec),1,'first'):end)}, ...
    anmlROIbyStim,'InputVariables',{'t_total','BPNsOnset'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF = rowfun(@(F,t,onset) ...
    {dFoFcalc(F,[find((t>onset-baselineSec),1,'first') ...
                 find((t<=onset),1,'last')],1)}, ...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','t_total','BPNsOnset'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% Merge same-sound, different-onset rows (also aligns raw F cells)
anmlROIbyStim = combineDiffOnset(anmlROIbyStim, baselineSec);

%% Trial-average dF/F per row
anmlROIbyStim.dFF_avg = rowfun(@(F) {mean(F,1,'omitnan')}, ...
    anmlROIbyStim,'InputVariables',{'dFF'}, ...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% Peak detection over the stim window
resultsTable = rowfun(@(F, t_dFF, frameRate, BPNsOnset, BPNmsPulseLen) ...
    pkFcalc(F, ...
        find(t_dFF >= BPNsOnset, 1, 'first'), ...
        ceil(BPNmsPulseLen/1000*frameRate) + nFramesPostPulse, ...
        pkBPNsigSD), ...
    anmlROIbyStim, ...
    'InputVariables',{'dFF_avg','t_dFF','frameRate','BPNsOnset','BPNmsPulseLen'}, ...
    'ExtractCellContents',true,'NumOutputs',5,'OutputFormat','cell');

anmlROIbyStim.sigPeak = resultsTable(:,1);
anmlROIbyStim.sig     = resultsTable(:,2);
anmlROIbyStim.pkResp  = resultsTable(:,3);

%% Save processed bundle (fresh save; never touches the _raw input)
save(fullfile(dataPath,[animal '_anmlROI_BPNstimTable.mat']), ...
    'anmlROIbyStim','stimTable','tifStimParamTable','dataPath','-v7.3');

%% --- Plot 1: single ROI, single dB, all reps + mean ---
targetROI = '1';
targetdB  = 30;

rowMask = (anmlROIbyStim.roiID == targetROI) & (anmlROIbyStim.BPNdBAmpl == targetdB);
if ~any(rowMask)
    error('No rows found for roiID=%s and dB=%g', targetROI, targetdB)
end

r = find(rowMask,1);

% t_dFF{r} is a 1 x nFrames row vector; dFF{r} is nReps x nFrames; dFF_avg{r} is 1 x nFrames.
timeMat = anmlROIbyStim.t_dFF{r};
dffMat  = anmlROIbyStim.dFF{r};
dffMean = anmlROIbyStim.dFF_avg{r};
stimOnSec  = anmlROIbyStim.BPNsOnset(r);
stimOffSec = stimOnSec + anmlROIbyStim.BPNmsPulseLen(r)/1000;

figure;
hold on;
plot(timeMat, dffMat, 'Color', [0.7 0.7 0.9]);
plot(timeMat, dffMean, '-k', 'LineWidth', 2);
xline(stimOnSec,'--','BPN');
xline(stimOffSec,'--');
hold off;
xlabel('Time (s)');
ylabel('dF/F');
title(sprintf('roiID=%s   dB=%g   %d repetitions', targetROI, targetdB, size(dffMat,1)));
legend({'Individual measurement','Average'});
grid on;

%% --- Plot 2: single ROI, all dB, mean +/- SEM per dB ---
targetROI = '2';

rows = find(anmlROIbyStim.roiID == targetROI);
if isempty(rows)
    error('No rows found for roiID=%s', targetROI)
end

dbVals = unique(anmlROIbyStim.BPNdBAmpl(rows));

figure; hold on;
cmap = jet(numel(dbVals));
colors = cmap;

for k = 1:numel(dbVals)
    db = dbVals(k);
    mask = (anmlROIbyStim.roiID == targetROI) & (anmlROIbyStim.BPNdBAmpl == db);
    r = find(mask);
    if isempty(r)
        continue
    end

    timeMat = anmlROIbyStim.t_dFF{r};
    dffMat  = anmlROIbyStim.dFF{r};

    mu  = mean(dffMat, 1, 'omitnan');
    sem = std(dffMat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(dffMat),1));

    [hPatch,hPlot] = fillSEMplot(timeMat,mu,sem,colors(k,:),colors(k,:));
    set(hPatch, 'FaceAlpha', 0.25)
    set(hPlot, 'LineWidth', 1.8)
    hLine(k) = plot(NaN, NaN, 'Color', colors(k,:), 'LineWidth', 1.8); %#ok<SAGROW>
end

% All rows for this ROI share the same BPNsOnset / BPNmsPulseLen post-merge.
stimOnSec  = anmlROIbyStim.BPNsOnset(rows(1));
stimOffSec = stimOnSec + anmlROIbyStim.BPNmsPulseLen(rows(1))/1000;
xline(stimOnSec,'--','BPN');
xline(stimOffSec,'--');
xlabel('Time (s)');
ylabel('dF/F');
title(sprintf('ROI %s: average dF/F per dB', targetROI));
legend(hLine, arrayfun(@(v) sprintf('%d dB', v), dbVals, 'UniformOutput', false), 'Location', 'best');
grid on;
hold off;

%% --- Plot 3: population, all ROIs, pooled-trial SEM per dB ---
dbVals = unique(anmlROIbyStim.BPNdBAmpl);

figure; hold on;
cmap = jet(numel(dbVals));
colors = cmap;
hLine = gobjects(numel(dbVals), 1);

for k = 1:numel(dbVals)
    db = dbVals(k);
    mask = (anmlROIbyStim.BPNdBAmpl == db);
    rows = find(mask);
    if isempty(rows)
        continue
    end

    timeMat = anmlROIbyStim.t_dFF{rows(1)};

    allDffCell = anmlROIbyStim.dFF(rows);
    dffMat = vertcat(allDffCell{:});

    mu  = mean(dffMat, 1, 'omitnan');
    sem = std(dffMat, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(dffMat), 1));

    [hPatch,hPlot] = fillSEMplot(timeMat,mu,sem,colors(k,:),colors(k,:));
    set(hPatch, 'FaceAlpha', 0.25)
    set(hPlot, 'LineWidth', 1.8)
    hLine(k) = plot(NaN, NaN, 'Color', colors(k,:), 'LineWidth', 1.8);
end

% All rows share BPNsOnset/BPNmsPulseLen after combineDiffOnset's min(onset)
% normalization; take from any row.
stimOnSec  = anmlROIbyStim.BPNsOnset(1);
stimOffSec = stimOnSec + anmlROIbyStim.BPNmsPulseLen(1)/1000;
xline(stimOnSec, '--', 'BPN');
xline(stimOffSec, '--');
xlabel('Time (s)');
ylabel('dF/F');
title('All ROIs: Population Average dF/F per dB');
legend(hLine, arrayfun(@(v) sprintf('%d dB', v), dbVals, 'UniformOutput', false), 'Location', 'best');
grid on;
hold off;

%% --- Plot 4: peak dF/F vs dB across cells ---
[G, dBvals] = findgroups(anmlROIbyStim.BPNdBAmpl);

raw_sigPeak = anmlROIbyStim.pkResp;

if iscell(raw_sigPeak)
    data = cellfun(@(x) double(reshape(x,[],1)), raw_sigPeak, 'UniformOutput', false);
    data = cell2mat(cellfun(@(x) [x; NaN(1-numel(x),1)], data, 'UniformOutput', false));
else
    data = double(raw_sigPeak);
end

mu = splitapply(@(x) mean(x,'omitnan'), data, G);
n  = splitapply(@(x) sum(~isnan(x)), data, G);
sd = splitapply(@(x) std(x,'omitnan'),  data, G);
sem = sd ./ sqrt(n);

[dBvals_sorted, idx] = sort(dBvals);
mu = mu(idx);
sem = sem(idx);

figure;
hb = bar(dBvals_sorted, mu, 'FaceColor',[0.8 0.8 0.8]);
hold on;
he = errorbar(dBvals_sorted, mu, sem, 'k', 'LineStyle', 'none', 'LineWidth',1);

pointsByGroup = splitapply(@(x) {x}, data, G);
pointsByGroup = pointsByGroup(idx);

markerColor = [0 0.2 0.6];
for k = 1:numel(dBvals_sorted)
    x0 = dBvals_sorted(k);
    pts = pointsByGroup{k}(:);
    mcount = numel(pts);
    if mcount == 0
        continue
    end
    if numel(dBvals_sorted) > 1
        dx = 0.15 * min(diff(dBvals_sorted));
    else
        dx = 0.15;
    end
    jitter = (rand(mcount,1)*2-1) * dx;
    xs = x0 + jitter;
    scatter(xs, pts, 36, 'MarkerEdgeColor', markerColor, ...
        'MarkerFaceColor', markerColor, 'MarkerFaceAlpha', 0.6);
end

xlabel('Sound level (dB)');
ylabel('peak dF/F');
box on;
hold off;
