% PROCESSCGC  CGC / pure-tone-in-contrast (contrast gain control) analysis.
%
%   Loads <animal>_anmlROI_CGCstimTable.mat, computes dF/F referenced to the
%   DRC baseline and then to the pre-pure-tone (PT) baseline, detects peak PT
%   responses + significance, and plots per-ROI and population summaries.
%
%   METHOD (matches matlabPAC_CGCplot/plotDataTable.m and the manuscript):
%     dFF_DRC          = (F - F0_DRC)/F0_DRC, F0_DRC = mean F over [-1.2 0] s (pre DRC onset).
%     dFF_PT  = dFF_DRC - F0_PT, F0_PT = mean of that dF/F over [1 2] s
%                        (an additive 2nd baseline subtraction, NOT a divisive
%                        re-normalization off raw F). Per-trial, as plotDataTable.
%
%   SIGNIFICANCE: per-ROI PT plot, and the population trace all use the
%   CELL-AVERAGE response (dFF_PT_avg = mean over reps), not individual
%   trials. pkFcalc is fed dFF_PT_avg.
%
%   See also dFoFcalc, pkFcalc, getContrastColors, fillSEMplot

%% ---- PARAMETERS ----
% dF/F baseline windows (seconds, re trial start after trigDelay correction)
tBaseDRC   = [-1.2 0];   % F0 window before DRC onset
tBasePT    = [1 2];      % F0_PT window before pure tone (valid for PTsOnset==2)
PTonsetSec = 2;          % pure-tone onset (s); used for plot markers/xlines

% peak-response detection
pkPTframeBin = 4;        % peak-search window length (frames) after PT onset
pkPTsigSD    = 2;        % significance threshold, in baseline-SD multiples

% plotting
ROIperFig       = 9;     % ROI subplots per figure (3x3)
colors          = getContrastColors();  % contrast color scheme (lohiPre/lohiTracePre)
avgTraceXlim    = [1 5]; % x-limits (s) for population average-trace plot
pkScatterLim    = [0 1]; % axis limits for low-vs-high peak dF/F scatter
jitterAmount    = 0.08;  % horizontal jitter for bar-graph scatter overlay
popTraceSigOnly = true;  % population avg trace: true = significant cells only,
                         % false = all cells. Scatter/bar/t-test always use
                         % significant-in-both-contrasts cells.

%% ---- LOAD DATA ----
if ~exist('dataPath','var')
    dataPath = uigetdir('','Select the animal data folder');
    if isequal(dataPath,0)
        error('processCGC:noDataPath','No data folder selected.');
    end
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
    load(fullfile(dataPath,[animal '_anmlROI_CGCstimTable.mat']))
end

%% Setup

% calculate time vector for each row based on frame rate and trigDelay
% 0 s corresponds to DRC stimulus start, time vector starts at -trigDelay (delay to Ephus stimulus trigger)
anmlROIbyStim.t_total = rowfun(@(F,fr,trigDelay) ...
    {(1:size(F,2))/fr-trigDelay},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.roiID = string(strtrim(cellstr(anmlROIbyStim.roiID)));
roiList = unique(anmlROIbyStim.roiID, 'stable');
nCell = numel(roiList);   % number of cells (unique ROIs); not stored in the .mat
remROIplotNo = rem(nCell,ROIperFig);
roiFigNo = floor(nCell/ROIperFig)+(remROIplotNo>=1);
dBdeltaList = unique(stimTable.dBdelta);
ndBdelta=length(dBdeltaList);

% Guard: the PT F0 window (tBasePT) must lie entirely before pure-tone onset,
% otherwise the "baseline" would include the tone response and every dFF_PT /
% dFF_PT / peak would be mis-normalized. tBasePT is hardcoded here
% (the source indexed it by round(PTsOnset)), so verify all rows are
% compatible rather than failing silently.
if tBasePT(2) > min(anmlROIbyStim.PTsOnset)
    error('processCGC:PTbaselineWindow',...
        ['tBasePT(2)=%.3g s is after the earliest PTsOnset (%.3g s); the PT ' ...
         'baseline would overlap the tone. Set tBasePT to a window before ' ...
         'onset (or branch on PTsOnset as in plotDataTable.m).'],...
        tBasePT(2), min(anmlROIbyStim.PTsOnset));
end
% (raw-F per-ROI diagnostic plot moved to the APPENDIX section at end of file)

%% dFF DRC
% dFF re pre DRC
% fluorescence traces (F(t) = ΔF/F ) first calculated as (F − F0_DRC)/ F0_DRC, 
% where F0_DRC is the average cell fluorescence intensity before DRC sound onset across −1.2 to 0 s. 
anmlROIbyStim.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1)),1,'first')...
    find((t<=tBaseDRC(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% cell-average dF/F re DRC across stim repetitions (1 x nFrames per row);
% used only for the per-ROI overview plot below. (The PT-response average that
% feeds significance is dFF_PT_avg, computed in the next section.)
anmlROIbyStim.dFF_DRC_avg = rowfun(@(F) ...
    {mean(F,1,'omitnan')},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

    %% PLOT dFF DRC
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure('Name','dF/F responses for each ROI');
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows),1);
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_DRC_avg{r};
                plot(x,y);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(PTonsetSec,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end

%
%% dFF re PT

% To quantify responses to pure tones that were preceded by 2 s of contrast DRC, 
% we calculated F(t)-F0_PT, where F0_PT is the average of F(t) across 1–2 s for a 2 s DRC duration.
% F(t) is already dF/F re pre-DRC, so this is a second baseline subtraction (re PT baseline) on top of the first (re DRC baseline). 

% PT response = F(t) - F0_PT, where F(t) is the existing pre-DRC dF/F (dFF_DRC)
% and F0_PT is the mean of that dF/F over the PT baseline window [1 2] s. This
% is an additive second baseline subtraction on the existing dF/F (NOT a fresh
% divisive (F-F0_PT)/F0_PT off raw F), and shares the t_dFF_DRC time axis.
% Per-trial (nReps x nFrames), exactly as plotDataTable.m.
anmlROIbyStim.dFF_PT = rowfun(@(dFF_DRC,t) ...
    {dFF_DRC  - ...
    nanmean(dFF_DRC(:,t>=tBasePT(1) & t<=tBasePT(2)),2)},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC','t_dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% CELL-AVERAGE PT response across reps (1 x nFrames). This is the ONLY intended
% deviation from plotDataTable.m: significance (pkFcalc), the per-ROI PT plot,
% and the population trace are all computed from each cell's averaged response,
% NOT from individual trials. Averaging commutes with the linear F0_PT
% subtraction, so this equals the per-trial dFF_PT averaged over reps.
anmlROIbyStim.dFF_PT_avg = rowfun(@(dFF) ...
    {mean(dFF,1,'omitnan')},...
    anmlROIbyStim,'InputVariables',{'dFF_PT'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% PLOT dFF re PT re DRCf0
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure('Name','dF/F responses re PT re DRCf0 for each ROI');
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows),1);
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_PT_avg{r};
                plot(x,y,'LineWidth', 2);
                if rows.dBdelta(r) == dBdeltaList(1)
                    label(r)='Low contrast';
                elseif rows.dBdelta(r) == dBdeltaList(2)
                    label(r)='High contrast';
                end
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(PTonsetSec,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end

%% Peak dFF PT response and significance
% pkFcalc receives the CELL-AVERAGE PT response (dFF_PT_avg, 1 x nFrames
% per row), so significance is thresholded on each cell's averaged trace, NOT on
% individual trials.
%
% SIGNIFICANCE BASELINE: pkFcalc uses frames 1:frameStart of the trace it is
% given as the baseline (mean + pkPTsigSD*SD). To make that baseline the 1-2 s
% pre-PT F0_PT window (sustained DRC just before the tone) we (1) crop the
% trace+time to t>=tBasePT(1), and (2) set frameStart at PT onset (NOT
% PTonset+1/fr). This makes baseline = [tBasePT(1), PT onset] exactly; using
% PTonset+1/fr would push the baseline one frame past onset, pulling the rising
% PT response into the baseline and inflating the SD. Including the DRC-onset
% transient (the whole pre-PT period) would likewise inflate the SD. The peak
% search runs nFrameWindow frames from PT onset; the onset frame sits in both
% baseline and peak window (the usual minor pkFcalc overlap) but precedes the
% Ca rise, so it does not affect the peak value.

tmp = rowfun(@(dFF,t,PTonset) ...
    pkFcalc(dFF(:, t>=tBasePT(1)),...
    find(t(t>=tBasePT(1))>=PTonset,1,'first'),...
    pkPTframeBin,pkPTsigSD),...
    anmlROIbyStim,'InputVariables',{'dFF_PT_avg','t_dFF_DRC','PTsOnset'},...
    'ExtractCellContents',true,'OutputFormat','cell','OutputVariableNames',{'sigPk','sig','pk'});
anmlROIbyStim.pkPT_sig = tmp(:,1);
anmlROIbyStim.sigPk = tmp(:,2);
anmlROIbyStim.pkPT = tmp(:,3);

%% Per-ROI peak + significance
% Per-ROI peak dF/F and significance flag: rows in roiList order, one column
% per contrast level (dBdeltaList order). Built by explicit roiID+dBdelta
% lookup so downstream indexing does NOT depend on the row order of
% anmlROIbyStim. pkByROI holds the peak of each cell-average trace (pkPT);
% sigByROI holds the pkFcalc significance flag for that same average trace.

pkByROI  = nan(nCell,ndBdelta);
sigByROI = false(nCell,ndBdelta);
for i = 1:nCell
    rows = anmlROIbyStim(anmlROIbyStim.roiID==roiList(i),:);
    for k = 1:ndBdelta
        sel = rows.dBdelta==dBdeltaList(k);
        v = cell2mat(rows.pkPT(sel));
        if ~isempty(v)
            pkByROI(i,k)  = v;                          % errors loudly if >1 row per (ROI,contrast)
            sigByROI(i,k) = logical(cell2mat(rows.sigPk(sel)));
        end
    end
end

% SIGNIFICANCE FILTER: keep only cells whose cell-AVERAGE peak is significant
% in BOTH contrasts (all(sigByROI,2)). This is a significance criterion from
% pkFcalc, NOT an amplitude (pkPT>0) threshold, and matches the source
% matlabPAC_CGCplot/plotDataTable.m (sigCellID = all(sig,2)). The scatter,
% bar graph and paired t-test below all use this mask.
valid = all(sigByROI,2);

%% Save output
save(fullfile(dataPath,[animal '_anmlROI_CGCstimTable.mat']),"anmlROIbyStim",'-append');

%% Plot avg of all ROIs
% Restrict to cells with significant responses in both contrasts when popTraceSigOnly is true, so this
% population trace uses the same cell set as the scatter/bar below; otherwise
% average across all ROIs.
if popTraceSigOnly
    Tpop = anmlROIbyStim(ismember(anmlROIbyStim.roiID,roiList(valid)),:);
else
    Tpop = anmlROIbyStim;
end

% per-cell PT response (dFF_PT_avg) -> mean +/- between-cell SEM per
% contrast. Same manuscript quantity that the peak/scatter/bar below quantify.
% group each cell's avg trace by contrast -> [nCells x nFrames] per contrast,
% then mean and between-cell SEM ACROSS CELLS. Dimension 1 is explicit so a
% single-cell group is averaged across cells, not collapsed across time.
groups = findgroups(Tpop.dBdelta);   % ascending dBdelta -> group r == dBdeltaList(r)
dFF_PT_mean = splitapply(@(x) {nanmean(vertcat(x{:}),1)},Tpop.dFF_PT_avg,groups);
dFF_PT_sem  = splitapply(@(x) {SEMcalc(vertcat(x{:}),1)},Tpop.dFF_PT_avg,groups);

figure;
hold on
for r=1:ndBdelta
    x=anmlROIbyStim(1,:).t_dFF_DRC{1,1};
    y=dFF_PT_mean{r,1};
    yerr=dFF_PT_sem{r,1};
    fillSEMplot(x,y,yerr,colors.lohiPre(r,:),colors.lohiTracePre(r,:));
end


xlabel('time/s')
ylabel('dF/F')
xline(PTonsetSec,'--','pure tone')
xlim(avgTraceXlim)
hold off;
title('Average across cell');
legend('Low contrast', 'High contrast')

%% Low vs High per ROI
% pkByROI / valid computed above (Per-ROI peak + significance section).
% valid = significant-in-both-contrasts cells.
x = pkByROI(valid,2);   % high contrast (dBdeltaList(2))
y = pkByROI(valid,1);   % low contrast  (dBdeltaList(1))
roiList_pos = roiList(valid);

% make scatter
figure;
scatter(x,y,45,'filled','MarkerFaceAlpha',0.8);
hold on;
% identity line
lims = pkScatterLim;
plot(lims, lims, '--k', 'LineWidth', 1);
hold off;

xlabel('High contrast');
ylabel('Low contrast');
title('peak dF/F per roi');
% grid on;
axis equal;
xlim(lims); ylim(lims);

%% Bar graph of peak dF/F per contrast, with scatter overlay
group=cell(ndBdelta,1);
pkResp_means=NaN(ndBdelta,1);
pkResp_sems=NaN(ndBdelta,1);

for k = 1:ndBdelta
    % pkByROI and valid share roiList order (see Per-ROI peak + significance)
    vals = pkByROI(valid,k);
    group{k}=vals;
    pkResp_means(k) = mean(vals,'omitnan');
    % SEM over the valid (filtered) cells only; SEMcalc divides by the count
    % of non-NaN entries, not total nCell, so it is not understated.
    pkResp_sems(k)= SEMcalc(vals);
end

% Create bar plot and error bars (SEM)
figure;
b = bar(1:ndBdelta, pkResp_means,'FaceColor','flat');
hold on;
errorbar(1:ndBdelta, pkResp_means, pkResp_sems, 'k.', 'LineWidth',1);

% Overlay individual points with jitter
rng(0); % for reproducible jitter
for k = 1:ndBdelta
    b.CData(k,:)=colors.lohiPre(k,:);
    vals = group{k};
    x = (k) + (rand(size(vals)) - 0.5) * 2 * jitterAmount;
    % Plot points
    scatter(x, vals, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

% Formatting

xticklabels({'Low contrast','High contrast'});             % works for numeric/categorical/string
ylabel('peak dF/F');
% title('pkPT by dBdelta (individual points and group mean)');
box on;

if ndBdelta == 2
    [h,p,ci,stats] = ttest(group{1}, group{2});
    means=[pkResp_means(1),pkResp_means(2)];
    sems=[pkResp_sems(1),pkResp_sems(2)];
    % Add significance star or text
    yMax = max([means + sems]) ;
    yStar = yMax + 0.5*range([means sems]) ;  % vertical position for star/line
    % Draw bar connecting line
    plot([1 2], [yStar yStar], '-k', 'LineWidth',1);
    % Draw short ticks
    plot([1 1], [yStar-0.1*range([means sems]) yStar], '-k', 'LineWidth',1);
    plot([2 2], [yStar-0.1*range([means sems]) yStar], '-k', 'LineWidth',1);

    if p < 0.001
        sigtxt = '***';
    elseif p < 0.01
        sigtxt = '**';
    elseif p < 0.05
        sigtxt = '*';
    else
        sigtxt = sprintf('p=%.3g', p);
    end
    text(1.5, yStar + 0.03*range([means sems]), sigtxt, 'HorizontalAlignment','center', 'FontSize',14);
end

hold off;

%% ========================== APPENDIX ==========================
% Optional QC plots, not part of the standard analysis. Opt-in via toggle so
% a normal script run does not spawn many figures.

%% APPENDIX: raw-F per-ROI diagnostic
% Plots raw (un-normalized) average fluorescence per ROI, ROIperFig per
% figure (3x3), one trace per contrast. Useful for spotting bleaching/motion
% artifacts that dF/F normalization would hide.
runRawFdiagnostic = false;   % set true to generate the raw-F QC figures

if runRawFdiagnostic
    % raw F averaged across stim repetitions (1 x nFrames per row)
    anmlROIbyStim.F_avg = rowfun(@(F) ...
        {mean(F,1,'omitnan')},...
        anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi'},...
        'ExtractCellContents',true,'OutputFormat','uniform');

    for roiFigN = 1:roiFigNo
        figure('Name','raw F responses for each ROI');
        for roiSubPlotN = 1:ROIperFig
            curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
            if curROIno <= nCell
                subplot(3,3,roiSubPlotN);
                roi=roiList(curROIno);
                rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
                hold on;
                label=strings(height(rows),1);
                for r=1:height(rows)
                    x=rows.t_total{r};
                    y=rows.F_avg{r};
                    plot(x,y);
                    if rows.dBdelta(r) == dBdeltaList(1)
                        label(r)='Low contrast';
                    elseif rows.dBdelta(r) == dBdeltaList(2)
                        label(r)='High contrast';
                    end
                end
                xlabel('time/s')
                ylabel('raw F')
                xline(PTonsetSec,'--','pure tone')
                hold off;
                title('ROI'+string(curROIno));
                legend(label)
            end
            clear curROIno
        end
        clear roiSubPlotN
    end
end
