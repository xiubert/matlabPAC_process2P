% PROCESSCGC  CGC / pure-tone-in-contrast (contrast gain control) analysis.
%
%   Loads <animal>_anmlROI_CGCstimTable.mat, computes dF/F referenced to the
%   DRC baseline and then to the pre-pure-tone (PT) baseline, detects peak PT
%   responses + significance, and plots per-ROI and population summaries.
%
%   SIGNIFICANCE: peak significance is computed from the cell-AVERAGE dF/F
%   trace (averaged across stim repetitions), not from individual trial
%   traces.  This is the accurate approach and is what pkFcalc expects (see
%   its docstring).  The average is formed in dFF_DRC_avg and carried through
%   dFF_PT_preDRCf0 into pkFcalc.  Do NOT switch the dFF_PT_preDRCf0 input
%   back to per-trial dFF_DRC.
%
%   Adapted from matlabPAC_CGCplot/plotDataTable.m, which thresholded
%   significance on per-trial traces.
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
%% PLOT RESP to CGC
% calculate average F
anmlROIbyStim.t_total = rowfun(@(F,fr,trigDelay) ...
    {(1:size(F,2))/fr-trigDelay},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
% anmlROIbyStim.F_avg = rowfun(@(F) ...
%     {mean(F)},...
%     anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi'},...
%     'ExtractCellContents',true,'OutputFormat','uniform');

%% Plot rawF
anmlROIbyStim.roiID = string(strtrim(cellstr(anmlROIbyStim.roiID)));
roiList = unique(anmlROIbyStim.roiID, 'stable');
ROIperFig = 9;
remROIplotNo = rem(nCell,ROIperFig);
roiFigNo = floor(nCell/ROIperFig)+(remROIplotNo>=1);
dBdeltaList = unique(stimTable.dBdelta);
ndBdelta=length(dBdeltaList);
% %initialize subplots for multiple ROI per fig
% for roiFigN = 1:roiFigNo
%     figure;
%     title('raw F responses for each ROI')
%     for roiSubPlotN = 1:ROIperFig
%         curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
%         if curROIno <= nCell
%             subplot(3,3,roiSubPlotN);
%             roi=roiList(curROIno);
%             rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
%             hold on;
%             label=strings(height(rows));
%             for r=1:height(rows)
%                 x=rows.t_total{r};
%                 y=rows.F_avg{r};
%                 plot(x,y);
%                 if rows.dBdelta(r) == dBdeltaList(1)
%                     label(r)='Low contrast';
%                 elseif rows.dBdelta(r) == dBdeltaList(2)
%                     label(r)='High contrast';
%                 end
%             end
%             xlabel('time/s')
%             ylabel('raw F')
%             xline(2,'--','pure tone')
%             hold off;
%             title('ROI'+string(curROIno));
%             legend(label(1), label(2))
%         end
%         clear curROIno
%     end
%     clear roiSubPlotN
% end
%% dFF re DRC
anmlROIbyStim.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1)),1,'first')...
    find((t<=tBaseDRC(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

% cell-average dF/F across stim repetitions (1 x nFrames per row).
% This average is what feeds significance detection (see dFF_PT_preDRCf0).
anmlROIbyStim.dFF_DRC_avg = rowfun(@(F) ...
    {mean(F,1,'omitnan')},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
%%
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure;
    title('dF/F responses for each ROI')
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows));
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
%% dFF re PT re DRCf0
anmlROIbyStim.t_dFF_PT = rowfun(@(t) ...
    {t(find((t>=tBasePT(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
anmlROIbyStim.dFF_PT = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBasePT(1)),1,'first')...
    find((t<=tBasePT(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
anmlROIbyStim.dFF_PT_avg = rowfun(@(F) ...
    {mean(F)},...
    anmlROIbyStim,'InputVariables',{'dFF_PT'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
% Re-baseline the cell-AVERAGE trace (dFF_DRC_avg) to the pre-PT window.
% Significance must be computed on this average, not per-trial dFF_DRC, so
% the input variable is dFF_DRC_avg (1 x nFrames). See header note.
anmlROIbyStim.dFF_PT_preDRCf0 = rowfun(@(dFF_DRC_avg,t) ...
    {dFF_DRC_avg  - ...
    nanmean(dFF_DRC_avg(:,t>=tBasePT(1) & t<=tBasePT(2)),2)},...
    anmlROIbyStim,'InputVariables',{'dFF_DRC_avg','t_dFF_DRC'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% 
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure;
    title('dF/F responses re PT re DRCf0 for each ROI');
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows));
            for r=1:height(rows)
                x=rows.t_dFF_DRC{r};
                y=rows.dFF_PT_preDRCf0{r};
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

%% Peak dFF
% pkFcalc receives the cell-average dFF_PT_preDRCf0 (1 x nFrames per row),
% so significance is thresholded on the averaged trace (correct approach).
tmp = rowfun(@(dFF,t,PTonset,fr) ...
    pkFcalc(dFF,...
    find(t>=(PTonset+(1/fr)),1,'first'),...
    pkPTframeBin,pkPTsigSD),...
    anmlROIbyStim,'InputVariables',{'dFF_PT_preDRCf0','t_dFF_DRC','PTsOnset','frameRate'},...
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

%% Save
save(fullfile(dataPath,[animal '_anmlROI_CGCstimTable.mat']),"anmlROIbyStim",'-append');
%% Plot avg of all ROIs
% Restrict to significant-in-both cells when popTraceSigOnly is true, so this
% population trace uses the same cell set as the scatter/bar below; otherwise
% average across all ROIs.
if popTraceSigOnly
    Tpop = anmlROIbyStim(ismember(anmlROIbyStim.roiID,roiList(valid)),:);
else
    Tpop = anmlROIbyStim;
end
[groups,idC] = findgroups(Tpop.dBdelta);
tmp = splitapply(@(x) {vertcat(x{:})},Tpop.dFF_PT_avg,groups);
[G0,idC0] = findgroups(idC);
dFF_PT_mean = splitapply(@(x) {cellfun(@nanmean,x,'uni',0)},tmp,G0);
dFF_PT_sem = splitapply(@(x) {cellfun(@SEMcalc,x,'uni',0)},tmp,G0);

figure;
hold on
for r=1:ndBdelta
    x=anmlROIbyStim(1,:).t_dFF_PT{1,1};
    y=dFF_PT_mean{r,1}{1,1};
    yerr=dFF_PT_sem{r,1}{1,1};
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

%% Bar graph
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
