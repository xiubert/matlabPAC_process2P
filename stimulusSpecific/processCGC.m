%%
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
%% 

tBaseDRC = [-1.2 0];


%dFF re DRC

anmlROIbyStim.t_dFF_DRC = rowfun(@(t) ...
    {t(find((t>=tBaseDRC(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF_DRC = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBaseDRC(1)),1,'first')...
    find((t<=tBaseDRC(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF_DRC_avg = rowfun(@(F) ...
    {mean(F)},...
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
            xline(2,'--')
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
tBasePT = [1 2];
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
anmlROIbyStim.dFF_PT_preDRCf0 = rowfun(@(dFF_DRC,t) ...
    {dFF_DRC  - ...
    nanmean(dFF_DRC(:,t>=tBasePT(1) & t<=tBasePT(2)),2)},...
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
            xline(2,'--')
            hold off;
            title('ROI'+string(curROIno));
            legend(label(1), label(2),'pure tone')
        end
        clear curROIno
    end
    clear roiSubPlotN
end

%% Peak dFF
pkPTframeBin = 4;
pkPTsigSD=2;

tmp = rowfun(@(dFF,t,PTonset,fr) ...
    pkFcalc(dFF,...
    find(t>=(PTonset+(1/fr)),1,'first'),...
    pkPTframeBin,pkPTsigSD),...
    anmlROIbyStim,'InputVariables',{'dFF_PT_preDRCf0','t_dFF_DRC','PTsOnset','frameRate'},...
    'ExtractCellContents',true,'OutputFormat','cell','OutputVariableNames',{'sigPk','sig','pk'});
anmlROIbyStim.pkPT_sig = tmp(:,1);
anmlROIbyStim.sigPk = tmp(:,2);
anmlROIbyStim.pkPT = tmp(:,3);
%% Save
save(fullfile(dataPath,[animal '_anmlROI_CGCstimTable.mat']),"anmlROIbyStim",'-append');
%% Plot avg of all ROIs

colors.lohi = [0,0.451000000000000,0.741200000000000;0.851000000000000,0.329400000000000,0.102000000000000];
colors.lohiTrace = [0.729400000000000,0.874500000000000,1;1,0.694100000000000,0.541200000000000];
[groups,idC] = findgroups(anmlROIbyStim.dBdelta);
tmp = splitapply(@(x) {vertcat(x{:})},anmlROIbyStim.dFF_PT_avg,groups);
[G0,idC0] = findgroups(idC);
dFF_PT_mean = splitapply(@(x) {cellfun(@nanmean,x,'uni',0)},tmp,G0);
dFF_PT_sem = splitapply(@(x) {cellfun(@SEMcalc,x,'uni',0)},tmp,G0);

figure;
hold on
for r=1:ndBdelta
    x=anmlROIbyStim(1,:).t_dFF_PT{1,1};
    y=dFF_PT_mean{r,1}{1,1};
    yerr=dFF_PT_sem{r,1}{1,1};
    fillSEMplot(x,y,yerr,colors.lohi(r,:),colors.lohiTrace(r,:));
end


xlabel('time/s')
ylabel('dF/F')
xline(2,'--','pure tone')
xlim([1 5])
hold off;
title('Average across cell');
legend('Low contrast', 'High contrast')
%% Low vs High per ROI

x = nan(nCell,1);
y = nan(nCell,1);
for i = 1:nCell
    roi=roiList(i);
    rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
    x(i)=cell2mat(rows.pkPT(rows.dBdelta==dBdeltaList(2)));
    y(i)=cell2mat(rows.pkPT(rows.dBdelta==dBdeltaList(1)));
end

% keep only roiIDs with both values positive
valid = (x>0) & (y>0);
x = x(valid); y = y(valid);
roiList_pos = roiList(valid);

% make scatter
figure;
scatter(x,y,45,'filled','MarkerFaceAlpha',0.8);
hold on;
% identity line
lims = [0 1];
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
    vals = cell2mat(anmlROIbyStim.pkPT(groups == k));
    vals=vals(valid);
    group{k}=vals;
    pkResp_means(k) = mean(vals,'omitnan');
    pkResp_sems(k)= std(vals)/sqrt(nCell);
end

% Create bar plot and error bars (SEM)
figure;
b = bar(1:ndBdelta, pkResp_means,'FaceColor','flat');
hold on;
errorbar(1:ndBdelta, pkResp_means, pkResp_sems, 'k.', 'LineWidth',1);

% Overlay individual points with jitter
rng(0); % for reproducible jitter
jitterAmount = 0.08; % tweak for spread
for k = 1:ndBdelta
    b.CData(k,:)=colors.lohi(k,:);
    vals = cell2mat(anmlROIbyStim.pkPT(groups == k));
    vals=vals(valid);
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
