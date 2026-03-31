%% RESP to BPN; one sound per tif
% calculate average F
anmlROIbyStim.t_total = rowfun(@(F,fr,trigDelay) ...
    {(1:size(F,2))/fr-trigDelay},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi','frameRate','trigDelay'},...
    'ExtractCellContents',true,'OutputFormat','uniform');
anmlROIbyStim.F_avg = rowfun(@(F) ...
    {mean(F)},...
    anmlROIbyStim,'InputVariables',{'SCALEDfissaFroi'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

%% Plot rawF
anmlROIbyStim.roiID = string(strtrim(cellstr(anmlROIbyStim.roiID)));
roiList = unique(anmlROIbyStim.roiID, 'stable');
ROIperFig = 9;
remROIplotNo = rem(nCell,ROIperFig);
roiFigNo = floor(nCell/ROIperFig)+(remROIplotNo>=1);
%initialize subplots for multiple ROI per fig
for roiFigN = 1:roiFigNo
    figure;
    title('raw F responses for each ROI')
    for roiSubPlotN = 1:ROIperFig
        curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
        if curROIno <= nCell
            subplot(3,3,roiSubPlotN);
            roi=roiList(curROIno);
            rows=anmlROIbyStim(anmlROIbyStim.roiID==roi,:);
            hold on;
            label=strings(height(rows));
            for r=1:height(rows)
                x=rows.t_total{r};
                y=rows.F_avg{r};
                plot(x,y);
                label(r)=string(rows.BPNdBampl(r))+'dB';
            end
            xlabel('time/s')
            ylabel('raw F')
            xline(2,'--','BPN')
            hold off;
            title('ROI'+string(curROIno));
            legend(label)
        end
        clear curROIno
    end
    clear roiSubPlotN
end
%% 

tBase = [1 2];

%dFF

anmlROIbyStim.t_dFF = rowfun(@(t) ...
    {t(find((t>=tBase(1)),1,'first'):end)},...
    anmlROIbyStim,'InputVariables',{'t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

anmlROIbyStim.dFF = rowfun(@(F,t) ...
    {dFoFcalc(F,[find((t>=tBase(1)),1,'first')...
    find((t<=tBase(2)),1,'last')],1)},...
    anmlROIbyStim,'InputVariables',{'F_avg','t_total'},...
    'ExtractCellContents',true,'OutputFormat','uniform');

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
                x=rows.t_dFF{r};
                y=rows.dFF{r};
                plot(x,y);
                label(r)=string(rows.BPNdBampl(r))+'dB';
            end
            xlabel('time/s')
            ylabel('dF/F')
            xline(2,'--','BPN')
            hold off;
            title('ROI'+string(curROIno));
            legend(label)
        end
        clear curROIno
    end
    clear roiSubPlotN
end

%% Plot avg of all ROIs

[groups,idC] = findgroups(anmlROIbyStim.BPNdBampl);
tmp = splitapply(@(x) {vertcat(x{:})},anmlROIbyStim.dFF,groups);
[G0,idC0] = findgroups(idC);
dFF_avg = splitapply(@(x) {cellfun(@nanmean,x,'uni',0)},tmp,G0);
% dFF_sem = splitapply(@(x) {cellfun(@SEMcalc,x,'uni',0)},tmp,G0);

figure;
hold on
for r=1:height(stimTable)
    x=anmlROIbyStim(1,:).t_dFF;
    y=dFF_avg{r,1};
    yerr=dFF_sem{r,1};
    % fillSEMplot(x,y,yerr);
    plot(x,y);
    label(r)=string(idC(r))+'dB';
end


xlabel('time/s')
ylabel('dF/F')
xline(2,'--','BPN')
xlim([1 5])
hold off;
title('Average across cell');
legend(label)
%% Peak dFF
pkPTframeBin = 4;
pkPTsigSD=1.5;

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
save(fullfile(dataPath,[animal '_anmlROI_stimTable.mat']),"anmlROIbyStim",'-append');