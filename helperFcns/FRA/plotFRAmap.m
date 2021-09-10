function [] = plotFRAmap(FRA,varargin)
% plotFRAmap: plot FRA output for an animal.
%   [] = plotFRAmap(FRA,varargin)
%
%       INPUT:
%           FRA, --> structure with FRAmap, BFuDB, etc for a given animal
%
%           varargin, --> 'meanAcrossROI': logical --> whether to plot mean
%                                                       responses across all ROI
%                         'plotAllROI':  logical --> whether to plot all
%                                                       ROI FRA
%                         'sigResp':  logical --> whether to plot sig resp
%                         'ROI':  scalar --> ROI id to plot
%                         'plotOctDiff':  logical --> whether to plot FRA
%                                                       in terms of oct diff
%
%   See also FRAmap.

p = inputParser;
addRequired(p,'FRA',@isstruct)
addParameter(p,'meanAcrossROI',true,@islogical)
addParameter(p,'plotAllROI',false,@islogical)
addParameter(p,'sigResp',true,@islogical)
addParameter(p,'ROI',0,@isscalar)
addParameter(p,'plotOctDiff',false,@islogical)

parse(p,FRA,varargin{:});
FRA = p.Results.FRA;
meanAcrossROI = p.Results.meanAcrossROI;
plotAllROI = p.Results.plotAllROI;
sigResp = p.Results.sigResp;
plotROIid = p.Results.ROI;
plotOctDiff = p.Results.plotOctDiff;

if sigResp==1
    uVar = 'uSigPkResp';
    roiVar = 'CellSigPkLinDBfreq';
else
    uVar = 'uPkResp';
    roiVar = 'CellPkRespLinDBfreq';
end
if plotROIid~=0
    meanAcrossROI = false;
end

if meanAcrossROI==1
    
    nCell = size(FRA.BFuDB,1);
    
    FRAfigs.(uVar).FH = figure('Name',sprintf('peak responses across cells | %s',uVar));
    imagesc(flipud(FRA.(uVar)))
    ax = gca;
    xticks(ax,find(rem(FRA.freqList,1000)==0))
    xticklabels(ax,...
        string(FRA.freqList(rem(FRA.freqList,1000)==0)./1000))
    yticks(ax,1:length(FRA.dBlist))
    yticklabels(ax,string(flipud(FRA.dBlist)))
    xlabel('Frequency (Hz)')
    ylabel('dB SPL')
    
    ax.FontSize = 12;
    c1 = colorbar(ax);
    c1.Label.String = '\DeltaF/F';
    c1.FontSize = 12;
    c1.Label.FontSize = 14;
    title(uVar)
    for fqNo = 1:length(FRA.freqList)
        text(fqNo,1.5,num2str(FRA.freqList(fqNo)),...
            'Rotation',90,'Color',[0.65 0.65 0.65],'FontWeight','bold')
    end
end
%% individual cell RF plots

if plotROIid ~= 0
    
    figure('Name',sprintf('ROI: %d | peak dF/F and BFuDB',plotROIid))
    imagesc(flipud(reshape(FRA.(roiVar)(plotROIid,:).*100,size(FRA.uPkResp))));
    if plotOctDiff==1
        octDiffROIbf = octBWfreq(FRA.BFuDB(plotROIid),FRA.freqList);
        evenOctIDX = find(rem(round(abs(octDiffROIbf),3),0.5)==0);
        xticks(evenOctIDX)
        xticklabels(string(round(octDiffROIbf(evenOctIDX),1)))
        xlabel('Octaves from BF')
    else
        xticks(find(rem(FRA.freqList,1000)==0))
        xticklabels(string(FRA.freqList(rem(FRA.freqList,1000)==0)./1000))
        xlabel('Frequency (Hz)')
    end
    yticks(1:length(FRA.dBlist))
    yticklabels(string(flipud(FRA.dBlist)))
    set(gca,'FontSize',16);
    ylabel('dB SPL')
    c = colorbar;
    c.Label.String = '% \DeltaF/F';
    c.FontSize = 16;
    c.Label.FontSize = 16;
    c.Box = 'off';
    c.LineWidth = 2;
    modPlotForPaper(false);
    title(sprintf('ROI: %d | BFuDB: %d Hz',plotROIid,FRA.BFuDB(plotROIid)),...
        'FontSize',10,'FontWeight','normal')
end

if plotAllROI==1
    ROIperFig = 9;
    remROIplotNo = rem(nCell,ROIperFig);
    roiFigNo = floor(nCell/ROIperFig)+(remROIplotNo>=1);
    
    %initialize subplots for multiple ROI per fig
    for roiFigN = 1:roiFigNo
        FRAfigs.ROIpkDelFoFcPlot(roiFigN).FH = figure('Name',...
            'peak dF/F responses and BFuDB for each ROI','Position',...
            [312 73 1304 934]);
        
        for roiSubPlotN = 1:ROIperFig
            curROIno = roiFigN*ROIperFig - ROIperFig + roiSubPlotN;
            if curROIno <= nCell
                figure(FRAfigs.ROIpkDelFoFcPlot(roiFigN).FH)
                FRAfigs.ROIpkDelFoFcPlot(roiFigN).(['ax'...
                    num2str(roiSubPlotN)]) = subplot(3,3,roiSubPlotN);
            end
            clear curROIno
        end
        clear roiSubPlotN
    end
    
    %plot RFs
    for roiN = 1:nCell
        roiFigN = ceil((roiN/ROIperFig));
        roiSubPlotN = rem(roiN,ROIperFig) +...
            (rem(roiN,ROIperFig)<1)*ROIperFig;
        ax = FRAfigs.ROIpkDelFoFcPlot(roiFigN).(['ax'...
            num2str(roiSubPlotN)]);
        
        imagesc(ax,...
            flipud(reshape(FRA.(roiVar)(roiN,:),size(FRA.(uVar)))));
        xticks(ax,find(rem(FRA.freqList,1000)==0))
        xticklabels(ax,...
            string(FRA.freqList(rem(FRA.freqList,1000)==0)./1000))
        yticks(ax,1:length(FRA.dBlist))
        yticklabels(ax,string(flipud(FRA.dBlist)))
        ax.FontSize = 16;
        xlabel(ax,'Frequency (Hz)')
        ylabel(ax,'dB SPL')
        c = colorbar(ax);
        c.Label.String = '\DeltaF/F';
        c.FontSize = 16;
        c.Label.FontSize = 16;
        title(ax,sprintf('ROI: %d | BFuDB: %d Hz',roiN,FRA.BFuDB(roiN)),...
            'FontSize',10,'FontWeight','normal')
        
        clear roiFigN roiSubPlotN xticklabelsH xticklabels_new yticklabelsH yticklabels_new
    end
end