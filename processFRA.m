%% Clear Workspace first if continuing AFTER GETTING tifFileList.map.SCALEDfissaFroi

if ~exist('dataPath','var')
    dataPath='C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\TO0003';
    animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
    load(fullfile(dataPath,[animal '_tifFileList.mat']))
end
if ~exist('nCell','var')
    nCell=size(tifFileList.stim(1).rawFroi,1);
end
if exist('ROIoutputTables','var')
    tifStimParamTable=ROIoutputTables{2};
    anmlROIbyStim=ROIoutputTables{4};
    stimTable=ROIoutputTables{6};
end
%% FRAmap
%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT IF NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkPTsigSD = 2;
nFramesPostPulse = 2;
plotAllROI = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MapOutputDir = fullfile(dataPath,'BFMap');
% if ~isfolder(MapOutputDir)
%     mkdir(MapOutputDir)
% end 

FRAmap = FRAmap(tifFileList,pkPTsigSD,nFramesPostPulse,'SCALEDfissaFroi');
save(fullfile(dataPath,[animal '_FRAmap.mat']),'dataPath','FRAmap','-v7.3')
%% PLOT OUTPUT | Sig Responses

plotFRAmap(FRAmap,'plotAllROI',plotAllROI)

figure;
semilogx(FRAmap.freqList,nanmean(FRAmap.uSigPkResp,1));
[~, maxID] = max(nanmean(FRAmap.uSigPkResp,1));

disp(FRAmap.freqList(maxID))

%% PLOT OUTPUT | All Responses

plotFRAmap(FRAmap,'sigResp',false,'plotAllROI',plotAllROI)

figure;
semilogx(FRAmap.freqList,nanmean(FRAmap.uPkResp,1));
xlabel('Frequency/Hz')
ylabel('dF/F')
[~, maxID] = max(nanmean(FRAmap.uPkResp,1));
disp(FRAmap.freqList(maxID))