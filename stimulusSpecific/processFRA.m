% processFRA: build a frequency response area (FRA) map for one animal.
%
% Workflow:
%   1. Resolve dataPath (prompts via uigetdir if not already in workspace)
%      and load <animal>_tifFileList.mat. Animal ID is parsed from the
%      folder name (pattern [A-Z]{2}\d{4}).
%   2. If ROIoutputTables exists in the workspace, unpack the per-stim
%      parameter table, anmlROIbyStim, and stimTable from it.
%   3. Call FRAmap on tifFileList using SCALEDfissaFroi traces, with
%      significance threshold pkPTsigSD (SD) and a post-pulse window of
%      nFramesPostPulse frames.
%   4. Save <animal>_FRAmap.mat into dataPath.
%   5. Plot per-ROI significant responses and the cohort-mean
%      significant-response tuning curve; print the peak frequency.
%   6. Plot per-ROI all-response maps and the cohort-mean all-response
%      tuning curve; print the peak frequency.
%
% Inputs are read from the workspace / disk:
%   dataPath        - folder containing <animal>_tifFileList.mat
%   tifFileList     - loaded from disk; supplies rawFroi/SCALEDfissaFroi
%   ROIoutputTables - optional; if present, tables 2/4/6 are unpacked
%
% Editable params (top of "FRAmap" section):
%   pkPTsigSD        - SD threshold for peak vs. pre-stim baseline
%   nFramesPostPulse - frames after stim onset used for the response window
%   plotAllROI       - if true, plotFRAmap renders every ROI tile
%
% Outputs:
%   FRAmap         - struct returned by FRAmap()
%   <animal>_FRAmap.mat saved into dataPath (contains dataPath + FRAmap)
%
% Note: re-running after clearing tifFileList.map.SCALEDfissaFroi requires
% clearing the workspace first so the load branch re-executes.
%
% TODO: compute significant responses from the trial-averaged trace per
% ROI/stim instead of testing each individual trial separately — current
% approach is noisier and biases significance toward high-trial-count
% conditions.

if ~exist('dataPath','var')
    dataPath = uigetdir(pwd,'Select animal data folder');
    if isequal(dataPath,0)
        error('No data folder selected.');
    end
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