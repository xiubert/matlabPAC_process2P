%% Gather traces from selected tif files and plot average cell responses
% must first run roiGUI to draw ROI on the first selected tif.

% define baseline idx for dFF
% see sROI.TriggerParams, sROI.Stim, sROI.t
% find(sROI.t==sROI.TriggerParams.stimDelay)
baselineIDX = [10 20];

% select list of tif files to average across
dataPath = uigetdir('D:');
allTifFiles = dir(fullfile(dataPath,'*.tif'));
selIDX = false(length(allTifFiles),1);
seqSelIDX = listdlg('PromptString',{'Select tif files for averaging','(cancel if none)'},...
    'ListString',{allTifFiles.name}','SelectionMode','multiple');
selIDX(seqSelIDX) = true;
clear seqSelIDX

tifFileList = allTifFiles(selIDX);

%use ROIs drawn on first tif for all tifs (for time efficiency)
try
    %write ROIoutput to tif list
    load(fullfile(tifFileList(1).folder,strrep(tifFileList(1).name,'.tif','_roiOutput.mat')))
    tifFileList(1).nFrames = fluo2p.nFrames;
    try
        tifFileList(1).frameRate = fluo2p.frameRate;
        tifFileList(1).rawFroi = fluo2p.rawFroi;
        tifFileList(1).dFF = dFoFcalc(fluo2p.rawFroi,...
            baselineIDX,1);
    catch
        tifFileList(1).frameRate = fluo2p.settings.FrameRate;
        tifFileList(1).rawFroi = fluo2p.rawintensityROI;
        tifFileList(1).dFF = dFoFcalc(fluo2p.rawintensityROI,...
            baselineIDX,1);
    end
    
    for nTif = 2:length(tifFileList)
        [tifFileList(nTif).rawFroi,tifFileList(nTif).nFrames,...
            tifFileList(nTif).frameRate] = TifROImask2rawFroi(...
            fullfile(tifFileList(nTif).folder,tifFileList(nTif).name),...
            fluo2p.roi);
         tifFileList(nTif).dFF = dFoFcalc(tifFileList(nTif).rawFroi,...
             baselineIDX,1);
    end
    
    % save adhoc data
    % save(fullfile(adhocMapOutputDir,'adHocMapTifList.mat'),'tifFileList')
catch
    error('could not find roiOutput for first tif file')
end

%% Plot average responses
% Calculate the standard deviation of the dFF responses for error shading
A = cat(3, tifFileList.dFF);
figure;
plot(mean(A,3)')
hold on
plot(mean(mean(A,3),1),'LineWidth',4)