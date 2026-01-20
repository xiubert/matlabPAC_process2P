%adhoc FRAmap
clearvars; close all; clc;
%use roiGUI to output fluo2p for first tif file
%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT IF NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkPTsigSD = 2;
nFramesPostPulse = 2;
plotAllROI = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataPath = uigetdir('C:');
adhocMapOutputDir = fullfile(dataPath,'adHocMap');

if ~isfolder(adhocMapOutputDir)
    mkdir(adhocMapOutputDir)
    
    try
        animal = regexp(dataPath,'[A-Z]{2}\d{4}','match','once');
        if isempty(animal)
            error('')
        end
    catch
        animal = char(inputdlg('Enter animal ID:','Animal Input',[1 35],{'AA0000'}));
    end
    
    allTifFiles = dir(fullfile(dataPath,[animal '*.tif']));
    FRAmapIDX = false(length(allTifFiles),1);
    FRAmapSelIDX = listdlg('PromptString',{'Select tif files for FRA mapping','(cancel if none)'},...
        'ListString',{allTifFiles.name}','SelectionMode','multiple');
    FRAmapIDX(FRAmapSelIDX) = true;
    clear FRAmapSelIDX
    treatment = repmat({'FRAmap'},sum(FRAmapIDX),1);
    
    tifFileList.map = allTifFiles(FRAmapIDX);
    [tifFileList.map.treatment] = treatment{:};
    
    %use ROIs drawn on first tif for all map tifs (for time efficiency)
    try
        %write ROIoutput to tif list
        load(fullfile(tifFileList.map(1).folder,strrep(tifFileList.map(1).name,'.tif','_roiOutput.mat')))
        tifFileList.map(1).nFrames = fluo2p.nFrames;
        try
            tifFileList.map(1).frameRate = fluo2p.frameRate;
            tifFileList.map(1).rawFroi = fluo2p.rawFroi;
        catch
            tifFileList.map(1).frameRate = fluo2p.settings.FrameRate;
            tifFileList.map(1).rawFroi = fluo2p.rawintensityROI;
        end
        
        for nTif = 2:length(tifFileList.map)
            [tifFileList.map(nTif).rawFroi,tifFileList.map(nTif).nFrames,...
                tifFileList.map(nTif).frameRate] = TifROImask2rawFroi(...
                fullfile(tifFileList.map(nTif).folder,tifFileList.map(nTif).name),...
                fluo2p.roi);
        end
        
        % save adhoc map data
        save(fullfile(adhocMapOutputDir,'adHocMapTifList.mat'),'tifFileList')
    catch
        error('could not find roiOutput for first tif map file')
    end
    
else
    load(fullfile(adhocMapOutputDir,'adHocMapTifList.mat'))
end

adHocFRAmap = FRAmap(tifFileList,pkPTsigSD,...
    nFramesPostPulse,'rawFroi');

%% PLOT OUTPUT | Sig Responses

plotFRAmap(adHocFRAmap,'plotAllROI',plotAllROI)

figure;semilogx(adHocFRAmap.freqList,nanmean(adHocFRAmap.uSigPkResp,1))
[~, maxID] = max(nanmean(adHocFRAmap.uSigPkResp,1));
[~, maxID] = max(nanmean(adHocFRAmap.uPkResp,1));

disp(adHocFRAmap.freqList(maxID))

%% PLOT OUTPUT | All Responses

plotFRAmap(adHocFRAmap,'sigResp',false)

figure;semilogx(adHocFRAmap.freqList,nanmean(adHocFRAmap.uPkResp,1))
[~, maxID] = max(nanmean(adHocFRAmap.uPkResp,1));
